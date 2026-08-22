#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-qc and stream-hvf contents have been combined here #

# cell/feat stats ####
# Calculate stats for parquetExprBase-backed expression. Dispatches
# via Giotto's analyzeData(x, param) generic
# 
# Computes, collects, and returns stats in memory
#
#   analyzeData(parquetExprBase, cellStatsParam) 
#     -> per-cell stats data.table
#   analyzeData(parquetExprBase, featStatsParam) 
#     -> per-feature stats data.table

## cellStatsParam ####
#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "cellStatsParam"),
    function(x, param, ...) {
        do.call(.pe_cellstats, c(list(pe = x), as.list(param@param)))
    }
)

## featStatsParam ####
#' @name analyzeData-featStatsParam
#' @rdname analyzeData-featStatsParam
#' @title Streaming per-feature statistics
#' @description
#' [GiottoClass::analyzeData()] method for a [Giotto::featStatsParam-class] on
#' a disk-backed expression store. One streamed pass over the triplet stream,
#' either over every cell or partitioned by a per-cell grouping.
#'
#' The grouped form is reusable beyond QC: per-cluster mean and
#' percent-detected is the input to a dot plot, and the group means are a
#' pseudobulk matrix.
#' @param x a `parquetExprBase` store.
#' @param param a [Giotto::featStatsParam-class].
#' @param groups optional vector of group assignments, one per cell of the
#'   current view, `NA` to exclude a cell. When supplied, the statistics are
#'   taken per (feature, group) instead of over every cell, and the result gains
#'   `group` and `n_cells` columns.
#' @param stats optional character vector of accumulators to compute, any of
#'   `"sum"`, `"sumsq"`, `"nnz"`, `"sum_det"`. Grouped path only. Emitted
#'   columns are whichever the requested accumulators support, so asking for
#'   less genuinely scans for less. Defaults to all four.
#' @param ... additional arguments (none used).
#' @returns A `data.table`, one row per feature, or per (feature, group) when
#'   `groups` is supplied.
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "featStatsParam"),
    function(x, param, ..., groups = NULL, stats = NULL) {
        do.call(.pe_featstats,
            c(list(pe = x, groups = groups, stats = stats),
                as.list(param@param)))
    }
)


## internals --- cell/feat stats ####

.pe_cellstats <- function(pe, detection_threshold = 0, ...) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[cell_qc_stats] pe must be a parquetExprBase.")

    n_genes <- as.integer(pe@n_genes)
    thr     <- as.numeric(detection_threshold)

    # Matching Giotto's cellStatsParam: nr_feats counts entries above the
    # threshold, total_expr is the unconditional colSums.
    acc <- .stream_expr_accum(pe,
        axis  = "cell",
        thr   = thr,
        stats = c("sum", "nnz")
    )

    # `@cell_ids` is the concatenation of the substores' own cell_ids in
    # substore order, which is the axis the accumulator indexes.
    data.table::data.table(
        cells      = pe@cell_ids,
        nr_feats   = acc$nnz,
        perc_feats = acc$nnz / n_genes * 100,
        total_expr = acc$sum
    )
}

.pe_featstats <- function(pe, detection_threshold = 0, groups = NULL,
    stats = NULL, ...) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[feat_qc_stats] pe must be a parquetExprBase.")

    thr <- as.numeric(detection_threshold)
    if (!is.null(groups)) {
        return(.pe_featstats_grouped(pe, thr, groups,
            stats = stats %null% c("sum", "sumsq", "nnz", "sum_det")))
    }
    if (!is.null(stats)) {
        stop("[feat_qc_stats] `stats` selection is only available on the ",
             "grouped path (pass `groups`). The ungrouped verb has a fixed ",
             "column contract.", call. = FALSE)
    }

    n_cells <- as.integer(pe@n_cells)

    # `sum` and `sum_det` are separate accumulators because the threshold
    # picks the population for `mean_expr_det` without touching the totals
    # (Giotto: total_expr / mean_expr are unconditional rowSums / rowMeans,
    # while mean_expr_det comes from .mean_expr_det_test over detected
    # entries only). Reusing one filtered sum for both would shrink the
    # totals.
    acc <- .stream_expr_accum(pe,
        axis  = "feat",
        thr   = thr,
        stats = c("sum", "nnz", "sum_det")
    )

    data.table::data.table(
        feats         = pe@feat_ids,
        nr_cells      = acc$nnz,
        perc_cells    = acc$nnz / n_cells * 100,
        total_expr    = acc$sum,
        mean_expr     = acc$sum / n_cells,
        mean_expr_det = ifelse(acc$nnz > 0L, acc$sum_det / acc$nnz, NaN)
    )
}


# Resolve a per-cell grouping against the cell axis of the current view.
#
# A grouping is a per-cell payload, so adr/0003 applies to it: keyed by view
# position it reads the wrong entries the moment `[` narrows the store, and
# nothing checks the keying. Named by cell ID it is keyed by identity, and the
# conversion to the on-disk key happens at the point of use -- the axis map in
# `.pe_featstats_grouped()`.
#
# Cells the grouping does not name resolve to NA and drop out of the aggregate.
# That is `.pe_axis_pos_map()`'s convention for an id with no place on the axis,
# and it is already what an NA label means here, so narrowing to a few groups by
# masking the rest is a supported way to call this. No overlap at all is a
# mistake rather than an empty selection, and says so.
#
# Unnamed input is positional against the view -- what a caller holding
# `pe@cell_ids` already has -- and warns.
.pe_cell_groups <- function(pe, groups) {
    ids <- pe@cell_ids

    if (!is.null(names(groups))) {
        ord <- match(ids, names(groups))
        if (all(is.na(ord))) {
            stop("[feat_qc_stats] `groups` is named but none of its names are ",
                 "cell IDs of the current view.", call. = FALSE)
        }
        return(groups[ord])
    }

    if (length(groups) != length(ids)) {
        stop("[feat_qc_stats] `groups` must have one entry per cell of the ",
             "current view (", length(ids), "), got ", length(groups), ".",
             call. = FALSE)
    }
    warning("`groups` is unnamed and is being matched to the store by ",
            "position. Pass a vector or factor named by cell ID to match on ",
            "identity instead -- cell metadata and expression are not ",
            "guaranteed to share an order.", call. = FALSE)
    groups
}


# Grouped feature statistics: the same accumulators, partitioned by a per-cell
# grouping instead of taken over every cell. One extra key on the aggregate,
# one extra pass over nothing -- it is the same single scan.
#
# The statistics are unchanged in meaning; only the population each one is
# taken over is narrower. So `n_cells` becomes a column (it varies by group,
# where ungrouped it is a scalar) and `mean_expr` / `sd` are over ALL cells of
# the group, absent entries included. adr/0009 still governs the threshold:
# `nr_cells` sees it, the sums do not.
#
# Emits the complete feats x groups cross product, zeros filled. A feature with
# no stored value in a group has mean 0 over that group's cells, not a missing
# row -- and a consumer plotting group means needs the zeros present.
#
# Reusable well beyond markers: per-cluster mean + percent-detected is the
# input to a dot plot, and the group means are a pseudobulk matrix.
#
# `stats` names ACCUMULATORS, not output columns, and passes straight through
# to `.pe_accum_raw()` -- the same vocabulary the ungrouped verbs already select
# from. Emitted columns are whichever the requested accumulators support:
#
#   sum             -> total_expr, mean_expr
#   sum + sumsq     -> sumsq, sd
#   nnz             -> nr_cells, perc_cells
#   sum_det + nnz   -> mean_expr_det
#   (always)        -> feats, group, n_cells
#
# `total_expr` and `sumsq` are the raw accumulators under their reporting
# names, emitted alongside the derived `mean_expr` / `sd` because they are
# already in hand. Keeping them visible is what lets a caller combine groups
# by addition (they are additive; means and standard deviations are not)
# instead of reconstructing them through a square root.
#
# Selection is at accumulator granularity only: columns sharing an accumulator
# are free, so they are emitted together rather than trimmed. A caller that
# wants only group moments asks for `c("sum", "sumsq")` and skips two of the
# four accumulators; it does not describe the columns it wants and leave this
# function to work backwards to a plan.
.pe_featstats_grouped <- function(pe, thr, groups,
    stats = c("sum", "sumsq", "nnz", "sum_det")) {
    k <- pos <- NULL   # NSE bindings

    stats <- match.arg(stats, several.ok = TRUE)

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)
    groups <- .pe_cell_groups(pe, groups)

    # `droplevels` so an unused level cannot surface as a group of zero cells.
    g <- droplevels(if (is.factor(groups)) groups else factor(groups))
    lvls <- levels(g)
    if (length(lvls) < 1L) {
        stop("[feat_qc_stats] `groups` has no non-empty levels.", call. = FALSE)
    }
    codes <- as.integer(g)

    # Cells per group over the view. NA-group cells are excluded here and by
    # the aggregate's inner join alike, so both sides see one population.
    nk <- as.numeric(tabulate(codes, nbins = length(lvls)))
    names(nk) <- lvls

    # View position -> on-disk cell key, then attach the integer code. Keyed
    # by on-disk id per adr/0003, which is what `.pe_accum_raw()` joins on.
    # Indexing by `pos` is sound because `.pe_cell_groups()` has already put the
    # codes in the view's cell order by identity; the map is the position ->
    # on-disk-key conversion the ADR asks consumers to do at the point of use.
    map <- data.table::copy(.pe_axis_pos_map(pe, "cell"))
    map[, k := codes[pos]]
    map <- map[!is.na(k)]
    data.table::setnames(map, "key_id", "row_id")
    by_cell <- map[, c(
        if ("source_id" %in% names(map)) "source_id", "row_id", "k"
    ), with = FALSE]

    acc <- lapply(stats, function(nm) matrix(0, n_genes, length(lvls)))
    names(acc) <- stats

    agg <- .pe_accum_raw(pe,
        axis = "feat", thr = thr, stats = stats, by_cell = by_cell
    )
    if (!is.null(agg)) {
        idx <- cbind(as.integer(agg$pos), as.integer(agg$k))
        for (nm in stats) acc[[nm]][idx] <- as.numeric(agg[[nm]])
    }

    # Long form, groups varying slowest so `feats` cycles within a group --
    # the order `matrix` unrolls in, so the accumulators drop straight in.
    nn <- rep(nk, each = n_genes)
    out <- data.table::data.table(
        feats   = rep(pe@feat_ids, times = length(lvls)),
        group   = rep(lvls, each = n_genes),
        n_cells = nn
    )

    if ("sum" %in% stats) {
        gene_sum <- as.numeric(acc$sum)
        out[, "total_expr" := gene_sum]
        out[, "mean_expr" := gene_sum / nn]

        if ("sumsq" %in% stats) {
            gene_sumsq <- as.numeric(acc$sumsq)
            out[, "sumsq" := gene_sumsq]
            # Sum-of-squares identity, clamped: the subtraction can go slightly
            # negative when the mean dominates the spread, and a negative
            # variance would surface downstream as an NaN standard deviation.
            gene_var <- ifelse(nn > 1,
                pmax((gene_sumsq - gene_sum * gene_sum / nn) / (nn - 1), 0),
                0
            )
            out[, "sd" := sqrt(gene_var)]
        }
    }
    if ("nnz" %in% stats) {
        gene_nnz <- as.numeric(acc$nnz)
        out[, "nr_cells" := as.integer(gene_nnz)]
        out[, "perc_cells" := ifelse(nn > 0, gene_nnz / nn * 100, NaN)]

        if ("sum_det" %in% stats) {
            out[, "mean_expr_det" := ifelse(
                gene_nnz > 0, as.numeric(acc$sum_det) / gene_nnz, NaN
            )]
        }
    }

    out[]
}

# One grouped accumulator pass over the triplet stream, shared by every
# per-axis statistic verb.
#
# Bounded, but for three different reasons. Ungrouped with `@post_ops` empty,
# Acero streams the aggregate to its own budget -- bounded by construction,
# whatever the store's size. Grouped, the join in front of the aggregate is
# not, so the input is windowed instead and the bound comes from
# `.recommend_chunk_size()`. With post-ops it must materialize rows, same
# window, weaker still. Assume the windowed shape when reasoning about a new
# caller: it is the weaker guarantee, and the one an incorrect window breaks.
#
# Neither guarantee survives leaving the framework. `storeRead(output =
# "query")` is lazy on purpose; `collect()` or `as.data.frame()` on it pulls
# the whole store into memory with none of the above applying.
#
# Reached by every statistic over expression values: QC stats, HVF, filtering,
# normalization scale factors, and (through the grouped feature-stats verb)
# marker moments. Mechanism -- the accumulator vocabulary, `by_cell` grouping,
# the two execution shapes, the union `source_id` key and the on-disk-id fold
# -- is in design.Rmd, "Expression Statistics". Only the local contract is
# repeated here.
#
# Returns vectors indexed by position along `axis` in the current view
# (n_genes long for "feat", n_cells for "cell"):
#
#   sum      sum(value)                       every stored entry
#   sumsq    sum(value^2)                     every stored entry
#   nnz      count(value > thr)               detection count
#   sum_det  sum(value where value > thr)     detected entries only
#
# `thr` is a detection predicate: it selects which entries COUNT, and never
# reduces a magnitude that participates. So `sum` / `sumsq` are unconditional
# and only `nnz` / `sum_det` see it (adr/0009). `inclusive = TRUE` switches
# those two to `>=`, which is what filtering means by a threshold.

.stream_expr_accum <- function(pe,
    axis = c("feat", "cell"),
    thr = 0,
    stats = c("sum", "sumsq", "nnz", "sum_det"),
    inclusive = FALSE) {

    axis  <- match.arg(axis)
    stats <- match.arg(stats, several.ok = TRUE)

    n_out <- as.integer(
        if (identical(axis, "feat")) pe@n_genes else pe@n_cells
    )
    out <- lapply(stats, function(nm) {
        if (identical(nm, "nnz")) integer(n_out) else numeric(n_out)
    })
    names(out) <- stats

    folded <- .pe_accum_raw(pe,
        axis = axis, stats = stats, thr = thr, inclusive = inclusive
    )
    if (is.null(folded)) return(out)

    for (nm in stats) {
        out[[nm]][folded$pos] <- if (identical(nm, "nnz")) {
            as.integer(folded[[nm]])
        } else {
            as.numeric(folded[[nm]])
        }
    }

    out
}


# Shared front of the accumulation: build the aggregate, run it on whichever
# execution shape the chain allows, and resolve on-disk ids to axis positions.
# Returns a LONG data.table -- `pos` (+ `k` when grouped), one column per
# statistic -- or NULL when nothing survived. Absent combinations are simply
# missing rows; it is the caller's job to decide what their absence means.
#
# `.stream_expr_accum()` layers the per-axis vector shape on top. Anything
# needing a second grouping key (marker detection: per-(feature, cluster)
# moments) reads this directly, because the vector shape cannot hold it.
#
# `by_cell` adds a cell-side grouping variable to the key: a data.table of
# `row_id` (ON-DISK, per ADR 0003), `k` (integer group code), and `source_id`
# on a union. It is joined in before the aggregate, so it narrows as well as
# groups -- cells absent from the table never reach the aggregate, which is how
# NA-group cells are excluded.
#
# Supplying it windows the INPUT -- see `.pe_accum_acero_windowed()`. Acero's
# bound on the ungrouped aggregate comes from the group key being small and the
# scan never being retained; the join is the one thing that breaks that,
# because the rows it emits are O(nonzeros) rather than O(groups). Windowing is
# not a lesser workaround for it: the accumulators are additive over cell
# windows, so a window is exact, its peak is chosen rather than discovered, and
# a contiguous `row_id` range prunes row groups where a gene-side narrowing
# cannot (the store is sorted cell-major). An engine that spilled instead would
# survive the join without ever avoiding it.
#
# A caller batching the FEATURE axis to dodge the same memory is paying for it
# the expensive way: gene ids are not the sort key, so every batch rescans the
# store in full and the cost is linear in batch count, not in genes per batch.
# Ask for all the features at once and let the window do the bounding.
#
# `k` is an integer code, never a label string. Aggregating or joining strings
# in Acero at scale leaves dangling `utf8_view` buffers; the caller maps codes
# back to labels after collect. (`source_id` is already a string in the union
# group key, unavoidably -- it is a hive partition value. The `by_cell` join
# therefore adds a string join column on unions and none at all on a single
# store, which is the common case. Worth measuring before assuming it is free
# on a wide union.)
.pe_accum_raw <- function(pe,
    axis = c("feat", "cell"),
    stats = c("sum", "sumsq", "nnz", "sum_det"),
    thr = 0,
    inclusive = FALSE,
    by_cell = NULL) {

    if (!inherits(pe, "parquetExprBase")) {
        stop("[.pe_accum_raw] pe must be a parquetExprBase.", call. = FALSE)
    }
    axis  <- match.arg(axis)
    stats <- match.arg(stats, several.ok = TRUE)
    thr   <- as.numeric(thr)

    value <- pos <- NULL   # NSE bindings

    key <- if (identical(axis, "feat")) "col_id" else "row_id"

    # `!!thr` folds the threshold in as a literal, so these stay plain R calls
    # -- they splice into an arrow `summarise()` and evaluate in a data.table
    # `j` unchanged.
    #
    # `inclusive` picks `>=` over `>`. The statistic verbs want strict (a
    # detection threshold of 0 counts anything nonzero), while filtering wants
    # inclusive -- Giotto's `expression_threshold = 1` means "expressed if the
    # count is at least 1".
    detected <- if (isTRUE(inclusive)) {
        rlang::expr(value >= !!thr)
    } else {
        rlang::expr(value > !!thr)
    }
    all_exprs <- list(
        sum     = rlang::expr(sum(value, na.rm = TRUE)),
        sumsq   = rlang::expr(sum(value * value, na.rm = TRUE)),
        nnz     = rlang::expr(sum(as.integer(!!detected), na.rm = TRUE)),
        sum_det = rlang::expr(
            sum(value * as.integer(!!detected), na.rm = TRUE)
        )
    )
    aggr_exprs <- all_exprs[stats]

    # `source_id` only earns its place in the group key on a union, where
    # row_id restarts per substore and col_id is meaningful only relative to
    # its own substore. A single store has exactly one, so including it just
    # widens the key with a constant string column -- measured at ~50% on the
    # aggregate, growing with the scan.
    is_union <- inherits(pe, "unionParquetExprStore")
    has_by   <- !is.null(by_cell)
    grp      <- c(if (is_union) "source_id", key, if (has_by) "k")
    join_by  <- if (is_union) c("source_id", "key_id") else "key_id"

    if (length(pe@post_ops) == 0L) {
        # Acero path: the aggregate runs as an arrow plan either way. `by_cell`
        # decides only whether it runs once over the whole store or once per
        # cell window -- what the join needs bounding, it does not need
        # materializing, so keep the execution in Acero and window the input.
        agg <- if (has_by) {
            .pe_accum_acero_windowed(pe,
                grp = grp, aggr_exprs = aggr_exprs, by_cell = by_cell)
        } else {
            .pe_agg_collect(pe, grp = grp, aggr_exprs = aggr_exprs)
        }
    } else {
        # R path: the chain has to run on materialized rows, so this one is
        # chunked rather than collected whole.
        agg <- .pe_accum_chunked_dt(pe,
            grp = grp, aggr_exprs = aggr_exprs, by_cell = by_cell
        )
    }
    if (is.null(agg) || nrow(agg) == 0L) return(NULL)

    # Arrow returns int64 for integer sums; normalize before any R arithmetic.
    for (nm in stats) {
        data.table::set(agg, j = nm, value = as.numeric(agg[[nm]]))
    }

    # Resolve on-disk ids to axis positions. This is a lookup, not a filter:
    # `storeRead()` applies exact axis predicates (`.pe_axis_pred_exprs()`), so
    # out-of-view entries never reach the aggregate and the inner join has
    # nothing to drop -- it stays inner only as a backstop against a store
    # whose on-disk ids fall outside its declared index. Both sides are small:
    # one row per (substore, axis position), never per stored entry.
    data.table::setnames(agg, key, "key_id")
    agg <- merge(agg, .pe_axis_pos_map(pe, axis), by = join_by)
    if (nrow(agg) == 0L) return(NULL)

    # Fold by position: real summation on the feature axis (a feature spans
    # substores), a no-op on the cell axis (cell axes are disjoint, so every
    # group is a singleton). `k` joins the fold key when present -- a feature
    # spans substores within a cluster, not across clusters.
    agg[, lapply(.SD, sum),
        by = c("pos", if (has_by) "k"), .SDcols = stats]
}

# Window for the chunked R-side pass. Not `.exprbase_chunk_size()`: this pass
# collects a triplet frame -- row_id + col_id + value + source_id, plus the
# arrow collect buffer and the data.table copy -- measured at 47-54 B/nonzero
# across three shapes, so it asks for 48 and splits ~4x sooner. `k = 0`
# because there is no PCA sketch to reserve against.
.pe_accum_chunk_size <- function(pe) {
    .pe_window_cells(pe, bytes_per_nz = 48, k = 0L)
}


# Cell-windowed Acero accumulation, for a GROUPED pass whose chain still
# lowers. One arrow plan per cell window instead of one over the whole store.
#
# Why the grouped pass needs windowing when the ungrouped one does not: the
# aggregate itself is O(groups) either way, but `by_cell` puts a join in front
# of it whose output is O(nonzeros). Acero has no spill, so at atlas scale that
# is the whole failure. Windowing bounds the join's input instead, which the
# accumulators permit for free -- they are additive over cell windows, so the
# caller's fold sums the partials and the answer is exact rather than
# approximate. Nothing here materializes; each window's aggregate is the only
# thing that crosses into R.
#
# Measured on 5k genes x 50k cells x 12 clusters, 20 windows: 1.0s here vs 2.1s
# routing the same windows through `.pe_accum_chunked_dt()`, which is correct
# but pays a collect + data.table group-by per window. Window the input, keep
# the execution in Acero.
#
# Windows are not free -- the same store is 0.09s in one window and 1.0s in
# twenty, so roughly 45ms of plan setup and fold per window on top of the scan
# it would have done anyway. That is an argument for taking the LARGEST window
# the budget allows, which is what `.pe_accum_chunk_size()` already returns; it
# is not an argument for a fixed window count. Do not tighten the window
# speculatively.
#
# A budget that covers the view yields one window, which is one plan over the
# whole store -- the pre-existing behaviour, so small and mid-size stores are
# not charged for the bound. Unions iterate substores for the same reason as
# the R path: `row_id` restarts per substore, so a global cell range is not a
# contiguous range and would not prune.
#
# Partials are folded as they arrive rather than collected and folded at the
# end. Keeping one per window would make the retained state O(groups x
# windows) -- n_genes x n_clusters rows each, so tightening the window to save
# memory would cost memory, which is backwards. Eager folding is the same
# arithmetic the caller's fold does, so the result is unchanged and what is
# held is O(groups) whatever the window count.
.pe_accum_acero_windowed <- function(pe, grp, aggr_exprs, by_cell) {
    by_cols <- setdiff(names(by_cell), "k")
    # Built once, not per window: it is the join's build side, one row per
    # grouped cell, and it is the same table for every window.
    bc_tab  <- arrow::as_arrow_table(as.data.frame(by_cell))

    acc <- NULL
    for (d in .pe_windows(pe, .pe_accum_chunk_size(pe))) {
        acc <- .pe_fold_partial(acc,
            .pe_agg_collect(.pe_window_store(d),
                grp = grp, aggr_exprs = aggr_exprs,
                join_tab = bc_tab, join_cols = by_cols),
            grp = grp, cols = names(aggr_exprs))
    }
    acc
}


# ---- shared by both windowed accumulators ----------------------------------

# One grouped aggregate as an arrow plan, collected. The whole-store branch of
# `.pe_accum_raw()` and each window of `.pe_accum_acero_windowed()` are the
# same plan differing only by the join, so they are the same function.
#
# `join_tab` is the join's BUILD side -- small by construction (one row per
# grouped cell), streamed against by the scan. Inner, so a cell absent from it
# drops, which is how NA-group cells are excluded without a predicate.
.pe_agg_collect <- function(pe, grp, aggr_exprs,
    join_tab = NULL, join_cols = NULL) {

    q <- storeRead(pe, output = "query")
    if (!is.null(join_tab)) {
        q <- dplyr::inner_join(q, join_tab, by = join_cols)
    }
    q |>
        dplyr::group_by(!!!rlang::syms(grp)) |>
        dplyr::summarise(!!!aggr_exprs) |>
        dplyr::collect() |>
        data.table::as.data.table()
}

# Fold one window's partial into the running total.
#
# Eager, not a list of partials collected and reduced at the end: a partial is
# one row per group, so retaining one per window would make the held state
# O(groups x windows) -- tightening the window to save memory would cost
# memory. Sound because the accumulators are additive and each window covers a
# disjoint set of cells; this is the same arithmetic the caller's
# fold-by-position does, applied earlier.
#
# Coerced to double on the way in: arrow hands back int64 for integer sums, and
# a running total should not depend on bit64 dispatch. Counts stay exact well
# past 2^53.
#
# NOT bitwise invariant, and cannot be. Folding on arrival reassociates the
# summation relative to reducing all the partials at once, so a float
# accumulator (`sum`, `sumsq`) can land 1 ULP away -- measured at 1.0-1.2 ULP,
# max relative 2.7e-16, on a union with a multi-window post-op chain, where the
# second fold level adds more reassociation. Integer accumulators (`nnz`) are
# exact. Do not write a test that demands bitwise equality across window counts;
# `expect_equal`'s default tolerance is the right strictness here.
.pe_fold_partial <- function(acc, part, grp, cols) {
    if (is.null(part) || nrow(part) == 0L) return(acc)
    for (nm in cols) {
        data.table::set(part, j = nm, value = as.numeric(part[[nm]]))
    }
    if (is.null(acc)) return(part)
    data.table::rbindlist(list(acc, part))[
        , lapply(.SD, sum), by = grp, .SDcols = cols]
}


# Chunked R-side accumulation, for chains that cannot be lowered.
#
# The Acero branch streams by construction; this one has to materialize rows
# to run the chain, so it reads a cell window at a time and aggregates each
# window, keeping peak memory at chunk_size cells rather than the whole store.
# Partial aggregates are safe to concatenate without reconciliation: the
# accumulators are additive, and the caller's fold-by-position sums whatever
# duplicate (source_id, key) rows the windowing produced.
#
# Windows come from `.pe_windows()`, the shared seam -- see its header for why
# they are cell ranges rather than row counts or a cell set. The one thing
# specific to this path: `[` also slices `@post_ops` down to the window's cells,
# so the chain applies with exactly the state that window needs.
#
# Note this path always narrows with `[`, where the Acero path can skip it on a
# single full-substore window. `.pe_window_store()` handles that, and the
# post-op slice is what makes taking it here correct rather than merely cheaper.
.pe_accum_chunked_dt <- function(pe, grp, aggr_exprs, by_cell = NULL) {
    dt_call <- as.call(c(quote(list), aggr_exprs))
    by_cols <- if (!is.null(by_cell)) setdiff(names(by_cell), "k")

    acc <- NULL
    for (d in .pe_windows(pe, .pe_accum_chunk_size(pe))) {
        w  <- .pe_window_store(d)
        df <- data.table::as.data.table(
            dplyr::collect(storeRead(w, output = "query")))
        if (nrow(df) > 0L) {
            df <- .pe_apply_post_ops_df(df, w@post_ops)
            # Inner, matching the Acero branch's join: a cell absent from
            # `by_cell` is out of the grouping and drops here too.
            if (!is.null(by_cell)) df <- merge(df, by_cell, by = by_cols)
        }
        if (nrow(df) == 0L) next
        acc <- .pe_fold_partial(acc, df[, eval(dt_call), by = grp],
            grp = grp, cols = names(aggr_exprs))
    }
    acc
}

# Build the `(source_id, key_id) -> pos` map for one axis. `key_id` is the
# on-disk id, `pos` the position in the view's axis.
#
# Slicing state lives entirely on the substores: `[` on a union pushes both
# axes down and rebuilds the parent, which carries no `@cell_idx`/`@gene_idx`
# of its own. So each substore's index vector plus (for cells) its offset into
# the union axis is the whole mapping.
.pe_axis_pos_map <- function(pe, axis) {
    pos <- NULL   # NSE binding

    is_union <- inherits(pe, "unionParquetExprStore")
    is_feat  <- identical(axis, "feat")
    axis_ids <- if (is_feat) pe@feat_ids else pe@cell_ids

    map <- data.table::rbindlist(lapply(.exprbase_substores(pe), function(se) {
        s <- se$store
        # `@feat_ids[k]` / `@cell_ids[k]` names the entry at `@gene_idx[k]` /
        # `@cell_idx[k]` -- `[` narrows the id vector and the index vector
        # together -- so ids and on-disk keys stay aligned per substore.
        if (is_feat) {
            ki  <- if (length(s@gene_idx) > 0L) s@gene_idx
                   else seq_len(as.integer(s@n_genes))
            ids <- s@feat_ids
        } else {
            ki  <- if (length(s@cell_idx) > 0L) s@cell_idx
                   else seq_len(as.integer(s@n_cells))
            ids <- s@cell_ids
        }

        # Resolve by identifier, not by position. On-disk `col_id` layout is
        # not guaranteed to agree across substores (same reason the op-slice
        # registry keys the gene axis on `feat_id`), and an id that has no
        # place on the union's axis resolves to NA and drops below -- which is
        # what partial alignment should do.
        #
        # A single store is its own substore, so its ids ARE the view's ids in
        # view order and the lookup is the identity; skip it, along with the
        # `source_id` column the single-store join does not key on.
        if (!is_union) {
            return(data.table::data.table(
                key_id = as.integer(ki), pos = seq_along(ki)))
        }
        data.table::data.table(
            source_id = rep_len(as.character(s@uid), length(ki)),
            key_id    = as.integer(ki),
            pos       = match(ids, axis_ids)
        )
    }))
    map[!is.na(pos)]
}

# hvf ####
# HVF-relevant stats for parquetExprStore-backed expression.
# Dispatches via Giotto's analyzeData(x, param) generic. Methods return
# the per-feature stats data.table.
#
#   analyzeData(parquetExprStore, covLoessParam)
#       -> per-feature stats including cov_diff (residual COV above
#          a LOESS fit of cov ~ log(mean_expr))
#   analyzeData(parquetExprStore, covGroupsParam)
#       -> per-feature stats including cov_group_zscore (within-bin
#          COV z-score)
#   analyzeData(parquetExprStore, varParam)
#       -> per-feature variance of analytic Poisson Pearson residuals,
#          computed from raw counts (see the varParam section below)
#
# The COV methods share one stats pass, `.stream_gene_stats()`, which
# reduces the triplet stream to n_genes-sized sum / sumsq / nnz vectors;
# LOESS and bin-zscore then run on those. That pass has two shapes: an
# arrow-native aggregate when every post-op can be lowered to arrow
# compute, and a generic per-substore R-side pass otherwise.

## covLoessParam ####

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "covLoessParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        stats <- .stream_gene_stats(x, expression_threshold = thr)

        # Match Giotto: drop zero-detection features before fitting
        nr_cells <- cov <- pred_cov <- cov_diff <- mean_expr <- NULL
        stats <- stats[nr_cells > 0]

        loess_fit <- stats::loess(cov ~ log(mean_expr), data = stats)
        stats[, pred_cov := stats::predict(loess_fit, newdata = stats)]
        stats[, cov_diff := cov - pred_cov]
        stats[, pred_cov := NULL]
        data.table::setorder(stats, -cov_diff)
        stats
    }
)


## covGroupsParam ####

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "covGroupsParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        stats <- .stream_gene_stats(x, expression_threshold = thr)

        # NSE bindings
        nr_cells <- cov <- expr_groups <- cov_group_zscore <- NULL

        # drop zero-detection features before binning
        stats <- stats[nr_cells > 0]

        # Quantile-bin by mean expression. If too many tied breaks (lots of
        # zero-mean genes), recompute on the strictly positive subset and
        # set the leading break to 0 so all-zero genes still bin into
        # group_1 -- matches Giotto's in-memory .calc_cov_group_hvf.
        n_groups <- as.integer(param$nr_expression_groups)
        prob_seq <- seq(0, 1, by = 1 / n_groups)
        prob_seq[length(prob_seq)] <- 1
        expr_group_breaks <- stats::quantile(stats$mean_expr, probs = prob_seq)
        if (any(duplicated(expr_group_breaks))) {
            m <- stats$mean_expr
            expr_group_breaks <- stats::quantile(m[m > 0], probs = prob_seq)
            expr_group_breaks[[1L]] <- 0
        }

        expr_groups_lbl <- cut(
            stats$mean_expr,
            breaks         = expr_group_breaks,
            labels         = paste0("group_", seq_len(n_groups)),
            include.lowest = TRUE
        )
        stats[, expr_groups := expr_groups_lbl]
        stats[, cov_group_zscore := scale(cov), by = expr_groups]
        stats[, expr_groups := NULL]
        stats
    }
)


## varParam ####
# analytic Pearson residual variance 
#
# `calculateHVF(method = "var_p_resid")` ranks features by the `var` column
# this returns.  In Giotto's in-memory path that is `rowVars()` of whatever
# matrix it was handed, which only equals Pearson residual variance if the
# caller arranged for the expression slot to hold residuals -- with the
# default `expression_values = "normalized"` it is the variance of
# log-normalized values instead.
#
# The streaming backend computes the residual variance directly from RAW
# counts, which is both what the method name means and the only form that
# works here: Pearson residuals are dense (every zero maps to a nonzero
# residual), so a residual matrix would be n_genes x n_cells -- 3.1e9 rows
# for Atera-scale data.  See `.stream_pearson_resid_var()` for how the zero
# block is folded in analytically instead of being materialized.
#
# Consequences worth knowing:
#   * `expression_values` defaults to "raw" here, not "normalized".  Pearson
#     residuals ARE the normalization -- the depth term lives inside
#     `mu = g_i * c_j / T` -- so normalizing first double-corrects.  Measured
#     on synthetic data with 20 injected overdispersed genes: raw counts
#     recover 20/20 with non-variable genes calibrated at 1.000, whereas
#     running the same formula on log2(1+libnorm) values recovers 0/20.
#   * the scale is absolute: 1.0 means "no more variable than Poisson
#     sampling noise", which is what makes `calculateHVF`'s default
#     `var_threshold = 1.5` meaningful rather than arbitrary.

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "varParam"),
    function(x, param, ...) {
        ev <- param$expression_values %null% "raw"
        if (!identical(ev, "raw")) {
            warning("[analyzeData(parquetExprBase, varParam)] ",
                    "expression_values = '", ev, "' is ignored: Pearson ",
                    "residual variance is defined on raw counts, so the ",
                    "streaming backend reads raw values and does not apply ",
                    "the normalization recipe.", call. = FALSE)
        }
        if (!inherits(x, "parquetExprStore")) {
            stop("[analyzeData(parquetExprBase, varParam)] union stores are ",
                 "not supported yet: per-cell totals key on ",
                 "(source_id, row_id) across substores and the gene-axis ",
                 "remap resolves against a single store's @gene_idx.",
                 call. = FALSE)
        }
        .stream_pearson_resid_var(x, size_factors = param$size_factors,
                                  theta = param$theta %null% 100)
    }
)

## internals --- hvf ####

# Per-gene variance of analytic Pearson residuals (Lause/Kobak), streamed.
#
#   mu_ij = g_i * c_j / T        d_ij = sqrt(mu_ij + mu_ij^2 / theta)
#   z_ij  = clamp((x_ij - mu_ij) / d_ij,  -sqrt(n), +sqrt(n))
#
# with g_i the gene total, c_j the cell total (or a supplied size factor), and
# T the grand total. Same formula as Giotto's `.prnorm()`, so the streaming
# and in-memory routes cannot drift.
#
# Residuals are dense, so the statistic is built as "assume every entry is
# zero, then correct the stored nonzeros":
#
#   zeros:   z0_ij = clamp(-mu_ij / d_ij, -sqrt(n), Inf)
#   stored:  dz  = z_ij - z0_ij,   dz2 = z_ij^2 - z0_ij^2
#
# With theta = Inf the zero block collapsed to `sum_j z0^2 = g_i`, a scalar per
# gene. Finite theta breaks that -- `mu/(1 + mu/theta)` is not separable into
# gene- and cell-side factors -- so the zero block is summed explicitly over
# cells, which is O(n_genes * n_cells) and hopeless at bin1 scale
# (2.6e4 * 5.0e6).
#
# The way out: the term depends on c_j only through mu_ij, so identical cell
# totals contribute identically. Collapsing the totals to unique values with
# multiplicities makes it O(n_genes * n_unique_totals). Measured on Stereo-seq
# tissue.gef: bin1 has 60 unique totals across 5,043,144 bins (1.6e6 products,
# an 80,000x reduction), cellbin 2,789 across 7,527. The collapse is largest
# exactly where the naive form is unaffordable, because a bin1 bin is a single
# DNB holding a handful of transcripts.
#
# @post_ops is deliberately NOT applied: the values must be counts.

.stream_pearson_resid_var <- function(pe, size_factors = NULL, theta = 100) {
    # NSE bindings
    row_id <- col_id <- value <- g <- cc <- mu <- dz <- dz2 <- NULL
    sum_dz <- sum_dz2 <- sum_z <- sum_z2 <- var <- feats <- NULL
    dd <- zc <- z0c <- s0 <- s02 <- w <- mean_expr <- NULL

    n_genes <- as.integer(pe@n_genes)

    # Pass 1: per-gene and per-cell totals (both pushed-down aggregates).
    gt <- storeRead(pe, output = "query") |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(g = sum(value, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()
    if (nrow(gt) == 0L) {
        return(data.table::data.table(feats = pe@feat_ids,
                                      var = numeric(n_genes),
                                      mean_expr = numeric(n_genes)))
    }
    ct <- storeRead(pe, output = "query") |>
        dplyr::group_by(row_id) |>
        dplyr::summarise(cc = sum(value, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()

    # Counts sanity: the model assumes integer counts, and nothing about the
    # store can guarantee that (a normalized matrix may have been ingested).
    if (any(abs(gt$g - round(gt$g)) > 1e-8)) {
        warning("[analyzeData(varParam)] expression values do not look like ",
                "integer counts; Pearson residual variance assumes raw ",
                "counts and will be misleading on transformed values.",
                call. = FALSE)
    }

    # Cells with no counts have mu = 0 for every gene, leaving the residual
    # undefined; drop them and shrink n accordingly.
    if (!is.null(size_factors)) {
        sf <- as.numeric(size_factors)
        if (length(sf) != as.integer(pe@n_cells)) {
            stop("[analyzeData(varParam)] size_factors must have one entry ",
                 "per cell (", pe@n_cells, "), got ", length(sf), ".",
                 call. = FALSE)
        }
        cvals <- sf[sf > 0]
    } else {
        cvals <- ct$cc[ct$cc > 0]
    }
    n_eff <- length(cvals)
    if (n_eff < 2L) {
        return(data.table::data.table(feats = pe@feat_ids,
                                      var = numeric(n_genes),
                                      mean_expr = numeric(n_genes)))
    }
    Tt   <- sum(gt$g)
    thta <- as.numeric(theta)
    if (!is.finite(thta) || thta <= 0) {
        stop("[analyzeData(varParam)] theta must be a positive finite number.",
             call. = FALSE)
    }
    # Clip bound follows `.prnorm()`: +/- sqrt(n) over the cells in play.
    clip <- sqrt(n_eff)

    # Pass 2a: the all-zero block, summed over UNIQUE cell totals.
    # z0 = -mu / sqrt(mu + mu^2/theta), clamped below at -clip. |z0| is at most
    # sqrt(theta), so the clamp only ever binds for absurdly small n; it is
    # applied anyway so the two backends agree by construction rather than by
    # a size argument.
    uc <- data.table::data.table(cc = cvals)[, .(w = .N), keyby = "cc"]
    gvec <- as.numeric(gt$g)
    n_u <- nrow(uc)
    s0 <- numeric(length(gvec))
    s02 <- numeric(length(gvec))
    # Chunk the gene axis so the gene x unique-total outer stays bounded
    # regardless of how many distinct totals a dataset has.
    blk <- max(1L, as.integer(5e6 %/% max(1L, n_u)))
    for (lo in seq.int(1L, length(gvec), by = blk)) {
        hi <- min(lo + blk - 1L, length(gvec))
        mu <- outer(gvec[lo:hi], uc$cc / Tt)          # gene x unique total
        z0 <- -mu / sqrt(mu + mu * mu / thta)
        z0[z0 < -clip] <- -clip
        s0[lo:hi]  <- as.vector(z0 %*% uc$w)
        s02[lo:hi] <- as.vector((z0 * z0) %*% uc$w)
    }
    z0dt <- data.table::data.table(col_id = gt$col_id, s0 = s0, s02 = s02)

    # Pass 2b: per-nonzero corrections, joined and aggregated in Acero. Each
    # stored value replaces its assumed-zero contribution, so accumulate the
    # difference.
    gt_a <- arrow::as_arrow_table(data.frame(
        col_id = as.integer(gt$col_id), g = as.numeric(gt$g)))
    ct_a <- arrow::as_arrow_table(data.frame(
        row_id = as.integer(ct$row_id), cc = as.numeric(ct$cc)))
    corr <- storeRead(pe, output = "query") |>
        dplyr::left_join(gt_a, by = "col_id") |>
        dplyr::left_join(ct_a, by = "row_id") |>
        dplyr::mutate(mu = g * cc / !!Tt) |>
        dplyr::mutate(dd = sqrt(mu + mu * mu / !!thta)) |>
        dplyr::mutate(
            zc  = pmin(pmax((value - mu) / dd, -!!clip), !!clip),
            z0c = pmax(-mu / dd, -!!clip)
        ) |>
        dplyr::mutate(dz = zc - z0c, dz2 = zc * zc - z0c * z0c) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(sum_dz  = sum(dz,  na.rm = TRUE),
                         sum_dz2 = sum(dz2, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()

    res <- merge(gt, corr, by = "col_id", all.x = TRUE)
    res <- merge(res, z0dt, by = "col_id", all.x = TRUE)
    res[is.na(sum_dz),  sum_dz  := 0]
    res[is.na(sum_dz2), sum_dz2 := 0]
    res[, sum_z  := s0 + sum_dz]
    res[, sum_z2 := s02 + sum_dz2]
    res[, var := (sum_z2 - sum_z * sum_z / n_eff) / (n_eff - 1L)]
    # Genes with no counts have undefined residuals; report 0 variance, which
    # is what rowVars() gives for an all-zero row in the in-memory path.
    res[g <= 0, var := 0]

    out <- numeric(n_genes)
    idx <- .pe_remap_col(res$col_id, pe)
    keep <- !is.na(idx)
    out[idx[keep]] <- as.numeric(res$var[keep])

    # Per-gene mean of the raw counts, for the mean-vs-variance diagnostic.
    # The gene totals are already in hand, so this costs nothing.
    mout <- numeric(n_genes)
    mout[idx[keep]] <- as.numeric(res$g[keep]) / as.numeric(pe@n_cells)

    dt <- data.table::data.table(
        feats = pe@feat_ids, var = out, mean_expr = mout
    )
    data.table::setorder(dt, -var)
    dt
}


# Shared tail for both stats paths: derive mean / sd / cov from the summed
# accumulators and assemble the result. Variance via the SS identity,
# σ² = (Σx² − n·μ²)/(n−1).
#
# Also the one place that can cheaply tell whether the values were normalized.
# There is no requirement that they be -- a store may hold values normalized on
# disk, with an empty op chain -- but per-gene totals that are all integers
# almost certainly mean raw counts reached a statistic that is normally taken
# on normalized values. `gene_sum` is already computed, so the check is free.
.gene_stats_dt <- function(pe, gene_sum, gene_sumsq, gene_nnz,
                            n_cells, n_genes, .warn_raw = TRUE) {
    if (isTRUE(.warn_raw)) {
        nzs <- gene_sum[gene_sum > 0]
        if (length(nzs) > 0L && all(abs(nzs - round(nzs)) < 1e-8)) {
            warning("[analyzeData] expression values look like raw integer ",
                    "counts. COV-based feature statistics are normally taken ",
                    "on normalized values; run normalizeGiotto() first, or ",
                    "ignore this if the store already holds normalized data.",
                    call. = FALSE)
        }
    }

    gene_mean <- gene_sum / n_cells
    gene_var  <- if (n_cells > 1L) {
        pmax((gene_sumsq - gene_sum * gene_sum / n_cells) / (n_cells - 1L), 0)
    } else {
        numeric(n_genes)
    }
    gene_sd  <- sqrt(gene_var)
    gene_cov <- ifelse(gene_mean > 0, gene_sd / gene_mean, NaN)

    data.table::data.table(
        feats      = pe@feat_ids,
        nr_cells   = gene_nnz,
        total_expr = gene_sum,
        mean_expr  = gene_mean,
        sd         = gene_sd,
        cov        = gene_cov
    )
}

# No normalization requirement and no required op chain.
# `.stream_feat_accum()` owns the pass; this is the HVF column set.

.stream_gene_stats <- function(pe, expression_threshold = 0) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[.stream_gene_stats] pe must be a parquetExprBase.")

    thr     <- as.numeric(expression_threshold)
    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)

    # `sum_det` is not requested: Giotto's COV statistics are built from the
    # unconditional mean and sd, with the threshold entering only the
    # detection count (see .calc_expr_general_stats).
    acc <- .stream_expr_accum(pe, 
        axis = "feat", 
        thr = thr,
        stats = c("sum", "sumsq", "nnz")
    )

    # The raw-counts check is only meaningful when nothing is queued in
    # either phase -- with any op the values are transformed by construction.
    .gene_stats_dt(pe, acc$sum, acc$sumsq, acc$nnz, n_cells, n_genes,
                   .warn_raw = length(pe@ops) == 0L &&
                               length(pe@post_ops) == 0L)
}