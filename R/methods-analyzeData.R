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
#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "featStatsParam"),
    function(x, param, ...) {
        do.call(.pe_featstats, c(list(pe = x), as.list(param@param)))
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

.pe_featstats <- function(pe, detection_threshold = 0, ...) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[feat_qc_stats] pe must be a parquetExprBase.")

    n_cells <- as.integer(pe@n_cells)
    thr     <- as.numeric(detection_threshold)

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

# One grouped accumulator pass over the triplet stream, shared by every
# per-axis statistic verb.
#
# Only use in memory-safe chunks. This is not memory-safe when there are
# any post-ops applied (which should be assumed for generalizability)
#
# Used by:
# - QC stats verbs (`analyzeData(featStatsParam / cellStatsParam)`)
# - HVF stats verbs (`analyzeData(covLoessParam / covGroupsParam)`)
#
# Accumulators, returned as vectors indexed by position along `axis` in the
# current view (n_genes long for "feat", n_cells for "cell"):
#
#   sum      sum(value)                       every stored entry
#   sumsq    sum(value^2)                     every stored entry
#   nnz      count(value > thr)               detection count
#   sum_det  sum(value where value > thr)     detected entries only
#
# `inclusive = TRUE` switches those two to `>=`, which is what filtering
# means by a threshold.
#
# `thr` is a detection predicate: it selects which entries COUNT, and never
# reduces a magnitude that participates. So `sum` / `sumsq` are unconditional
# and only `nnz` / `sum_det` see it. adr/0009.
#
# One query, unions included. `storeRead()` on a union already opens every
# substore as a single Dataset and composes their per-substore subset filters
# and the union's @ops into one plan, so there is nothing to iterate: Acero
# schedules all fragments together instead of running N plans in sequence.
#
# `source_id` joins the group key, which is what makes that safe. On the cell
# axis it is mandatory -- `row_id` restarts per substore, so grouping on it
# alone would merge cells from different samples. On the feature axis it keeps
# the identifier remap *after* the aggregate: each group knows its substore,
# so positions resolve against that substore's on-disk index without a join
# over the full stream (and without depending on substores agreeing about
# what a given `col_id` means).
#
# Two execution shapes:
#   * `@post_ops` empty -- the whole aggregate is pushed into Acero and only
#     the grouped result crosses into R.
#   * otherwise -- collect, apply the R-side chain, aggregate with data.table.
# Both hand off to the same join-and-fold tail.

# Build the `(source_id, key_id) -> pos` map for one axis. `key_id` is the
# on-disk id, `pos` the position in the view's axis.
#
# Slicing state lives entirely on the substores: `[` on a union pushes both
# axes down and rebuilds the parent, which carries no `@cell_idx`/`@gene_idx`
# of its own. So each substore's index vector plus (for cells) its offset into
# the union axis is the whole mapping.
.stream_expr_accum <- function(pe,
    axis = c("feat", "cell"),
    thr = 0,
    stats = c("sum", "sumsq", "nnz", "sum_det"),
    inclusive = FALSE) {

    if (!inherits(pe, "parquetExprBase")) {
        stop("[.stream_expr_accum] pe must be a parquetExprBase.",
             call. = FALSE)
    }
    axis  <- match.arg(axis)
    stats <- match.arg(stats, several.ok = TRUE)
    thr   <- as.numeric(thr)

    value <- pos <- NULL   # NSE bindings

    key   <- if (identical(axis, "feat")) "col_id" else "row_id"
    n_out <- as.integer(
        if (identical(axis, "feat")) pe@n_genes else pe@n_cells
    )

    out <- lapply(stats, function(nm) {
        if (identical(nm, "nnz")) integer(n_out) else numeric(n_out)
    })
    names(out) <- stats

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
    grp      <- if (is_union) c("source_id", key) else key
    join_by  <- if (is_union) c("source_id", "key_id") else "key_id"

    if (length(pe@post_ops) == 0L) {
        # Acero path: one plan over the whole store (union included), streamed
        # by the hash aggregate. Only the grouped result crosses into R.
        agg <- storeRead(pe, output = "query") |>
            dplyr::group_by(!!!rlang::syms(grp)) |>
            dplyr::summarise(!!!aggr_exprs) |>
            dplyr::collect() |>
            data.table::as.data.table()
    } else {
        # R path: the chain has to run on materialized rows, so this one is
        # chunked rather than collected whole.
        agg <- .pe_accum_chunked_dt(pe, grp = grp, aggr_exprs = aggr_exprs)
    }
    if (is.null(agg) || nrow(agg) == 0L) return(out)

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
    if (nrow(agg) == 0L) return(out)

    # Fold by position: real summation on the feature axis (a feature spans
    # substores), a no-op on the cell axis (cell axes are disjoint, so every
    # group is a singleton).
    folded <- agg[, lapply(.SD, sum), by = pos, .SDcols = stats]
    for (nm in stats) {
        out[[nm]][folded$pos] <- if (identical(nm, "nnz")) {
            as.integer(folded[[nm]])
        } else {
            as.numeric(folded[[nm]])
        }
    }

    out
}

# Window for the chunked R-side pass. Not `.exprbase_chunk_size()`: this pass
# collects a triplet frame -- row_id + col_id + value + source_id, plus the
# arrow collect buffer and the data.table copy -- measured at 47-54 B/nonzero
# across three shapes, so it asks for 48 and splits ~4x sooner. `k = 0`
# because there is no PCA sketch to reserve against.
.pe_accum_chunk_size <- function(pe) {
    .pe_window_cells(pe, bytes_per_nz = 48, k = 0L)
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
# Windows are cell ranges, not row counts, for two reasons -- a contiguous
# `row_id` range is the gapless case in `.pe_axis_pred()` so it prunes row
# groups, and `[` slices `@post_ops` down to the window's cells so the chain
# applies with exactly the state that window needs.
#
# Unions iterate substores here, unlike the Acero branch. A global cell range
# is not a contiguous `row_id` range across a union (row_id restarts per
# substore), so the window has to be taken per substore.
.pe_accum_chunked_dt <- function(pe, grp, aggr_exprs) {
    chunk_size <- .pe_accum_chunk_size(pe)
    dt_call    <- as.call(c(quote(list), aggr_exprs))
    is_union   <- inherits(pe, "unionParquetExprStore")

    parts <- list()
    for (sub_entry in .exprbase_substores(pe)) {
        sub <- sub_entry$store
        if (is_union) {
            sub <- .exprbase_inject_parent_ops(sub, pe@ops, pe@post_ops)
        }
        n_sub <- as.integer(sub@n_cells)
        cs <- 1L
        while (cs <= n_sub) {
            ce <- min(cs + chunk_size - 1L, n_sub)
            w  <- sub[, cs:ce]
            df <- data.table::as.data.table(
                dplyr::collect(storeRead(w, output = "query")))
            if (nrow(df) > 0L) {
                df <- .pe_apply_post_ops_df(df, w@post_ops)
                parts[[length(parts) + 1L]] <- df[, eval(dt_call), by = grp]
            }
            cs <- ce + 1L
        }
    }
    if (length(parts) == 0L) return(NULL)
    data.table::rbindlist(parts)
}

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
        .stream_pearson_resid_var(x, size_factors = param$size_factors)
    }
)

## internals --- hvf ####

# Per-gene variance of analytic (Poisson) Pearson residuals, streamed.
#
#   mu_ij = g_i * c_j / T        z_ij = (x_ij - mu_ij) / sqrt(mu_ij)
#
# with g_i the gene total, c_j the cell total (or a supplied size factor),
# and T the grand total.  Residuals are dense, but the all-zero block folds
# into closed form: treat every cell as zero, then correct only the stored
# nonzeros.
#
#   zeros:  sum_j z_ij^2 = sum_j mu_ij = g_i
#           sum_j z_ij   = -sqrt(g_i / T) * S,      S = sum_j sqrt(c_j)
#   stored nonzero (x > 0):
#           dz  = x / sqrt(mu_ij)
#           dz2 = x^2 / mu_ij - 2x
#
# So the whole statistic needs the gene totals, the cell totals, two scalars,
# and ONE joined pass over the existing triplet stream -- all pushed into
# Acero.  Validated against a dense reference to 2.7e-14; the same algebra
# runs on a dgCMatrix via @i/@p/@x if an in-memory version is ever wanted.
#
# @post_ops is deliberately NOT applied: the values must be counts.

.stream_pearson_resid_var <- function(pe, size_factors = NULL) {
    # NSE bindings
    row_id <- col_id <- value <- g <- cc <- mu <- dz <- dz2 <- NULL
    sum_dz <- sum_dz2 <- sum_z <- sum_z2 <- var <- feats <- NULL

    n_genes <- as.integer(pe@n_genes)

    # Pass 1: per-gene and per-cell totals (both pushed-down aggregates).
    gt <- storeRead(pe, output = "query") |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(g = sum(value, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()
    if (nrow(gt) == 0L) {
        return(data.table::data.table(feats = pe@feat_ids,
                                      var = numeric(n_genes)))
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
                                      var = numeric(n_genes)))
    }
    Tt <- sum(gt$g)
    S  <- sum(sqrt(cvals))

    # Pass 2: per-nonzero corrections, joined and aggregated in Acero.
    gt_a <- arrow::as_arrow_table(data.frame(
        col_id = as.integer(gt$col_id), g = as.numeric(gt$g)))
    ct_a <- arrow::as_arrow_table(data.frame(
        row_id = as.integer(ct$row_id), cc = as.numeric(ct$cc)))
    corr <- storeRead(pe, output = "query") |>
        dplyr::left_join(gt_a, by = "col_id") |>
        dplyr::left_join(ct_a, by = "row_id") |>
        dplyr::mutate(mu = g * cc / !!Tt) |>
        dplyr::mutate(dz  = value / sqrt(mu),
                      dz2 = value * value / mu - 2 * value) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(sum_dz  = sum(dz,  na.rm = TRUE),
                         sum_dz2 = sum(dz2, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()

    res <- merge(gt, corr, by = "col_id", all.x = TRUE)
    res[is.na(sum_dz),  sum_dz  := 0]
    res[is.na(sum_dz2), sum_dz2 := 0]
    res[, sum_z  := -sqrt(g / Tt) * S + sum_dz]
    res[, sum_z2 := g + sum_dz2]
    res[, var := (sum_z2 - sum_z * sum_z / n_eff) / (n_eff - 1L)]
    # Genes with no counts have undefined residuals; report 0 variance, which
    # is what rowVars() gives for an all-zero row in the in-memory path.
    res[g <= 0, var := 0]

    out <- numeric(n_genes)
    idx <- .pe_remap_col(res$col_id, pe)
    keep <- !is.na(idx)
    out[idx[keep]] <- as.numeric(res$var[keep])

    dt <- data.table::data.table(feats = pe@feat_ids, var = out)
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