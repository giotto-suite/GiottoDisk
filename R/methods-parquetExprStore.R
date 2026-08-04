#' @include class-parquetExprStore.R
NULL

# storeRead ####
# Two-phase op execution:
#   @ops       arrow-lazy: folded into the composed query via @read_fun
#              wrap + .pe_apply_ops before materialization.
#   @post_ops  R-side post-collect: applied AFTER the arrow query has been
#              collected, mutating `value` in place. Only applied for
#              output modes that materialize (tibble / data.table /
#              dgcmatrix). `query` and `duckdb` outputs return the
#              lazy-only view — callers using those outputs are
#              responsible for applying @post_ops themselves if they
#              want fully-baked values.
#
# Index remapping is a SEPARATE bake step from @post_ops, and only the
# dgcmatrix output currently performs it:
#   dgcmatrix  row_id / col_id are remapped to subset positions (via
#              .pe_remap_row / .pe_remap_col, always AFTER @post_ops --
#              the executors key on on-disk ids) and surfaced as
#              @feat_ids / @cell_ids dimnames. Rows whose ids fall
#              outside the current subset are dropped.
#   tibble     row_id / col_id come back as raw ON-DISK ids, not yet
#   data.table remapped, and no unmappable-row filtering is applied.
#              Callers needing subset positions call .pe_remap_row /
#              .pe_remap_col themselves -- normally AFTER aggregating,
#              since remap is a per-row match() against @cell_idx /
#              @gene_idx and is far cheaper on an 18k-row aggregate than
#              on a full triplet stream.
#   query      raw on-disk ids (no materialization, so no remap).
#   duckdb
#
# The tabular outputs are expected to converge on dgcmatrix's remapped
# coordinates; until they do, treat their id columns as on-disk.

# ---- shared axis predicates -------------------------------------------------
#
# `cell_idx` / `gene_idx` narrowing is expressed as an arrow predicate. The
# naive form, `%in% <large int vec>`, does not row-group-prune and pays a hash
# probe per returned row; `row_id >= a & row_id <= b` uses the parquet
# row-group min/max statistics to skip untouched groups entirely.
#
# Strategy per axis:
#   1. gapless (unique count == range span) -> range predicate ALONE, which
#      exactly matches the kept set and beats a hash probe per row.
#   2. gaps, few dropped -> range + `!(x %in% dropped)`. The range is REQUIRED
#      for correctness here, or rows outside [lo, hi] pass the negated test.
#   3. gaps, few kept -> `x %in% kept` ALONE. That is already correct and
#      complete, and the range is deliberately NOT emitted: it can only help
#      via row-group pruning, and on the gene axis it never does -- rows are
#      written in cell batches, so every row group spans most col_ids.
#      Measured on Atera (2k of 18k genes, HVG-ranked): stacking a range onto
#      `is_in(kept)` cost +6.4% for zero pruning. On the cell axis a range
#      over a contiguous chunk was 3.9x faster than the `%in%` form.
#
# Bounds come from min/max, never the first/last element: `idx` is NOT
# guaranteed sorted (`feats_to_use` may be HVG-rank ordered), and these are
# set predicates -- row order is restored downstream by `.pe_remap_row` /
# `.pe_remap_col`. Gap detection works on unique values so duplicated entries
# cannot make `n == span` accidentally true while gaps remain, and `dropped`
# is only materialized when the anti-set branch actually wins.

#' @keywords internal
#' @noRd
.pe_axis_pred <- function(idx) {
    if (length(idx) == 0L) return(NULL)
    u  <- unique(as.integer(idx))
    lo <- min(u)
    hi <- max(u)
    span <- hi - lo + 1L
    gapless  <- (length(u) == span)
    use_anti <- !gapless && (span - length(u)) < length(u)
    list(lo = lo, hi = hi,
         gapless   = gapless,
         use_anti  = use_anti,
         use_range = gapless || use_anti,
         dropped   = if (use_anti) setdiff(lo:hi, u) else integer(0),
         kept      = u)
}

# Emit a plan as quoted predicate expressions over column `col`, with values
# inlined as literals. One emitter for both consumers: `storeRead` chains them
# onto a lazy query, `.union_substore_filter_expr()` ANDs them into a
# per-substore clause. Returns list() when there is nothing to filter.

#' @keywords internal
#' @noRd
.pe_axis_pred_exprs <- function(plan, col) {
    if (is.null(plan)) return(list())
    sym <- as.name(col)
    out <- list()
    if (plan$use_range) {
        out[[length(out) + 1L]] <-
            bquote(.(sym) >= .(plan$lo) & .(sym) <= .(plan$hi))
    }
    if (!plan$gapless) {
        out[[length(out) + 1L]] <- if (plan$use_anti) {
            bquote(!(.(sym) %in% .(plan$dropped)))
        } else {
            bquote(.(sym) %in% .(plan$kept))
        }
    }
    out
}


# * pestore ####

#' @rdname storeRead
#' @param max_rows,max_cols integer or `Inf` (optional). Set a dimension
#'   guard for `output = "dgcmatrix"`. A materialized sparseMatrix must
#'   have at least one axis narrowed to within the cap — both exceeding
#'   errors. Defaults via `getOption("giottodisk.dgc_max_rows", 100L)` /
#'   `getOption("giottodisk.dgc_max_cols", 100L)`. Asymmetric so
#'   "narrow slice" usage (100 × Inf or Inf × 100) is allowed; the
#'   intent is `[`-subset along one axis before materializing. Pass
#'   `Inf` to disable the cap on an axis (`Inf`/`Inf` disables the
#'   guard entirely).
#' @section Index coordinates by output mode:
#' `parquetExprStore` / `unionParquetExprStore` outputs differ on *two*
#' independent axes, and switching `output` changes both.
#'
#' Beyond whether `@post_ops` are applied (see Details), the `row_id` /
#' `col_id` columns are not in the same coordinate system across modes:
#'
#' * `output = "dgcmatrix"` returns values indexed by **subset position**,
#'   with `@feat_ids` / `@cell_ids` attached as dimnames. Rows whose ids
#'   fall outside the store's current `[`-subset are dropped.
#' * `output = "tibble"` and `output = "data.table"` return the raw
#'   **on-disk** `row_id` / `col_id`, without remapping and without
#'   unmappable-row filtering, even though `@post_ops` *are* applied. A
#'   caller that needs subset positions applies `.pe_remap_row()` /
#'   `.pe_remap_col()` itself, and is usually better off doing so after
#'   aggregating -- remapping is a per-row `match()` against `@cell_idx` /
#'   `@gene_idx`, so it is far cheaper on a per-gene aggregate than on a
#'   full triplet stream.
#' * `output = "query"` and `output = "duckdb"` never materialize, so ids
#'   are likewise raw on-disk values.
#'
#' The practical consequence: moving from `"dgcmatrix"` to a tabular
#' output does not error, but `row_id` / `col_id` silently change meaning.
#' The tabular modes are expected to converge on the remapped coordinates
#' in a later release; until then, treat their id columns as on-disk.
#' @export
setMethod("storeRead", signature("parquetExprStore"), function(store,
    output = c("query", "tibble", "duckdb", "dgcmatrix"),
    max_rows = NULL, max_cols = NULL, ...) {
    if (is.character(output)) output <- match.arg(output)

    # Phase 1: wrap @read_fun with subset filters + arrow-side @ops
    has_subset <- length(store@cell_idx) > 0L || length(store@gene_idx) > 0L
    has_ops    <- length(store@ops) > 0L
    if (has_subset || has_ops) {
        orig_rf <- store@read_fun
        ci  <- store@cell_idx
        gi  <- store@gene_idx
        ops <- store@ops
        # Axis predicates are built by `.pe_axis_pred()` / emitted by
        # `.pe_axis_pred_exprs()`, shared with the union read path so both
        # get the same row-group pruning. See those functions for the
        # strategy and the measurements behind it.
        preds <- c(
            .pe_axis_pred_exprs(.pe_axis_pred(ci), "row_id"),
            .pe_axis_pred_exprs(.pe_axis_pred(gi), "col_id")
        )

        store@read_fun <- function(x, ...) {
            ds <- orig_rf(x, ...)
            for (p in preds) ds <- dplyr::filter(ds, !!p)
            if (length(ops) > 0L) ds <- .pe_apply_ops(ds, ops)
            ds
        }
    }

    # dgcmatrix intercept: collect + apply @post_ops + build sparseMatrix
    if (identical(output, "dgcmatrix")) {
        atab <- callNextMethod(store = store, output = "query", ...)
        return(.pe_to_dgcmatrix(store, atab,
            max_rows = max_rows, max_cols = max_cols))
    }

    # query / duckdb: return lazy without applying @post_ops (documented)
    if (output %in% c("query", "duckdb")) {
        return(callNextMethod(store = store, output = output, ...))
    }

    # tibble path: collect + apply @post_ops R-side. @post_ops mutates
    # via data.table `:=`, so convert to data.table for the mutation, then
    # back to tibble to preserve the "tibble" output contract.
    df <- data.table::as.data.table(
        callNextMethod(store = store, output = "tibble", ...))
    df <- .pe_apply_post_ops_df(df, store@post_ops)
    dplyr::as_tibble(df)
})

# Build a Bioconductor-convention sparseMatrix (rows = genes, cols = cells)
# from the lazy arrow triplet query. The `x` slot is filled from `value`,
# which @post_ops has already mutated in place -- there is no separate
# normalized column. dimnames = (feat_ids, cell_ids) reflect any pending
# subset state since `[` already narrows those slots.
#
# Asymmetric dimension guard: a materialized matrix is intended as a
# "slice" of the dataset. At least one axis must be within
# max_rows / max_cols; both unbounded is rejected with a clear hint to
# `[`-subset first. Defaults from options
# `giottodisk.dgc_max_rows` / `giottodisk.dgc_max_cols` (100 each).
#
# TODO ScaledMatrix: when zscoreScaleParam / per-gene centering arrives
# (cell-population op kind), the dgCMatrix here can be wrapped in
# ScaledMatrix::ScaledMatrix(., center = gene_means, scale = gene_sds)
# so per-gene centering applies lazily on top of the sparse libsize-log
# values without densifying. The op record's gene-axis lookup table
# would supply `center` / `scale` directly.

.pe_check_dgc_dims <- function(n_rows, n_cols, max_rows, max_cols) {
    if (is.null(max_rows)) max_rows <-
        getOption("giottodisk.dgc_max_rows", 100L)
    if (is.null(max_cols)) max_cols <-
        getOption("giottodisk.dgc_max_cols", 100L)
    # Accept Inf as "no cap on this axis"; otherwise require a whole number.
    .as_cap <- function(x, nm) {
        if (is.numeric(x) && length(x) == 1L && !is.na(x) &&
            (is.infinite(x) || x %% 1 == 0)) {
            return(if (is.infinite(x)) x else as.integer(x))
        }
        stop("`", nm, "` must be a whole number or `Inf`.", call. = FALSE)
    }
    max_rows <- .as_cap(max_rows, "max_rows")
    max_cols <- .as_cap(max_cols, "max_cols")
    # Asymmetric guard: at least one axis must be ≤ its cap.
    if (n_rows > max_rows && n_cols > max_cols) {
        stop("[storeRead] dgcmatrix output would materialize ",
             format(n_rows, big.mark = ","), " genes x ",
             format(n_cols, big.mark = ","),
             " cells, with neither axis below the cap (max_rows = ",
             max_rows, ", max_cols = ", max_cols, "). ",
             "`[`-subset along one axis first, or override via ",
             "`max_rows` / `max_cols` args or options ",
             "`giottodisk.dgc_max_rows` / `giottodisk.dgc_max_cols`.",
             call. = FALSE)
    }
    invisible(NULL)
}

.pe_to_dgcmatrix <- function(store, atab,
    max_rows = NULL,
    max_cols = NULL) {
    n_rows <- as.integer(store@n_genes)
    n_cols <- as.integer(store@n_cells)
    .pe_check_dgc_dims(n_rows, n_cols, max_rows, max_cols)

    # `atab` is the lazy query with @ops + @cell_idx / @gene_idx already
    # composed (built by callNextMethod(output = "query") in the caller).
    df <- data.table::as.data.table(dplyr::collect(atab))
    row_id <- col_id <- value <- NULL  # NSE

    # Apply @post_ops R-side; mutates df$value in place.
    df <- .pe_apply_post_ops_df(df, store@post_ops)

    # Map original parquet ids -> 1..n_rows / 1..n_cols positions.
    # row_id (cell axis) -> j (matrix col); col_id (gene axis) -> i (matrix row).
    i_pos <- .pe_remap_col(df$col_id, store)
    j_pos <- .pe_remap_row(df$row_id, store)
    keep  <- !is.na(i_pos) & !is.na(j_pos)
    if (!all(keep)) {
        i_pos <- i_pos[keep]; j_pos <- j_pos[keep]; df <- df[keep]
    }
    x_col <- df$value

    Matrix::sparseMatrix(
        i = i_pos,
        j = j_pos,
        x = as.double(x_col),
        dims = c(n_rows, n_cols),
        dimnames = list(store@feat_ids, store@cell_ids),
        repr = "C"
    )
}

# * upestore ####

# storeRead: fuses substores via Arrow's UnionDataset (purely virtual —
# no filesystem ops). Per-substore subset state (@cell_idx / @gene_idx)
# is preserved by wrapping each substore's read_fun the same way
# parquetExprStore::storeRead does. Output dispatch (query / tibble /
# duckdb) handled inline since unionParquetExprStore doesn't inherit
# queryableStore.

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("unionParquetExprStore"), function(store,
    fields = NULL,
    output = c("query", "tibble", "duckdb", "dgcmatrix"),
    callback = NULL,
    duckdb_params = list(),
    max_rows = NULL, max_cols = NULL,
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "duckdb",
        "dgcmatrix"))

    # Open each substore as a Dataset (no read_fun wrapping), then union
    # them via `arrow::open_dataset(list(...))`. Per-substore
    # `@cell_idx`/`@gene_idx` filters from `[`-subset are applied as a
    # single composite source_id-aware expression below -- previously we
    # wrapped each substore's read_fun with `dplyr::filter`, which
    # converted Dataset to arrow_dplyr_query and broke the union step
    # (open_dataset(list(...)) can't accept arrow_dplyr_query objects).
    substore_dsets <- lapply(store@stores, .store_simple_read)
    atab <- arrow::open_dataset(substore_dsets)

    filt_expr <- .union_substore_filter_expr(store@stores)
    if (!is.null(filt_expr)) {
        atab <- dplyr::filter(atab, !!filt_expr)
    }
    # Apply union-level @ops chain. Composite (source_id, orig_row_id)
    # cell-axis keys span all substores in a single join, identical to
    # the single-store path. No per-substore dispatch needed.
    if (length(store@ops) > 0L) {
        atab <- .pe_apply_ops(atab, store@ops)
    }

    if (!is.null(fields)) atab <- dplyr::select(atab, dplyr::all_of(fields))
    if (!is.null(callback)) atab <- callback(atab)
    switch(output,
        # query / duckdb: return lazy without applying @post_ops
        "query"  = atab,
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params),
        # tibble: collect + apply @post_ops R-side; return tibble
        "tibble" = {
            df <- data.table::as.data.table(dplyr::collect(atab))
            df <- .pe_apply_post_ops_df(df, store@post_ops)
            dplyr::as_tibble(df)
        },
        # dgcmatrix: collect + apply @post_ops + build sparseMatrix
        "dgcmatrix" = .union_pe_to_dgcmatrix(store, atab,
            max_rows = max_rows, max_cols = max_cols)
    )
})


# dgCMatrix materialization for unionParquetExprStore — same shape as the
# single-store path. The union's storeRead has already composed @ops +
# substore filters into the lazy query; collect and build sparseMatrix.
# Position mapping is done via union-level @feat_ids / @cell_ids, which
# reflect any pending [`-subset state.

.union_pe_to_dgcmatrix <- function(store, atab,
                                    max_rows = NULL, max_cols = NULL) {
    n_rows <- as.integer(store@n_genes)
    n_cols <- as.integer(store@n_cells)
    .pe_check_dgc_dims(n_rows, n_cols, max_rows, max_cols)

    df <- data.table::as.data.table(dplyr::collect(atab))
    source_id <- row_id <- col_id <- value <- NULL  # NSE

    # Apply @post_ops R-side; mutates df$value in place. Post-op union
    # payload (e.g., norm_libsize's scalef table) carries composite
    # (source_id, orig_row_id) keys — the R-side update-join in
    # .pe_apply_post_op_norm_libsize_df resolves across substores.
    df <- .pe_apply_post_ops_df(df, store@post_ops)

    # j: union-global cell position from (source_id, row_id). Build a
    # lookup per substore: orig_row_id within the substore's @cell_idx
    # (or 1..n_cells if no subset) -> union-global position via offsets.
    j_pos <- integer(nrow(df))
    offset <- 0L
    for (s in store@stores) {
        ci <- if (length(s@cell_idx) > 0L) s@cell_idx
              else seq_len(as.integer(s@n_cells))
        n_keep <- length(ci)
        in_sub <- which(df$source_id == s@uid)
        if (length(in_sub) > 0L) {
            j_local <- match(df$row_id[in_sub], ci)
            j_pos[in_sub] <- offset + j_local
        }
        offset <- offset + n_keep
    }
    # i: gene position. Substores have aligned feat_ids by union
    # invariant; the union's @gene_idx reflects any pending gene subset.
    # Use the first substore's mapping for col_id -> gene position (same
    # in every substore).
    rep_store <- store@stores[[1L]]
    gi <- if (length(rep_store@gene_idx) > 0L) rep_store@gene_idx
          else seq_len(as.integer(rep_store@n_genes))
    i_pos <- match(df$col_id, gi)

    keep <- !is.na(i_pos) & !is.na(j_pos)
    if (!all(keep)) {
        i_pos <- i_pos[keep]; j_pos <- j_pos[keep]; df <- df[keep]
    }
    x_col <- df$value

    Matrix::sparseMatrix(
        i = i_pos,
        j = j_pos,
        x = as.double(x_col),
        dims = c(n_rows, n_cols),
        dimnames = list(store@feat_ids, store@cell_ids),
        repr = "C"
    )
}

# Build the source_id-aware composite filter expression for a union read.
#
# Returns an unevaluated R call suitable for `dplyr::filter(atab, !!expr)`,
# or NULL if no substore has a `@cell_idx`/`@gene_idx` subset (no filter
# needed, all rows pass).
#
# Shape:
#   (source_id == uid_a & row_id %in% ci_a & col_id %in% gi_a) |
#   (source_id == uid_b & row_id %in% ci_b)                    |
#   ...                                                        |
#   (source_id %in% <no_subset_uids>)
#
# A substore without any subset contributes its uid to the trailing
# `source_id %in% <...>` OR'd clause, so all rows from that substore pass.
#
# Efficiency: arrow's expression engine compiles `%in%` to `is_in()` (hash
# set, O(1) per row) and applies parquet stats pushdown on `row_id` /
# `col_id` row-group min/max, so the engine skips entire row groups that
# can't contain a match. Cost scales with `length(@cell_idx) +
# length(@gene_idx)` per substore (inline IN list size). The natural
# ceiling is the cell/feature ID universe; for parquetExprStore the
# slot-resident char-vector limit (~10^8 cells) keeps IN lists tractable.
# Atlas-scale graphs (transcript edges, 10^9+) would need a sidecar
# parquet for the ID universe + join -- not applicable to expr stores.
.union_substore_filter_expr <- function(stores) {
    needs_filter <- vapply(stores, function(s) {
        length(s@cell_idx) > 0L || length(s@gene_idx) > 0L
    }, logical(1L))
    if (!any(needs_filter)) return(NULL)

    no_subset_uids <- vapply(stores[!needs_filter],
        function(s) s@uid, character(1L))

    clauses <- list()
    for (i in which(needs_filter)) {
        s <- stores[[i]]
        parts <- c(
            list(bquote(source_id == .(s@uid))),
            .pe_axis_pred_exprs(.pe_axis_pred(s@cell_idx), "row_id"),
            .pe_axis_pred_exprs(.pe_axis_pred(s@gene_idx), "col_id")
        )
        clause <- Reduce(function(a, b) bquote(.(a) & .(b)), parts)
        clauses[[length(clauses) + 1L]] <- clause
    }
    if (length(no_subset_uids) > 0L) {
        clauses[[length(clauses) + 1L]] <- bquote(source_id %in% .(no_subset_uids))
    }
    Reduce(function(a, b) bquote(.(a) | .(b)), clauses)
}

# storeWrite ####

# * pestore ####

# Direct parquet -> parquet path. Reads the input's lazy triplet stream
# (with `[`-subset filters already applied), remaps row_id / col_id from
# original-substore positions to local positions in the input's narrowed
# universe via a small in-mem arrow lookup, then streams to the output's
# `source_id=<new_uid>/` hive partition. No matrix materialization.

#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("parquetExprStore", "parquetExprStore"),
    function(store, data, ...) {
        .pestore_write_from_parquet_input(store, data)
    }
)

# * upestore ####

#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("parquetExprStore", "unionParquetExprStore"),
    function(store, data, ...) {
        .pestore_write_from_parquet_input(store, data)
    }
)

# Shared body for both signatures above. `data` is the input store; `store`
# is the fresh target parquetExprStore (path + uid).
.pestore_write_from_parquet_input <- function(store, data) {
    if (file.exists(store@path) && !dir.exists(store@path)) {
        stop("[storeWrite] output path exists as a file: ", store@path,
            "\n  pre-allocated store path must be a directory or absent.",
            call. = FALSE)
    }

    # Dispatch on `@post_ops`:
    #   * empty: existing fast arrow-lazy remap. Cell/gene narrowing is
    #     handled via left-join in `.pestore_remap_query`, no R-side pass.
    #   * non-empty: chunk-stream bake via the mat-shape @post_ops
    #     executor.  Each chunk becomes a sparseMatrix, @post_ops
    #     applies in place via positional scalef indexing (no arrow
    #     left-join), triplets are extracted and written to a shard.
    #     Row ids get remapped to the narrowed output cell axis.
    if (length(data@post_ops) > 0L) {
        .pestore_write_baked(store, data)
    } else {
        q <- .pestore_remap_query(data)
        .write_parquet(store, q)
    }

    # Slots copied from input -- already correctly narrowed by any
    # `[`-subset. `cell_idx`/`gene_idx` stay empty on the fresh output
    # (no subset queued).  `@post_ops` is also empty on the fresh output
    # since values are baked in.
    store@cell_ids <- data@cell_ids
    store@feat_ids <- data@feat_ids
    store@n_cells  <- as.numeric(data@n_cells)
    store@n_genes  <- as.numeric(data@n_genes)
    .pestore_finalize_chunk_size(store)
}


# Chunk-stream bake: iterate substores + chunks, read raw triplets, build
# a sparseMatrix band, apply `@post_ops` in place via the mat-shape
# executor (positional scalef — no arrow left-join), extract normalized
# triplets, remap row_id to the narrowed output cell axis, write a shard.
#
# Complexity: O(nnz) per chunk for build + apply + extract; global sort
# of shards is implicit because chunks are visited in cell order and each
# shard is sorted by (row_id, col_id) at write time.
.pestore_write_baked <- function(store, data) {
    n_out      <- as.integer(data@n_cells)
    P_out      <- as.integer(data@n_genes)
    chunk_size <- as.integer(.exprbase_chunk_size(data))
    hvg_idx    <- seq_len(P_out)   # identity on the narrowed feat axis
    post_ops   <- data@post_ops

    # Build sub_infos, one entry per substore.
    parent_ops <- if (inherits(data, "unionParquetExprStore")) data@ops
                  else list()
    # Both phase chains are transplanted in one step, so each substore is
    # self-sufficient for reading: `[` then slices @post_ops per chunk and
    # `storeRead` applies them.
    sub_infos <- lapply(.exprbase_substores(data), function(se) {
        sub <- .exprbase_inject_parent_ops(se$store, parent_ops, post_ops)
        # scalef_vecs is vestigial here for the same reason as on the PCA
        # paths: `[` slices @post_ops per chunk and `storeRead` applies them,
        # so no caller reads a pre-extracted per-cell vector. Kept empty so
        # `info` is one shape across readers.
        list(sub         = sub,
             offset      = as.integer(se$cell_offset),
             n_sub       = as.integer(sub@n_cells),
             hvg_orig    = .pe_orig_col(hvg_idx, sub),
             scalef_vecs = list())
    })

    partition_dir <- .idpath(store@path, store@uid)
    if (!dir.exists(partition_dir)) {
        dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)
    }

    part_idx <- 0L
    for (info in sub_infos) {
        n_sub  <- info$n_sub
        offset <- info$offset
        cs <- 1L
        while (cs <= n_sub) {
            ce <- min(cs + chunk_size - 1L, n_sub)
            M  <- .pe_read_chunk_sub(info, cs, ce, post_ops, P_out)
            if (!is.null(M) && length(M@x) > 0L) {
                # M is P_out × chunk_n (genes × cells, normalized), so the
                # CSC structure gives gene per nonzero via @i and cell per
                # nonzero by expanding the column pointers.
                col_pos    <- M@i + 1L                # gene pos, 1-based
                cell_local <- rep.int(seq_len(ncol(M)), diff(M@p))
                row_global <- offset + cs + cell_local - 1L
                dt <- data.table::data.table(
                    row_id = as.integer(row_global),
                    col_id = as.integer(col_pos),
                    value  = M@x
                )
                data.table::setorder(dt, row_id, col_id)
                .write_parquet_file(dt,
                    file.path(partition_dir,
                        sprintf("part-%d.parquet", part_idx)))
                part_idx <- part_idx + 1L
            }
            cs <- ce + 1L
        }
    }
    invisible(NULL)
}

# Build the lazy remapped triplet query for a (possibly subset) input.
# Returns an arrow_dplyr_query whose schema is (row_id, col_id, value),
# row_id / col_id renumbered to local positions in the input's narrowed
# `@cell_ids` / `@feat_ids` universe, sorted by (row_id, col_id).
.pestore_remap_query <- function(pe) {
    is_union <- inherits(pe, "unionParquetExprStore")

    # ---- row_id (cell) remap table ----
    # Use data.table for LUT construction: data.table::rbindlist over the
    # per-substore frames is materially faster than do.call(rbind, ...) when
    # the union has many substores. Same constructor style as the
    # (parquetExprStore, memoryMatrix) write path above.
    if (is_union) {
        n_per <- vapply(pe@stores, function(s) as.integer(s@n_cells),
            integer(1L))
        offsets <- c(0L, cumsum(n_per)[-length(n_per)])
        cell_remap_dt <- data.table::rbindlist(lapply(seq_along(pe@stores),
            function(k) {
                s <- pe@stores[[k]]
                local_orig <- if (length(s@cell_idx) > 0L) {
                    s@cell_idx
                } else {
                    seq_len(s@n_cells)
                }
                data.table::data.table(
                    source_id   = s@uid,
                    row_id_orig = as.integer(local_orig),
                    row_id_new  = as.integer(seq_along(local_orig) +
                        offsets[k])
                )
            }))
    } else {
        local_orig <- if (length(pe@cell_idx) > 0L) {
            pe@cell_idx
        } else {
            seq_len(pe@n_cells)
        }
        cell_remap_dt <- data.table::data.table(
            row_id_orig = as.integer(local_orig),
            row_id_new  = as.integer(seq_along(local_orig))
        )
    }
    cell_remap <- arrow::as_arrow_table(cell_remap_dt)

    # ---- col_id (gene) remap table ----
    # Gene-axis subset applies uniformly across union substores (the union's
    # `[` method calls `s[i, ]` on each substore), so any substore's
    # `@gene_idx` is representative.
    gene_idx <- if (is_union) pe@stores[[1L]]@gene_idx else pe@gene_idx
    do_gene_remap <- length(gene_idx) > 0L
    if (do_gene_remap) {
        gene_remap <- arrow::as_arrow_table(data.table::data.table(
            col_id_orig = as.integer(gene_idx),
            col_id_new  = as.integer(seq_along(gene_idx))
        ))
    }

    # ---- join + remap ----
    # `by` must be fully named -- arrow's dplyr join handler trips on the
    # mixed-named form `c("source_id", "row_id" = "row_id_orig")` because
    # the unnamed element parses with an empty name on the right side.
    q <- storeRead(pe)
    if (is_union) {
        q <- dplyr::left_join(q, cell_remap,
            by = c("source_id" = "source_id",
                   "row_id" = "row_id_orig"))
    } else {
        q <- dplyr::left_join(q, cell_remap,
            by = c("row_id" = "row_id_orig"))
    }
    q <- dplyr::mutate(q, row_id = row_id_new)
    q <- dplyr::select(q, -dplyr::any_of(c("row_id_new", "source_id")))
    if (do_gene_remap) {
        q <- dplyr::left_join(q, gene_remap,
            by = c("col_id" = "col_id_orig"))
        q <- dplyr::mutate(q, col_id = col_id_new)
        q <- dplyr::select(q, -col_id_new)
    }
    q <- dplyr::arrange(q, row_id, col_id)
    dplyr::select(q, row_id, col_id, value)
}

# from a dgCMatrix / Matrix / matrix.  Convenience path: useful for tests
# and small datasets that already live in memory.  The streaming Input
# classes (e.g. `mtxInput()`, `tenxH5Input()`) plus
# `storeWrite(parquetExprStore, exprInput)` are the production entry
# point for raw inputs — that path never materializes a dgCMatrix.

# * memoryMatrix ####

#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("parquetExprStore", "memoryMatrix"),
    function(store, data, ...) {
        # Coerce to dgCMatrix for uniform .summary access.
        if (!inherits(data, "dgCMatrix")) {
            data <- methods::as(data, "CsparseMatrix")
        }
        sm <- Matrix::summary(data)
        # In Giotto convention, expression matrices are gene x cell
        # (rows = genes, cols = cells), whereas the store's triplet schema is
        # cell-major: row_id = cell index, col_id = gene index (row_id is the
        # sort key, so cell-major is what makes row-group pruning work). Hence
        # the flip: row_id <- sm$j, col_id <- sm$i.
        dt <- data.table::data.table(
            row_id = as.integer(sm$j),
            col_id = as.integer(sm$i),
            value  = as.double(sm$x)
        )
        data.table::setorder(dt, row_id, col_id)

        # source_id=<uid>/ hive partition layout — shared with parquetStore
        # via .write_parquet (calls arrow::write_dataset, produces
        # part-N.parquet naming). A union store can hardlink substore
        # partition dirs without renaming or rewriting files.
        if (file.exists(store@path) && !dir.exists(store@path)) {
            unlink(store@path)
        }
        .write_parquet(store, dt)

        store@n_cells <- as.numeric(ncol(data))
        store@n_genes <- as.numeric(nrow(data))
        if (length(store@cell_ids) == 0L && !is.null(colnames(data)))
            store@cell_ids <- as.character(colnames(data))
        if (length(store@feat_ids) == 0L && !is.null(rownames(data)))
            store@feat_ids <- as.character(rownames(data))
        .pestore_finalize_chunk_size(store)
    }
)



# dim / nrow / ncol ####
# Bioconductor convention: expression matrices are gene x cell, so
# nrow = genes and ncol = cells.

#' @export
setMethod("nrow", "parquetExprStore", function(x) x@n_genes)

#' @export
setMethod("ncol", "parquetExprStore", function(x) x@n_cells)

#' @export
setMethod("dim", "parquetExprStore", function(x) c(x@n_genes, x@n_cells))

#' @export
setMethod("nrow", "unionParquetExprStore", function(x) x@n_genes)

#' @export
setMethod("ncol", "unionParquetExprStore", function(x) x@n_cells)

#' @export
setMethod("dim", "unionParquetExprStore",
    function(x) c(x@n_genes, x@n_cells)
)

# dimnames / rownames / colnames ####
# `rownames()` and `colnames()` in base R consult `dimnames()` first; defining
# `dimnames` here makes them work uniformly without separate methods.

#' @export
setMethod("dimnames", "unionParquetExprStore",
    function(x) list(x@feat_ids, x@cell_ids)
)

#' @export
setMethod("dimnames", "parquetExprStore",
    function(x) list(x@feat_ids, x@cell_ids)
)

# `rownames<-` and `colnames<-` fall back to `dimnames<-`. Define the
# setter so downstream Giotto code that does `rownames(x) <- ...` after
# normalize works transparently with our class.
#' @export
setMethod("dimnames<-",
    signature(x = "parquetExprStore", value = "list"),
    function(x, value) {
        if (length(value) >= 1L && !is.null(value[[1L]])) {
            x@feat_ids <- as.character(value[[1L]])
        }
        if (length(value) >= 2L && !is.null(value[[2L]])) {
            x@cell_ids <- as.character(value[[2L]])
        }
        x
    }
)

# show methods live in methods-show.R alongside the other store types.


# [ subset ####

# `pe[i, j]` returns a new parquetExprStore narrowed to the kept rows
# (genes) / columns (cells). The Parquet file on disk is unchanged;
# the @cell_idx / @gene_idx slots record the original-parquet positions
# of the kept entries so storeRead can filter via Arrow.
#
# Bioconductor convention: rows = genes, cols = cells.
# Supported index types: integer, logical, character, missing.

.resolve_subset_idx <- function(idx, all_ids, axis_name) {
    if (is.logical(idx)) {
        if (length(idx) != length(all_ids)) {
            stop("[parquetExprStore subset] logical ", axis_name,
                 " index length (", length(idx), ") != n (", length(all_ids),
                 ").", call. = FALSE)
        }
        return(which(idx))
    }
    if (is.character(idx)) {
        m <- match(idx, all_ids)
        if (anyNA(m)) {
            bad <- idx[is.na(m)]
            stop("[parquetExprStore subset] character ", axis_name,
                 " index has unknown IDs: ",
                 toString(head(bad, 5L)),
                 if (length(bad) > 5L) ", ..." else "",
                 call. = FALSE)
        }
        return(m)
    }
    if (is.numeric(idx)) {
        return(as.integer(idx))
    }
    stop("[parquetExprStore subset] unsupported ", axis_name, " index type: ",
         class(idx)[1L], call. = FALSE)
}

# Helpers used by streaming methods after storeRead(pe) returns rows
# whose row_id / col_id are still in the ORIGINAL parquet coordinate
# system. After collect(), call these to map to subset positions
# (1..n_cells / 1..n_genes of the current view).

#' @keywords internal
#' @noRd
.pe_remap_row <- function(orig_row_ids, pe) {
    if (length(pe@cell_idx) == 0L) return(as.integer(orig_row_ids))
    as.integer(match(orig_row_ids, pe@cell_idx))
}

#' @keywords internal
#' @noRd
.pe_remap_col <- function(orig_col_ids, pe) {
    if (length(pe@gene_idx) == 0L) return(as.integer(orig_col_ids))
    as.integer(match(orig_col_ids, pe@gene_idx))
}

# Translate a vector of subset positions to the original parquet col_ids --
# used when a method needs the on-disk gene indices rather than view positions.

#' @keywords internal
#' @noRd
.pe_orig_col <- function(subset_pos, pe) {
    if (length(pe@gene_idx) == 0L) return(as.integer(subset_pos))
    as.integer(pe@gene_idx[subset_pos])
}

#' @export
setMethod("[",
    signature(x = "parquetExprStore", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        # i = genes (rows); j = cells (cols)
        if (!missing(i)) {
            i_int <- .resolve_subset_idx(i, x@feat_ids, "row (gene)")
            new_gene_idx <- if (length(x@gene_idx) == 0L) {
                as.integer(i_int)
            } else {
                x@gene_idx[i_int]
            }
            new_feat_ids <- x@feat_ids[i_int]
            x@feat_ids <- new_feat_ids
            x@gene_idx <- as.integer(new_gene_idx)
            x@n_genes  <- as.numeric(length(x@feat_ids))
            # Slice gene-axis op tables on both phase chains. Current op
            # kinds have no gene-axis tables; this is a no-op for
            # the current kinds but generic for future ones.
            surviving_genes <- data.table::data.table(
                feat_id = as.character(new_feat_ids))
            if (length(x@ops) > 0L) {
                x@ops <- .pe_slice_ops(x@ops, axis = "gene",
                    surviving_keys = surviving_genes)
            }
            if (length(x@post_ops) > 0L) {
                x@post_ops <- .pe_slice_ops(x@post_ops, axis = "gene",
                    surviving_keys = surviving_genes)
            }
        }
        if (!missing(j)) {
            j_int <- .resolve_subset_idx(j, x@cell_ids, "col (cell)")
            new_cell_idx <- if (length(x@cell_idx) == 0L) {
                as.integer(j_int)
            } else {
                x@cell_idx[j_int]
            }
            x@cell_ids <- x@cell_ids[j_int]
            x@cell_idx <- as.integer(new_cell_idx)
            x@n_cells  <- as.numeric(length(x@cell_ids))
            # Slice cell-axis op tables on both phase chains. Composite
            # key (source_id, orig_row_id); source_id is constant (=
            # x@uid) for a single store.
            surviving_cells <- data.table::data.table(
                source_id   = rep_len(as.character(x@uid),
                                      length(new_cell_idx)),
                orig_row_id = as.integer(new_cell_idx))
            if (length(x@ops) > 0L) {
                x@ops <- .pe_slice_ops(x@ops, axis = "cell",
                    surviving_keys = surviving_cells)
            }
            if (length(x@post_ops) > 0L) {
                x@post_ops <- .pe_slice_ops(x@post_ops, axis = "cell",
                    surviving_keys = surviving_cells)
            }
        }
        x
    }
)

# i (genes) — applied uniformly to all substores (feat_ids are shared).
# j (cells) — mapped from union positions to per-substore positions via
# cumulative offsets; substores that get zero cells after the subset
# are dropped. Result is rebuilt through the constructor for invariant
# checks. Union-level @ops survive and get axis-sliced (cell axis uses
# the composite (source_id, orig_row_id) key spanning surviving
# substores; gene axis uses surviving feat_ids).

#' @export
setMethod("[",
    signature(x = "unionParquetExprStore", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        if (!missing(i)) {
            new_stores <- lapply(x@stores, function(s) s[i, ])
        } else {
            new_stores <- x@stores
        }
        if (!missing(j)) {
            j_int <- .resolve_subset_idx(j, x@cell_ids, "col (cell)")
            offsets <- c(0L, cumsum(vapply(new_stores,
                function(s) s@n_cells, numeric(1L))))
            kept <- list()
            for (k in seq_along(new_stores)) {
                lo <- as.integer(offsets[k]) + 1L
                hi <- as.integer(offsets[k + 1L])
                in_range <- j_int >= lo & j_int <= hi
                if (any(in_range)) {
                    local_j <- j_int[in_range] - as.integer(offsets[k])
                    kept[[length(kept) + 1L]] <- new_stores[[k]][, local_j]
                }
            }
            if (length(kept) == 0L) {
                stop("[unionParquetExprStore] cell subset selected no ",
                     "cells from any substore", call. = FALSE)
            }
            new_stores <- kept
        }
        new_union <- unionParquetExprStore(new_stores)

        # Inherit + slice union-level phase chains along the axes being
        # narrowed. Substore @ops / @post_ops stay empty by constraint;
        # only the union carries ops, and slicing applies on cell axis
        # with the composite (source_id, orig_row_id) key spanning
        # surviving substores.
        new_union@ops      <- x@ops
        new_union@post_ops <- x@post_ops
        has_any_ops <- length(new_union@ops) > 0L ||
                       length(new_union@post_ops) > 0L
        if (has_any_ops) {
            if (!missing(i)) {
                surviving_genes <- data.table::data.table(
                    feat_id = as.character(new_union@feat_ids))
                if (length(new_union@ops) > 0L)
                    new_union@ops <- .pe_slice_ops(new_union@ops,
                        axis = "gene", surviving_keys = surviving_genes)
                if (length(new_union@post_ops) > 0L)
                    new_union@post_ops <- .pe_slice_ops(new_union@post_ops,
                        axis = "gene", surviving_keys = surviving_genes)
            }
            if (!missing(j)) {
                surviving_cells <- data.table::rbindlist(lapply(
                    new_stores, function(s) {
                        ci <- if (length(s@cell_idx) > 0L) s@cell_idx
                              else seq_len(as.integer(s@n_cells))
                        data.table::data.table(
                            source_id   = rep_len(as.character(s@uid),
                                                  length(ci)),
                            orig_row_id = as.integer(ci))
                    }))
                if (length(new_union@ops) > 0L)
                    new_union@ops <- .pe_slice_ops(new_union@ops,
                        axis = "cell", surviving_keys = surviving_cells)
                if (length(new_union@post_ops) > 0L)
                    new_union@post_ops <- .pe_slice_ops(new_union@post_ops,
                        axis = "cell", surviving_keys = surviving_cells)
            }
        }
        new_union
    }
)

# cbind ####

# cbind2: pairwise combination producing a unionParquetExprStore. Higher
# arity (cbind(a, b, c, d)) lands here pairwise via base R's cbind/Matrix
# dispatch — left-fold builds a chain unionParquetExprStore(list(a, b)),
# then unionParquetExprStore(c(<existing union>@stores, list(c))).

#' @rdname unionParquetExprStore
#' @export
setMethod("cbind2",
    signature("parquetExprStore", "parquetExprStore"),
    function(x, y, ...) unionParquetExprStore(list(x, y))
)

#' @rdname unionParquetExprStore
#' @export
setMethod("cbind2",
    signature("unionParquetExprStore", "parquetExprStore"),
    function(x, y, ...) unionParquetExprStore(c(x@stores, list(y)))
)

#' @rdname unionParquetExprStore
#' @export
setMethod("cbind2",
    signature("parquetExprStore", "unionParquetExprStore"),
    function(x, y, ...) unionParquetExprStore(c(list(x), y@stores))
)

#' @rdname unionParquetExprStore
#' @export
setMethod("cbind2",
    signature("unionParquetExprStore", "unionParquetExprStore"),
    function(x, y, ...) unionParquetExprStore(c(x@stores, y@stores))
)
