#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-normalize ####
# Streaming JIT normalization for parquetExprBase. Plugs into Giotto's
# existing processData(x, param) dispatch via two setMethod calls on the
# shared `parquetExprBase` virtual:
#
#   processData(parquetExprBase, libraryNormParam)
#       -> appends a `norm_libsize_log` op (log = FALSE) to x@post_ops
#   processData(parquetExprBase, logNormParam)
#       -> if a `norm_libsize_log` op is already on @post_ops, flips its
#          `log = TRUE` flag in place (libsize+log fuse into one op record;
#          the R-side executor applies scale then log in a single pass).
#          Otherwise errors — log-only on raw counts isn't a documented
#          Giotto streaming path.
#
# Single (`parquetExprStore`) and union (`unionParquetExprStore`) collapse
# to one algorithm via `.exprbase_substores()`: each substore's per-cell
# libsize is computed via `.stream_colsums` and its slice of the
# `(source_id, orig_row_id, scalef)` lookup table is appended; the final
# table spans all substores in one op record. For a single store the loop
# runs once.
#
# Neither method rewrites the Parquet file. The recipe lives as a pure-data
# record on @post_ops and is applied R-side after collect (see
# .pe_apply_post_op_norm_libsize_log_df / _mat). The recipe survives
# saveRDS / load cycles without special handling — no closures.
#
# zscoreScaleParam is intentionally NOT implemented for parquetExprBase:
# per-cell / per-gene centering+scaling densifies the sparse matrix and
# breaks the O(N*k) streaming guarantee. normalizeGiotto already errors
# upstream when scale_cells / scale_feats = TRUE on a streaming backend.

# Build a per-substore `(source_id, orig_row_id, scalef)` slice for the
# `norm_libsize_log` op record. `scalef` is the per-cell scale factor
# vector in the SAME positional order as `store@cell_idx` (or
# 1..n_cells if no subset). Both single and union paths concatenate the
# per-substore slices into one composite-keyed table — row_id restarts
# per substore so source_id is the disambiguator.
.pe_norm_libsize_scalef_slice <- function(store, scalef) {
    orig_row_id <- if (length(store@cell_idx) > 0L) {
        as.integer(store@cell_idx)
    } else {
        seq_len(as.integer(store@n_cells))
    }
    data.table::data.table(
        source_id   = rep_len(as.character(store@uid), length(orig_row_id)),
        orig_row_id = orig_row_id,
        scalef      = as.numeric(scalef)
    )
}

# ---- libraryNormParam: compute & store JIT scale factors -------------------

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "libraryNormParam"),
    function(x, param, ...) {
        scalefactor <- param$scalefactor %null% 6e3

        subs <- .exprbase_substores(x)
        slices <- lapply(subs, function(sub) {
            store <- sub$store
            libsizes <- .stream_colsums(store)
            libsizes[libsizes == 0] <- 1   # guard div-by-zero
            scalef <- as.numeric(scalefactor) / libsizes
            .pe_norm_libsize_scalef_slice(store, scalef)
        })
        scalef_dt <- data.table::rbindlist(slices)

        # If a libsize-log op already exists on @post_ops (re-running
        # normalize), preserve its log flag and base. Otherwise default
        # to log=FALSE.
        existing <- .pe_find_op_type(x@post_ops, "norm_libsize_log")
        log_flag <- if (is.na(existing)) FALSE else
                    isTRUE(x@post_ops[[existing]]$log)
        log_base <- if (is.na(existing)) 2 else x@post_ops[[existing]]$base
        new_op <- list(
            type   = "norm_libsize_log",
            scalef = scalef_dt,
            log    = log_flag,
            base   = log_base
        )
        if (is.na(existing)) {
            .pe_push_op(x, new_op, phase = "post")
        } else {
            x@post_ops[[existing]] <- new_op
            x
        }
    }
)

# ---- logNormParam: set JIT log flag ----------------------------------------

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "logNormParam"),
    function(x, param, ...) {
        base   <- param$base   %null% 2
        offset <- param$offset %null% 1

        if (!isTRUE(offset == 1)) {
            stop("[processData(parquetExprBase, logNormParam)] ",
                 "offset != 1 is not supported for streaming because it ",
                 "would densify the sparse representation. Use offset = 1 ",
                 "(log1p) to preserve sparsity.", call. = FALSE)
        }

        existing <- .pe_find_op_type(x@post_ops, "norm_libsize_log")
        if (is.na(existing)) {
            stop("[processData(parquetExprBase, logNormParam)] no ",
                 "library-size normalization op present. Run ",
                 "processData(libraryNormParam) first; log-only on raw ",
                 "counts is not a supported streaming path.",
                 call. = FALSE)
        }
        # Fuse log flag onto the existing libsize op.
        x@post_ops[[existing]]$log  <- TRUE
        x@post_ops[[existing]]$base <- as.numeric(base)
        x
    }
)


# ---- internal: streaming colSums for parquetExprStore ----------------------
# Per-substore helper. Returns a full-length numeric vector of length
# `store@n_cells` (positions for cells with no non-zero entries are 0;
# callers guard with `[libsizes == 0] <- 1` as needed). Used by the
# union path via the substore iterator.

.stream_colsums <- function(pe) {
    if (!inherits(pe, "parquetExprStore"))
        stop("[.stream_colsums] pe must be a parquetExprStore.")

    n_cells <- as.integer(pe@n_cells)
    row_id <- value <- s <- NULL  # NSE bindings

    ds <- storeRead(pe)
    agg <- ds |>
        dplyr::group_by(row_id) |>
        dplyr::summarise(s = sum(value, na.rm = TRUE)) |>
        dplyr::collect() |>
        data.table::as.data.table()

    cs <- numeric(n_cells)
    if (nrow(agg) > 0L) {
        idx <- .pe_remap_row(agg$row_id, pe)
        keep <- !is.na(idx)
        cs[idx[keep]] <- as.numeric(agg$s[keep])
    }
    cs
}


# ---- internal: single-pass streaming per-cell stats ------------------------
# Pulls sum, sumsq, nnz, min_nz, max_nz from a single grouped aggregation
# (one Arrow scan, all stats from the same record-batch traversal — the
# BPCells pattern). Derives mean / var / sd in R using n_genes as the
# denominator so implicit zeros are counted.
#
# Returns a list of n_cells-long vectors indexed by cell position in the
# current view (gene_idx / cell_idx subsetting is applied via storeRead and
# .pe_remap_row, matching .stream_colsums).
#
#   sum, sumsq      doubles, 0 for cells with no non-zero entries
#   nnz             integer, count of non-zero entries per cell
#   min_nz, max_nz  doubles, min/max of *non-zero* values; NA when nnz == 0
#   mean            sum / n_genes  (zeros counted)
#   var             sample variance, (sumsq - sum^2 / n_genes) / (n_genes - 1)
#                   with a small-negative guard for cancellation near zero
#   sd              sqrt(var)

.stream_colstats <- function(pe) {
    if (!inherits(pe, "parquetExprStore"))
        stop("[.stream_colstats] pe must be a parquetExprStore.")

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)
    row_id <- value <- s <- s2 <- nnz <- vmin <- vmax <- NULL  # NSE

    ds <- storeRead(pe)
    agg <- ds |>
        dplyr::group_by(row_id) |>
        dplyr::summarise(
            s    = sum(value, na.rm = TRUE),
            s2   = sum(value * value, na.rm = TRUE),
            nnz  = dplyr::n(),
            vmin = min(value, na.rm = TRUE),
            vmax = max(value, na.rm = TRUE)
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    sum_v    <- numeric(n_cells)
    sumsq_v  <- numeric(n_cells)
    nnz_v    <- integer(n_cells)
    min_nz_v <- rep(NA_real_, n_cells)
    max_nz_v <- rep(NA_real_, n_cells)

    if (nrow(agg) > 0L) {
        idx <- .pe_remap_row(agg$row_id, pe)
        keep <- !is.na(idx)
        pos <- idx[keep]
        sum_v[pos]    <- as.numeric(agg$s[keep])
        sumsq_v[pos]  <- as.numeric(agg$s2[keep])
        nnz_v[pos]    <- as.integer(agg$nnz[keep])
        min_nz_v[pos] <- as.numeric(agg$vmin[keep])
        max_nz_v[pos] <- as.numeric(agg$vmax[keep])
    }

    mean_v <- if (n_genes > 0L) sum_v / n_genes else numeric(n_cells)
    if (n_genes > 1L) {
        var_v <- (sumsq_v - sum_v * sum_v / n_genes) / (n_genes - 1L)
        # cancellation can push true-zero variance slightly negative
        var_v[var_v < 0] <- 0
    } else {
        var_v <- numeric(n_cells)
    }
    sd_v <- sqrt(var_v)

    list(
        sum    = sum_v,
        sumsq  = sumsq_v,
        nnz    = nnz_v,
        min_nz = min_nz_v,
        max_nz = max_nz_v,
        mean   = mean_v,
        var    = var_v,
        sd     = sd_v
    )
}
