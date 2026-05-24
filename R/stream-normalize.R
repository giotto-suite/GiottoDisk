#' @include class-parquetExprStore.R
NULL

# stream-normalize ####
# Streaming JIT normalization for parquetExprStore. Plugs into Giotto's
# existing processData(x, param) dispatch via two setMethod calls:
#
#   processData(parquetExprStore, libraryNormParam)
#       -> parquetExprStore with @params$norm$scale_factors set
#   processData(parquetExprStore, logNormParam)
#       -> parquetExprStore with @params$norm$log set
#
# Neither method rewrites the Parquet file. The recipe lives on the store
# and is applied on-the-fly by downstream streaming readers (sc_hvg,
# sc_pca dispatch in later Phase 2 steps). This preserves the JIT
# normalization design from scstream.
#
# zscoreScaleParam is intentionally NOT implemented for parquetExprStore:
# per-cell / per-gene centering+scaling densifies the sparse matrix and
# breaks the O(N*k) streaming guarantee. normalizeGiotto already errors
# upstream when scale_cells / scale_feats = TRUE on a streaming backend.

# ---- libraryNormParam: compute & store JIT scale factors -------------------

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprStore", param = "libraryNormParam"),
    function(x, param, ...) {
        scalefactor <- param$scalefactor %null% 6e3

        libsizes <- .stream_colsums(x)
        libsizes[libsizes == 0] <- 1   # guard against div-by-zero

        norm <- x@params$norm %null% list()
        norm$method        <- "library_size"
        norm$scalefactor   <- as.numeric(scalefactor)
        norm$scale_factors <- as.numeric(scalefactor) / libsizes

        x@params$norm <- norm
        x
    }
)

# ---- logNormParam: set JIT log flag ----------------------------------------

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprStore", param = "logNormParam"),
    function(x, param, ...) {
        base   <- param$base   %null% 2
        offset <- param$offset %null% 1

        if (!isTRUE(offset == 1)) {
            stop("[processData(parquetExprStore, logNormParam)] ",
                 "offset != 1 is not supported for streaming because it ",
                 "would densify the sparse representation. Use offset = 1 ",
                 "(log1p) to preserve sparsity.", call. = FALSE)
        }

        norm <- x@params$norm %null% list()
        norm$log    <- TRUE
        norm$base   <- as.numeric(base)
        norm$offset <- 1
        x@params$norm <- norm
        x
    }
)


# ---- internal: streaming colSums for parquetExprStore ----------------------

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
