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
