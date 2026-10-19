#' @include pkg_imports.R class-parquetExprStore.R
NULL

# crossprod(W, x) = t(W) %*% x for an expression store x (features x cells)
# and a dense features x k matrix W: k weighted sums of features per cell, in
# one windowed pass over the cells. GiottoClass's createMetafeats() scores
# every module this way (mean and sum), so a store is read once rather than
# once per module.

#' @title Weighted feature sums of an expression store
#' @name crossprod-parquetExprBase
#' @description `crossprod(W, x)` with `x` an expression store (features x
#' cells) and `W` a numeric matrix with one row per feature of `x`: the
#' `ncol(W) x cells` matrix `t(W) %*% x`, computed in one pass over cell
#' windows.
#' @param x a numeric matrix, features x k
#' @param y a `parquetExprBase` store
#' @returns a dense `ncol(x) x cells` matrix, cells named
#' @export
setMethod("crossprod", signature(x = "matrix", y = "parquetExprBase"),
    function(x, y) {
        if (nrow(x) != length(y@feat_ids)) {
            stop("[crossprod] `x` needs one row per feature of the store (",
                length(y@feat_ids), ").", call. = FALSE)
        }
        out <- matrix(0, ncol(x), length(y@cell_ids),
            dimnames = list(colnames(x), y@cell_ids))
        # reading a window into a sparse matrix holds its arrow buffers, the
        # triplet frame and the matrix at once: measured ~140 bytes a value
        # at 10M cells, so windows are sized at 150
        cells <- .pe_window_cells(y, bytes_per_nz = 150, k = 0L)
        for (d in .pe_windows(y, cells)) {
            m <- storeRead(.pe_window_store(d), output = "dgcmatrix",
                max_rows = Inf, max_cols = Inf)
            if (!length(m@x)) next
            pos <- d$offset + (d$cs:d$ce)
            j <- pos[match(colnames(m), y@cell_ids[pos])]
            out[, j] <- as.matrix(Matrix::crossprod(
                x[match(rownames(m), y@feat_ids), , drop = FALSE], m))
        }
        out
    }
)
