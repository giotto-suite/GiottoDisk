#' @include class-parquetExprStore.R
NULL

# storeRead ####

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetExprStore"), function(store, ...) {
    store@read_fun(store@path, ...)
})

# storeWrite ####
# from a dgCMatrix / Matrix / matrix.  Convenience path: useful for tests
# and small datasets that already live in memory.  The streaming converter
# `mtx_to_parquetExprStore()` is the production entry point for raw inputs
# (10x / Xenium MatrixMarket) — that path never materializes a dgCMatrix.

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
        # In Giotto convention, expression matrices are gene × cell
        # (rows = genes, cols = cells).  scstream's row_id = cell index,
        # col_id = gene index.  So row_id <- sm$j, col_id <- sm$i.
        dt <- data.table::data.table(
            row_id = as.integer(sm$j),
            col_id = as.integer(sm$i),
            value  = as.double(sm$x)
        )
        data.table::setorder(dt, row_id, col_id)
        arrow::write_parquet(dt, store@path)

        store@n_cells <- as.numeric(ncol(data))
        store@n_genes <- as.numeric(nrow(data))
        if (length(store@cell_ids) == 0L && !is.null(colnames(data)))
            store@cell_ids <- as.character(colnames(data))
        if (length(store@feat_ids) == 0L && !is.null(rownames(data)))
            store@feat_ids <- as.character(rownames(data))
        store
    }
)

# dim / nrow / ncol ####
# Bioconductor convention: expression matrices are gene × cell, so
# nrow = genes and ncol = cells.

#' @export
setMethod("nrow", "parquetExprStore", function(x) x@n_genes)

#' @export
setMethod("ncol", "parquetExprStore", function(x) x@n_cells)

#' @export
setMethod("dim", "parquetExprStore", function(x) c(x@n_genes, x@n_cells))

# dimnames / rownames / colnames ####
# `rownames()` and `colnames()` in base R consult `dimnames()` first; defining
# `dimnames` here makes them work uniformly without separate methods.

#' @export
setMethod("dimnames", "parquetExprStore",
    function(x) list(x@feat_ids, x@cell_ids)
)

# show ####

#' @export
setMethod("show", "parquetExprStore", function(object) {
    cat("An object of class parquetExprStore\n")
    cat(sprintf("  path     : %s\n", object@path))
    cat(sprintf("  uid      : %s\n", object@uid))
    cat(sprintf("  dim      : %s genes × %s cells\n",
        format(object@n_genes, big.mark = ",", scientific = FALSE),
        format(object@n_cells, big.mark = ",", scientific = FALSE)))
    if (length(object@feat_ids) > 0L) {
        n_show <- min(3L, length(object@feat_ids))
        cat(sprintf("  feat_ids : %s%s\n",
            paste(object@feat_ids[seq_len(n_show)], collapse = ", "),
            if (length(object@feat_ids) > 3L) ", ..." else ""))
    }
    if (length(object@cell_ids) > 0L) {
        n_show <- min(3L, length(object@cell_ids))
        cat(sprintf("  cell_ids : %s%s\n",
            paste(object@cell_ids[seq_len(n_show)], collapse = ", "),
            if (length(object@cell_ids) > 3L) ", ..." else ""))
    }
    cat(sprintf("  chunk    : %s cells\n",
        format(object@chunk_size, big.mark = ",", scientific = FALSE)))
    invisible(object)
})
