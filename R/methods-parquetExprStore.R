#' @include class-parquetExprStore.R
NULL

# storeRead ####
# Subset state lives in @cell_idx / @gene_idx (populated by `[`). Injected
# into the lazy Arrow query by wrapping @read_fun, then delegated to
# queryableStore::storeRead which handles fields / callback / output
# dispatch (query / tibble / duckdb).

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetExprStore"), function(store, ...) {
    if (length(store@cell_idx) > 0L || length(store@gene_idx) > 0L) {
        orig_rf <- store@read_fun
        ci <- store@cell_idx
        gi <- store@gene_idx
        store@read_fun <- function(x, ...) {
            row_id <- col_id <- NULL  # NSE bindings
            ds <- orig_rf(x, ...)
            if (length(ci) > 0L) ds <- dplyr::filter(ds, row_id %in% !!ci)
            if (length(gi) > 0L) ds <- dplyr::filter(ds, col_id %in% !!gi)
            ds
        }
    }
    callNextMethod(store = store, ...)
})

# storeWrite ####
# from a dgCMatrix / Matrix / matrix.  Convenience path: useful for tests
# and small datasets that already live in memory.  The streaming Input
# classes (e.g. `mtxInput()`, `tenxH5Input()`) plus
# `storeWrite(parquetExprStore, exprInput)` are the production entry
# point for raw inputs — that path never materializes a dgCMatrix.

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
        # (rows = genes, cols = cells). scstream's row_id = cell index,
        # col_id = gene index. So row_id <- sm$j, col_id <- sm$i.
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
        store
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

# dimnames / rownames / colnames ####
# `rownames()` and `colnames()` in base R consult `dimnames()` first; defining
# `dimnames` here makes them work uniformly without separate methods.

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

# Translate a vector of subset positions to the original parquet
# row_ids / col_ids -- used when methods filter Arrow by HVG genes or
# cell bands and need the on-disk integer indices.

#' @keywords internal
#' @noRd
.pe_orig_row <- function(subset_pos, pe) {
    if (length(pe@cell_idx) == 0L) return(as.integer(subset_pos))
    as.integer(pe@cell_idx[subset_pos])
}

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
            x@feat_ids <- x@feat_ids[i_int]
            x@gene_idx <- as.integer(new_gene_idx)
            x@n_genes  <- as.numeric(length(x@feat_ids))
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
        }
        x
    }
)


# unionParquetExprStore methods ####

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
    output = c("query", "tibble", "duckdb"),
    callback = NULL,
    duckdb_params = list(),
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "duckdb"))
    wrapped <- lapply(store@stores, function(s) {
        if (length(s@cell_idx) > 0L || length(s@gene_idx) > 0L) {
            orig_rf <- s@read_fun
            ci <- s@cell_idx
            gi <- s@gene_idx
            s@read_fun <- function(x, ...) {
                row_id <- col_id <- NULL  # NSE bindings
                ds <- orig_rf(x, ...)
                if (length(ci) > 0L) ds <- dplyr::filter(ds, row_id %in% !!ci)
                if (length(gi) > 0L) ds <- dplyr::filter(ds, col_id %in% !!gi)
                ds
            }
        }
        s
    })
    atab <- arrow::open_dataset(lapply(wrapped, .store_simple_read))
    if (!is.null(fields)) atab <- dplyr::select(atab, dplyr::all_of(fields))
    if (!is.null(callback)) atab <- callback(atab)
    switch(output,
        "query"  = atab,
        "tibble" = dplyr::collect(atab),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params)
    )
})

# dim / dimnames / nrow / ncol — same conventions as parquetExprStore.

#' @export
setMethod("nrow", "unionParquetExprStore", function(x) x@n_genes)

#' @export
setMethod("ncol", "unionParquetExprStore", function(x) x@n_cells)

#' @export
setMethod("dim", "unionParquetExprStore",
    function(x) c(x@n_genes, x@n_cells)
)

#' @export
setMethod("dimnames", "unionParquetExprStore",
    function(x) list(x@feat_ids, x@cell_ids)
)


# [ subset
# i (genes) — applied uniformly to all substores (feat_ids are shared).
# j (cells) — mapped from union positions to per-substore positions via
# cumulative offsets; substores that get zero cells after the subset
# are dropped. Result is rebuilt through the constructor for invariant
# checks.

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
        unionParquetExprStore(new_stores)
    }
)


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
