#' @include class-parquetExprStore.R
NULL

# storeRead ####
# Subset state lives in @cell_idx / @gene_idx (populated by `[`). Op chain
# lives on @ops as pure-data records. Both fold into the lazy Arrow query
# by wrapping @read_fun: subset filter first, then `.pe_apply_ops()`
# translates every op record into composed arrow steps (see .pe_do_op).
# The result is one composed arrow_dplyr_query executed once at output
# dispatch time. Delegated to queryableStore::storeRead for the
# query / tibble / duckdb output switch.

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetExprStore"), function(store, ...) {
    has_subset <- length(store@cell_idx) > 0L || length(store@gene_idx) > 0L
    has_ops    <- length(store@ops) > 0L
    if (has_subset || has_ops) {
        orig_rf <- store@read_fun
        ci  <- store@cell_idx
        gi  <- store@gene_idx
        ops <- store@ops
        store@read_fun <- function(x, ...) {
            row_id <- col_id <- NULL  # NSE bindings
            ds <- orig_rf(x, ...)
            if (length(ci) > 0L) ds <- dplyr::filter(ds, row_id %in% !!ci)
            if (length(gi) > 0L) ds <- dplyr::filter(ds, col_id %in% !!gi)
            if (length(ops) > 0L) ds <- .pe_apply_ops(ds, ops)
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
# storeWrite from a (possibly subset) parquetExprStore ####
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
    q <- .pestore_remap_query(data)
    .write_parquet(store, q)

    # Slots copied from input -- already correctly narrowed by any
    # `[`-subset. `cell_idx`/`gene_idx` stay empty on the fresh output
    # (no subset queued).
    store@cell_ids <- data@cell_ids
    store@feat_ids <- data@feat_ids
    store@n_cells  <- as.numeric(data@n_cells)
    store@n_genes  <- as.numeric(data@n_genes)
    .pestore_finalize_chunk_size(store)
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
            new_feat_ids <- x@feat_ids[i_int]
            x@feat_ids <- new_feat_ids
            x@gene_idx <- as.integer(new_gene_idx)
            x@n_genes  <- as.numeric(length(x@feat_ids))
            # Slice gene-axis op tables. Current op kinds have no
            # gene-axis tables; this is a no-op for norm_libsize_log but
            # generic for future kinds (e.g. center_scale_genes).
            if (length(x@ops) > 0L) {
                surviving_genes <- data.table::data.table(
                    feat_id = as.character(new_feat_ids))
                x@ops <- .pe_slice_ops(x@ops, axis = "gene",
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
            # Slice cell-axis op tables. Composite key (source_id,
            # orig_row_id); source_id is constant (= x@uid) for a single
            # store.
            if (length(x@ops) > 0L) {
                surviving_cells <- data.table::data.table(
                    source_id   = rep_len(as.character(x@uid),
                                          length(new_cell_idx)),
                    orig_row_id = as.integer(new_cell_idx))
                x@ops <- .pe_slice_ops(x@ops, axis = "cell",
                    surviving_keys = surviving_cells)
            }
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
        "query"  = atab,
        "tibble" = dplyr::collect(atab),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params)
    )
})

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
        parts <- list(bquote(source_id == .(s@uid)))
        if (length(s@cell_idx) > 0L) {
            parts[[length(parts) + 1L]] <- bquote(row_id %in% .(s@cell_idx))
        }
        if (length(s@gene_idx) > 0L) {
            parts[[length(parts) + 1L]] <- bquote(col_id %in% .(s@gene_idx))
        }
        clause <- Reduce(function(a, b) bquote(.(a) & .(b)), parts)
        clauses[[length(clauses) + 1L]] <- clause
    }
    if (length(no_subset_uids) > 0L) {
        clauses[[length(clauses) + 1L]] <- bquote(source_id %in% .(no_subset_uids))
    }
    Reduce(function(a, b) bquote(.(a) | .(b)), clauses)
}

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

        # Inherit + slice union-level @ops along the axes being narrowed.
        # Substore @ops stay empty by constraint; only the union carries
        # ops, and slicing applies on cell axis with the composite
        # (source_id, orig_row_id) key spanning surviving substores.
        new_union@ops <- x@ops
        if (length(new_union@ops) > 0L) {
            if (!missing(i)) {
                surviving_genes <- data.table::data.table(
                    feat_id = as.character(new_union@feat_ids))
                new_union@ops <- .pe_slice_ops(new_union@ops,
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
                new_union@ops <- .pe_slice_ops(new_union@ops,
                    axis = "cell", surviving_keys = surviving_cells)
            }
        }
        new_union
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
