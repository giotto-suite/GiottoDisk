#' @include class-fileInputs.R class-parquetEdgeStore.R methods-fileInputs.R
NULL

# AnnData .h5ad ####
# Streaming inputs for `sourceWrite()`, no python. Decoding AnnData elements
# (obs/var frames, obsm, uns) and assembling an object is GiottoConnect's
# AnnData adapter; these inputs are what it calls for GiottoDisk targets.
#
# Contents:
#   accessor layer   .h5ad_*: the only code that touches HDF5 (hdf5r); a zarr
#                    backend can supply the same helpers later
#   h5adInput        exprInput over X, raw/X or layers/<name> (CSR, CSC, dense)
#   h5adEdgeInput    edgeInput over obsp graphs, streamed to parquetEdgeStore
#
# Tuning: batches are sized in records (stored non-zeros), as for a 10x .h5:
# per call with `batch_records`; per session with
# options(giottodisk.h5_batch_rows) or options(giottodisk.chunk_ram_frac).
#
# Orientation: AnnData is obs x var, row-major. hdf5r reports dims in R order,
# reversed from the on-disk shape. Sparse groups carry their (n_obs, n_var)
# shape as an attribute.


# accessor layer ####

.h5ad_open <- function(path) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("[h5ad] hdf5r is required to read .h5ad files. ",
             "Install with: install.packages(\"hdf5r\").", call. = FALSE)
    }
    hdf5r::H5File$new(path, mode = "r")
}

# Walk the path one level at a time: H5Lexists errors on a missing
# intermediate group rather than returning FALSE.
.h5ad_exists <- function(h5, p) {
    parts <- strsplit(p, "/", fixed = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    cur <- character(0L)
    for (s in parts) {
        cur <- paste(c(cur, s), collapse = "/")
        if (!h5$exists(cur)) return(FALSE)
    }
    TRUE
}

.h5ad_attr <- function(node, name) {
    if (!node$attr_exists(name)) return(NULL)
    hdf5r::h5attr(node, name)
}

.h5ad_encoding <- function(node) {
    enc <- .h5ad_attr(node, "encoding-type")
    if (!is.null(enc)) return(enc)
    if (inherits(node, "H5Group")) "dict" else "array"
}

# obs / var index: the group's `_index` attribute names the index dataset.
.h5ad_index <- function(h5, group) {
    g <- h5[[group]]
    as.character(g[[.h5ad_attr(g, "_index") %||% "_index"]][])
}

.h5ad_layers <- function(h5) {
    c("X",
      if (.h5ad_exists(h5, "raw/X")) "raw",
      if (.h5ad_exists(h5, "layers")) names(h5[["layers"]]))
}

# Matrices with their encoding, e.g. "X (array), counts (csr_matrix)".
.h5ad_layers_desc <- function(h5) {
    toString(vapply(.h5ad_layers(h5), function(l) {
        sprintf("%s (%s)", l, .h5ad_encoding(h5[[.h5ad_matrix_info(h5, l)$root]]))
    }, character(1L)))
}

.h5ad_sparse_layers <- function(h5) {
    Filter(function(l) {
        .h5ad_encoding(h5[[.h5ad_matrix_info(h5, l)$root]]) %in%
            c("csr_matrix", "csc_matrix")
    }, setdiff(.h5ad_layers(h5), "X"))
}

# Resolve a matrix: NULL/"X" -> /X over var, "raw" -> /raw/X over raw/var,
# anything else -> /layers/<layer> over var.
.h5ad_matrix_info <- function(h5, layer = NULL) {
    if (is.null(layer) || identical(layer, "X")) {
        root <- "X"; feat_group <- "var"
    } else if (identical(layer, "raw")) {
        root <- "raw/X"; feat_group <- "raw/var"
    } else {
        root <- paste0("layers/", layer); feat_group <- "var"
    }
    if (!.h5ad_exists(h5, root)) {
        stop("[h5ad] no matrix at '/", root, "'. Available: ",
             toString(.h5ad_layers(h5)), call. = FALSE)
    }
    node <- h5[[root]]
    shape <- if (inherits(node, "H5Group")) {
        as.numeric(.h5ad_attr(node, "shape"))
    } else {
        rev(as.numeric(node$dims))
    }
    list(root = root, encoding = .h5ad_encoding(node), shape = shape,
         feat_group = feat_group)
}

# A call's `batch_records`, else the 10x .h5 pin or budget: CSR ranges are
# read with the same `.tenxh5_read_range()`.
.h5ad_max_rows <- function(batch_records) {
    if (is.null(batch_records) || is.na(batch_records)) {
        .ingest_batch_rows(128, pinned = getOption("giottodisk.h5_batch_rows"))
    } else {
        as.numeric(batch_records)
    }
}

# Rows a..b of a CSR group as (from_id, to_id, value), 1-based. `ds` holds the
# open data/indices datasets and the full indptr.
.h5ad_csr_rows <- function(ds, a, b) {
    lo <- ds$indptr[a] + 1
    hi <- ds$indptr[b + 1L]
    if (hi < lo) return(NULL)
    # `lo:hi`, never `seq.int(lo, hi)` -- adr/0007.
    data.table::data.table(
        from_id = rep.int(seq.int(a, b), diff(ds$indptr[a:(b + 1L)])),
        to_id   = as.integer(ds$indices[lo:hi]) + 1L,
        value   = as.double(ds$data[lo:hi])
    )
}

.h5ad_csr_handle <- function(h5, root) {
    list(indptr  = as.numeric(h5[[paste0(root, "/indptr")]][]),
         data    = h5[[paste0(root, "/data")]],
         indices = h5[[paste0(root, "/indices")]])
}


# h5adInput ####

#' @name h5adInput-class
#' @title AnnData .h5ad expression input
#' @description
#' Wraps one matrix of an `.h5ad` file: `X`, `raw/X` or `layers/<name>`.
#' Cell and feature IDs are read at construction; the matrix is streamed by
#' [storeRead()] in record-budgeted cell ranges. A cell-major `csr_matrix` is
#' read directly; a gene-major `csc_matrix` is reordered through a temporary
#' spill, as for Stereo-seq GEF files; a dense array is read in blocks of
#' cells and only its non-zero values are stored.
#' @slot batch_records numeric. Records per batch; `NA` (default) sizes it
#'   from free RAM.
#' @family store types
NULL

setClass("h5adInput",
    contains = "exprInput",
    slots = list(batch_records = "numeric"),
    prototype = list(batch_records = NA_real_)
)

#' @name h5adInput
#' @title Create an AnnData .h5ad expression input
#' @param h5ad_path character. Path to the `.h5ad` file.
#' @param layer character or NULL. `NULL` or `"X"` (default) for `X`,
#'   `"raw"` for `raw/X`, otherwise a name under `layers/`.
#' @param batch_records numeric or NULL. Records (stored values) per batch.
#'   `NULL` uses `options(giottodisk.h5_batch_rows)` when set, otherwise
#'   a budget from free RAM (`giottodisk.chunk_ram_frac`); see
#'   [storeChunkInfo()].
#' @return An `h5adInput` object.
#' @family store constructors
#' @export
h5adInput <- function(h5ad_path, layer = NULL, batch_records = NULL) {
    stopifnot(file.exists(h5ad_path))
    h5 <- .h5ad_open(h5ad_path)
    on.exit(h5$close_all())

    info <- .h5ad_matrix_info(h5, layer)
    if (!info$encoding %in% c("csr_matrix", "csc_matrix", "array")) {
        stop("[h5adInput] '/", info$root, "' is a '", info$encoding,
             "'; expected csr_matrix, csc_matrix or a dense array. ",
             "Matrices in this file: ", .h5ad_layers_desc(h5), call. = FALSE)
    }
    if (identical(info$encoding, "array") && is.null(layer)) {
        sparse <- .h5ad_sparse_layers(h5)
        if (length(sparse)) {
            message("[h5adInput] X is dense; sparse layer(s) ",
                    toString(sQuote(sparse, FALSE)), " are available ",
                    "(`layer = `). Reading X in row blocks.")
        }
    }
    cell_ids <- .h5ad_index(h5, "obs")
    feat_ids <- .disambiguate_feat_ids(.h5ad_index(h5, info$feat_group))
    if (length(cell_ids) != info$shape[1L] ||
        length(feat_ids) != info$shape[2L]) {
        stop("[h5adInput] '/", info$root, "' has shape ",
             paste(info$shape, collapse = " x "), " but obs/",
             info$feat_group, " hold ", length(cell_ids), " x ",
             length(feat_ids), " ids.", call. = FALSE)
    }

    new("h5adInput",
        path          = normalizePath(h5ad_path),
        params        = list(root = info$root, encoding = info$encoding,
                             layer = layer %||% "X"),
        cell_ids      = cell_ids,
        feat_ids      = feat_ids,
        n_cells       = as.integer(info$shape[1L]),
        n_genes       = as.integer(info$shape[2L]),
        batch_records = as.numeric(batch_records %||% NA_real_)
    )
}

# storeRead -- h5adInput ####

# A dense matrix is stored obs x var, row-major, so a block of cells is one
# contiguous hyperslab. hdf5r reports it var x obs (see header): block
# [, a:b] is genes x cells, and which() walks it cell by cell, gene by gene,
# so the records come out in store order. Blocks hold `max_rows` values,
# zeros included, since that is what is read.
.h5ad_dense_reader <- function(store, max_rows) {
    path <- store@path
    root <- store@params$root
    n_cells <- as.integer(store@n_cells)
    per_block <- max(1L, as.integer(floor(max_rows / max(1L, store@n_genes))))
    starts <- seq.int(1L, n_cells, by = per_block)
    i <- 0L
    closed <- FALSE
    close_fn <- function() { closed <<- TRUE; invisible(NULL) }
    next_batch <- function() {
        repeat {
            if (closed || i >= length(starts)) { close_fn(); return(NULL) }
            i <<- i + 1L
            a <- starts[[i]]
            b <- min(n_cells, a + per_block - 1L)
            h5 <- .h5ad_open(path)
            m <- h5[[root]][, a:b, drop = FALSE]   # `a:b`: one hyperslab
            h5$close_all()
            nz <- which(m != 0, arr.ind = TRUE)
            if (!nrow(nz)) next
            return(data.table::data.table(
                row_id = as.integer(nz[, 2L] + a - 1L),
                col_id = as.integer(nz[, 1L]),
                value  = as.double(m[nz])
            ))
        }
    }
    list(
        next_batch = next_batch,
        close      = close_fn,
        cell_ids   = function() store@cell_ids,
        feat_ids   = function() store@feat_ids,
        n_cells    = function() store@n_cells,
        n_genes    = function() store@n_genes
    )
}

# CSR stores each cell's (gene, value) records contiguously, the same arrays a
# 10x h5 holds by cell, so contiguous cell ranges are single hyperslabs read
# by `.tenxh5_read_range()`. CSC stores them by gene: pass 1 streams gene
# ranges into cell buckets on disk, pass 2 returns the buckets in cell order
# (`.ingest_spill_iterator()`), so the store is cell-major either way.

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("h5adInput"), function(store, ...) {
    row_id <- col_id <- NULL  # NSE bindings
    path     <- store@path
    root     <- store@params$root
    max_rows <- .h5ad_max_rows(store@batch_records)

    if (identical(store@params$encoding, "array")) {
        return(.h5ad_dense_reader(store, max_rows))
    }

    h5 <- .h5ad_open(path)
    indptr <- as.numeric(h5[[paste0(root, "/indptr")]][])
    h5$close_all()

    if (identical(store@params$encoding, "csr_matrix")) {
        ranges <- .ingest_cell_ranges(diff(indptr), max_rows)
        i <- 0L
        closed <- FALSE
        close_fn <- function() { closed <<- TRUE; invisible(NULL) }
        next_batch <- function() {
            repeat {
                if (closed || i >= length(ranges)) { close_fn(); return(NULL) }
                i <<- i + 1L
                r <- ranges[[i]]
                dt <- .tenxh5_read_range(path, root, r[1L], r[2L], indptr)
                if (!is.null(dt) && nrow(dt)) return(dt)
            }
        }
    } else {
        n_cells <- as.integer(store@n_cells)
        # .ingest_cell_ranges cuts any axis by record count; here, genes
        chunks <- .ingest_cell_ranges(diff(indptr), max_rows)
        n_buckets <- .ingest_n_buckets(indptr[length(indptr)], max_rows)
        per_bucket <- max(1L, as.integer(ceiling(n_cells / n_buckets)))
        n_buckets <- as.integer(ceiling(n_cells / per_bucket))
        produce <- function(k) {
            g_lo <- chunks[[k]][1L]; g_hi <- chunks[[k]][2L]
            h5 <- .h5ad_open(path)
            on.exit(h5$close_all())
            dt <- .h5ad_csr_rows(.h5ad_csr_handle(h5, root), g_lo, g_hi)
            if (is.null(dt)) return(NULL)
            # rows of a CSC are genes, its indices are cells
            data.table::data.table(
                b      = (dt$to_id - 1L) %/% per_bucket + 1L,
                row_id = dt$to_id,
                col_id = dt$from_id,
                value  = dt$value
            )
        }
        finish <- function(dt) data.table::setorder(dt, row_id, col_id)[]
        itr <- .ingest_spill_iterator(produce, length(chunks), n_buckets, finish, max_rows)
        next_batch <- itr$next_batch
        close_fn <- itr$close
    }

    list(
        next_batch = next_batch,
        close      = close_fn,
        cell_ids   = function() store@cell_ids,
        feat_ids   = function() store@feat_ids,
        n_cells    = function() store@n_cells,
        n_genes    = function() store@n_genes
    )
})


# h5adEdgeInput ####

#' @name h5adEdgeInput-class
#' @title AnnData .h5ad graph input
#' @description
#' Wraps cell x cell graphs under `obsp/` of an `.h5ad` file for streaming
#' into a [parquetEdgeStore-class]. With a distances matrix, its sparsity
#' pattern is the edge set (a directed kNN) and the connectivity at the same
#' cell pair is the weight. With connectivities alone, the upper triangle is
#' the edge set of an undirected graph.
#' @family store types
NULL

setClass("h5adEdgeInput", contains = "edgeInput")

#' @name h5adEdgeInput
#' @title Create an AnnData .h5ad graph input
#' @param h5ad_path character. Path to the `.h5ad` file.
#' @param connectivities_key,distances_key character or NULL. Keys under
#'   `obsp/`. At least one is required.
#' @param batch_records numeric or NULL. Edges per batch, sized as in
#'   [h5adInput()].
#' @return An `h5adEdgeInput` object.
#' @family store constructors
#' @export
h5adEdgeInput <- function(h5ad_path,
                          connectivities_key = NULL,
                          distances_key      = NULL,
                          batch_records      = NULL) {
    stopifnot(file.exists(h5ad_path))
    if (is.null(connectivities_key) && is.null(distances_key)) {
        stop("[h5adEdgeInput] give connectivities_key, distances_key or both.",
             call. = FALSE)
    }
    h5 <- .h5ad_open(h5ad_path)
    on.exit(h5$close_all())
    n_obs <- length(.h5ad_index(h5, "obs"))
    for (k in c(connectivities_key, distances_key)) {
        p <- paste0("obsp/", k)
        if (!.h5ad_exists(h5, p)) {
            stop("[h5adEdgeInput] no graph at '/", p, "'.", call. = FALSE)
        }
        enc <- .h5ad_encoding(h5[[p]])
        if (!identical(enc, "csr_matrix")) {
            stop("[h5adEdgeInput] '/", p, "' is '", enc,
                 "'; only csr_matrix graphs are streamed.", call. = FALSE)
        }
        if (!all(as.numeric(.h5ad_attr(h5[[p]], "shape")) == n_obs)) {
            stop("[h5adEdgeInput] '/", p, "' is not ", n_obs, " x ", n_obs,
                 ".", call. = FALSE)
        }
    }
    directed <- !is.null(distances_key)
    new("h5adEdgeInput",
        path         = normalizePath(h5ad_path),
        params       = list(
            pattern  = if (directed) distances_key else connectivities_key,
            other    = if (directed) connectivities_key else NULL,
            directed = directed,
            batch_records = batch_records),
        weight_col   = if (!is.null(connectivities_key)) "weight",
        distance_col = if (directed) "distance"
    )
}

# storeWrite -- parquetEdgeStore <- h5adEdgeInput ####

# Same layout as `.edge_storewrite_core()` -- nodes/nodes.parquet plus
# edges/ -- but streamed. Nodes are the obs index in file order, so a CSR row
# and its column indices are already the dense 1-based ids and no match() is
# needed. Ranges ascend, so the edge parts are globally sorted by
# (from_id, to_id). Self-loops are dropped: some writers keep each cell as its
# own first neighbour at distance 0. A kNN pair missing from connectivities
# (UMAP's fuzzy graph can prune it) gets weight 0, its sparse value.

#' @rdname storeWrite
setMethod("storeWrite",
    signature(store = "parquetEdgeStore", data = "h5adEdgeInput"),
    function(store, data, type = c("kNN", "sNN", "spatial"),
             directed = data@params$directed, ...) {
        type <- match.arg(type)
        from_id <- to_id <- weight <- NULL  # NSE bindings
        p <- data@params

        h5 <- .h5ad_open(data@path)
        on.exit(h5$close_all())
        nodes <- .h5ad_index(h5, "obs")
        n <- length(nodes)
        nodes_dt <- data.table::data.table(
            row_index = seq_len(n),
            node_id   = nodes,
            int_id    = seq_len(n)
        )
        edge_dir <- file.path(store@path, "edges")
        node_dir <- file.path(store@path, "nodes")
        dir.create(edge_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(node_dir, recursive = TRUE, showWarnings = FALSE)
        .write_parquet_file(nodes_dt, file.path(node_dir, "nodes.parquet"))
        rm(nodes, nodes_dt)

        pat_col <- if (p$directed) "distance" else "weight"
        pat <- .h5ad_csr_handle(h5, paste0("obsp/", p$pattern))
        oth <- if (!is.null(p$other)) {
            .h5ad_csr_handle(h5, paste0("obsp/", p$other))
        }
        ranges <- .ingest_cell_ranges(diff(pat$indptr), .h5ad_max_rows(p$batch_records))

        n_edges <- 0
        n_parts <- 0L
        for (r in ranges) {
            e <- .h5ad_csr_rows(pat, r[1L], r[2L])
            if (is.null(e)) next
            e <- if (p$directed) e[from_id != to_id] else e[from_id < to_id]
            data.table::setnames(e, "value", pat_col)
            if (!is.null(oth)) {
                w <- .h5ad_csr_rows(oth, r[1L], r[2L])
                if (!is.null(w)) {
                    data.table::setnames(w, "value", "weight")
                    e <- w[e, on = c("from_id", "to_id")]
                    # absent from a sparse matrix is a stored zero
                    e[is.na(weight), weight := 0]
                } else {
                    e[, weight := 0]
                }
                data.table::setcolorder(e, c("from_id", "to_id", "weight", "distance"))
            }
            data.table::setorder(e, from_id, to_id)
            .write_parquet_file(e, file.path(
                edge_dir, sprintf("part-%d.parquet", n_parts)))
            n_parts <- n_parts + 1L
            n_edges <- n_edges + nrow(e)
        }

        store@nodes    <- parquetStore(path = node_dir)
        store@n_cells  <- as.numeric(n)
        store@n_edges  <- as.numeric(n_edges)
        store@type     <- type
        store@directed <- as.logical(directed)
        store
    }
)
