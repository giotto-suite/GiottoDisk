#' @include class-fileInputs.R class-parquetEdgeStore.R methods-fileInputs.R
NULL

# AnnData .h5ad ####
# Pure-R reading of AnnData files, no python. Building blocks, not a reader:
# each returns an input for `sourceWrite()` or a Giotto subobject, and an
# object is assembled with ordinary `setGiotto()` calls.
#
# Contents:
#   accessor layer   .h5ad_*: the only code that touches HDF5 (hdf5r); a zarr
#                    backend can supply the same helpers later
#   h5adInput        exprInput over X, raw/X or layers/<name> (CSR or CSC)
#   h5adEdgeInput    edgeInput over obsp graphs, streamed to parquetEdgeStore
#   loaders          h5adCellMeta, h5adFeatMeta, h5adSpatLocs, h5adDimReds,
#                    h5adNetworks: one per slot
#
# Tuning: batches are sized in records (stored non-zeros) by
# `.ingest_batch_records()`. Per call with `batch_records`; per session with
# options(giottodisk.ingest_batch_records) or options(giottodisk.chunk_ram_frac).
#
# Orientation: AnnData is obs x var, row-major. hdf5r reports dims in R order,
# reversed from the on-disk shape, so a (n_obs, k) array reads as k x n_obs
# and is transposed back here. Sparse groups carry their (n_obs, n_var) shape
# as an attribute.


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

.h5ad_read_array <- function(node) {
    x <- node$read()
    if (length(node$dims) == 2L) x <- t(x)      # reversed dims, see header
    if (is.factor(x) && all(levels(x) %in% c("FALSE", "TRUE"))) {
        x <- as.logical(as.character(x))         # bool enum on older hdf5r
    }
    x
}

# Read one AnnData element by its encoding-type.
.h5ad_read_elem <- function(node, path = node$get_obj_name()) {
    enc <- .h5ad_encoding(node)
    switch(enc,
        "array" = , "string-array" = {
            if (node$attr_exists("categories")) {
                # anndata < 0.8: codes here, categories behind a reference
                warning("[h5ad] '", path, "' is a pre-0.8 categorical; ",
                        "returning its integer codes", call. = FALSE)
            }
            .h5ad_read_array(node)
        },
        "categorical" = {
            codes <- as.integer(node[["codes"]][])
            cats  <- as.character(node[["categories"]][])
            codes[codes < 0L] <- NA_integer_         # -1 is missing
            factor(cats[codes + 1L], levels = cats)
        },
        "nullable-integer" = , "nullable-boolean" = ,
        "nullable-string-array" = {
            v <- node[["values"]][]
            v[as.logical(node[["mask"]][])] <- NA
            v
        },
        stop("[h5ad] '", path, "' has unsupported encoding '", enc, "'",
             call. = FALSE)
    )
}

# obs / var -> data.table, columns in `column-order`, index first as `id_col`.
.h5ad_dataframe <- function(h5, group, id_col, columns = NULL) {
    g <- h5[[group]]
    cols <- as.character(.h5ad_attr(g, "column-order") %||% character(0L))
    if (!is.null(columns)) cols <- intersect(cols, columns)
    out <- data.table::data.table(.h5ad_index(h5, group))
    data.table::setnames(out, id_col)
    for (cn in setdiff(cols, id_col)) {
        data.table::set(out, j = cn,
            value = .h5ad_read_elem(g[[cn]], paste0(group, "/", cn)))
    }
    out[]
}

.h5ad_layers <- function(h5) {
    c("X",
      if (.h5ad_exists(h5, "raw/X")) "raw",
      if (.h5ad_exists(h5, "layers")) names(h5[["layers"]]))
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

# A call's `batch_records`, else the session pin or budget.
.h5ad_max_rows <- function(batch_records) {
    if (is.null(batch_records) || is.na(batch_records)) {
        .ingest_batch_records()
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
#' Wraps one sparse matrix of an `.h5ad` file: `X`, `raw/X` or
#' `layers/<name>`. Cell and feature IDs are read at construction; the matrix
#' is streamed by [storeRead()] in record-budgeted cell ranges. A cell-major
#' `csr_matrix` is read directly; a gene-major `csc_matrix` is reordered
#' through a temporary spill, as for Stereo-seq GEF files.
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
#'   `NULL` uses `options(giottodisk.ingest_batch_records)` when set, otherwise
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
    if (!info$encoding %in% c("csr_matrix", "csc_matrix")) {
        stop("[h5adInput] '/", info$root, "' is a dense '", info$encoding,
             "'; only sparse csr_matrix and csc_matrix are streamed.",
             call. = FALSE)
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

# CSR stores each cell's (gene, value) records contiguously, the same arrays a
# 10x h5 holds by cell, so contiguous cell ranges are single hyperslabs read
# by `.tenxh5_read_range()`. CSC stores them by gene: pass 1 streams gene
# ranges into cell buckets on disk, pass 2 returns the buckets in cell order
# (`.gef_spill_iterator()`), so the store is cell-major either way.

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("h5adInput"), function(store, ...) {
    row_id <- col_id <- NULL  # NSE bindings
    path     <- store@path
    root     <- store@params$root
    max_rows <- .h5ad_max_rows(store@batch_records)

    h5 <- .h5ad_open(path)
    indptr <- as.numeric(h5[[paste0(root, "/indptr")]][])
    h5$close_all()

    if (identical(store@params$encoding, "csr_matrix")) {
        ranges <- .gef_cell_ranges(diff(indptr), max_rows)
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
        # .gef_cell_ranges cuts any axis by record count; here, genes
        chunks <- .gef_cell_ranges(diff(indptr), max_rows)
        n_buckets <- .gef_n_buckets(indptr[length(indptr)], max_rows)
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
        itr <- .gef_spill_iterator(produce, length(chunks), n_buckets, finish, max_rows)
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
        ranges <- .gef_cell_ranges(diff(pat$indptr), .h5ad_max_rows(p$batch_records))

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


# loaders ####

# Each loader returns one Giotto subobject (or a named list of them) for
# `setGiotto()`. Reading goes through the accessor layer only.

#' @name h5adCellMeta
#' @title Cell metadata from an AnnData .h5ad
#' @description
#' Reads `obs` in its stored column order. Categoricals become factors (code
#' -1 as `NA`), nullable columns use their mask as `NA`, string arrays become
#' character and HDF5 boolean enums logical. The obs index becomes `cell_ID`.
#' @param h5ad_path character. Path to the `.h5ad` file.
#' @param columns character or NULL. Columns to read; `NULL` (default) reads
#'   all.
#' @param spat_unit,feat_type character. Spatial unit and feature type of the
#'   returned object.
#' @return A `cellMetaObj`.
#' @family h5ad readers
#' @export
h5adCellMeta <- function(h5ad_path, columns = NULL,
                         spat_unit = "cell", feat_type = "rna") {
    h5 <- .h5ad_open(h5ad_path)
    on.exit(h5$close_all())
    GiottoClass::createCellMetaObj(
        .h5ad_dataframe(h5, "obs", "cell_ID", columns),
        spat_unit = spat_unit, feat_type = feat_type,
        provenance = spat_unit, verbose = FALSE)
}

#' @name h5adFeatMeta
#' @title Feature metadata from an AnnData .h5ad
#' @description
#' Reads `var`, or `raw/var` for `layer = "raw"`, decoded as in
#' [h5adCellMeta()]. `feat_ID` matches [h5adInput()] for the same layer,
#' duplicate names included.
#' @inheritParams h5adCellMeta
#' @param layer character or NULL. As in [h5adInput()].
#' @return A `featMetaObj`.
#' @family h5ad readers
#' @export
h5adFeatMeta <- function(h5ad_path, layer = NULL, columns = NULL,
                         spat_unit = "cell", feat_type = "rna") {
    h5 <- .h5ad_open(h5ad_path)
    on.exit(h5$close_all())
    feat_ID <- NULL  # NSE binding
    info <- .h5ad_matrix_info(h5, layer)
    dt <- .h5ad_dataframe(h5, info$feat_group, "feat_ID", columns)
    dt[, feat_ID := .disambiguate_feat_ids(feat_ID)]
    GiottoClass::createFeatMetaObj(dt,
        spat_unit = spat_unit, feat_type = feat_type,
        provenance = spat_unit, verbose = FALSE)
}

#' @name h5adSpatLocs
#' @title Spatial locations from an AnnData .h5ad
#' @description
#' Reads a 2- or 3-column `obsm` array as `sdimx`, `sdimy` (and `sdimz`).
#' Coordinates are taken as stored.
#' @inheritParams h5adCellMeta
#' @param key character. Key under `obsm/` (default `"spatial"`).
#' @param name character. Name of the spatial locations (default `"raw"`).
#' @param flip_y logical. Negate y, for image-style (y-down) coordinates.
#' @param spat_unit character. Spatial unit of the returned object.
#' @return A `spatLocsObj`.
#' @family h5ad readers
#' @export
h5adSpatLocs <- function(h5ad_path, key = "spatial", name = "raw",
                         flip_y = FALSE, spat_unit = "cell") {
    sdimy <- cell_ID <- NULL  # NSE bindings
    h5 <- .h5ad_open(h5ad_path)
    on.exit(h5$close_all())
    p <- paste0("obsm/", key)
    if (!.h5ad_exists(h5, p)) {
        stop("[h5adSpatLocs] no coordinates at '/", p, "'.", call. = FALSE)
    }
    m <- .h5ad_read_elem(h5[[p]], p)
    if (!is.matrix(m) || !ncol(m) %in% 2:3) {
        stop("[h5adSpatLocs] '/", p, "' must be a 2- or 3-column array.",
             call. = FALSE)
    }
    dt <- data.table::as.data.table(m)
    data.table::setnames(dt, c("sdimx", "sdimy", "sdimz")[seq_len(ncol(m))])
    if (flip_y) dt[, sdimy := -sdimy]
    dt[, cell_ID := .h5ad_index(h5, "obs")]
    GiottoClass::createSpatLocsObj(dt, name = name, spat_unit = spat_unit,
        provenance = spat_unit, verbose = FALSE)
}

#' @name h5adDimReds
#' @title Dimension reductions from an AnnData .h5ad
#' @description
#' Reads `obsm` arrays other than spatial ones. `X_<name>` becomes a reduction
#' named `<name>`, with method `"pca"`, `"umap"` or `"tsne"` when the name
#' says so and the name otherwise. A PCA also takes its loadings from
#' `varm/PCs` (`varm/PCs_<name>` for another name) and its eigenvalues from
#' `uns/<name>/variance`, when present. Entries that are not numeric arrays
#' are skipped with a message.
#' @inheritParams h5adCellMeta
#' @param keys character or NULL. Keys under `obsm/`; `NULL` (default) reads
#'   every key without "spatial" in its name.
#' @param verbose logical. Report skipped keys.
#' @return A named list of `dimObj`.
#' @family h5ad readers
#' @export
h5adDimReds <- function(h5ad_path, keys = NULL,
                        spat_unit = "cell", feat_type = "rna",
                        verbose = TRUE) {
    h5 <- .h5ad_open(h5ad_path)
    on.exit(h5$close_all())
    if (!.h5ad_exists(h5, "obsm")) return(list())
    keys <- keys %||% grep("spatial", names(h5[["obsm"]]), value = TRUE,
                           invert = TRUE, ignore.case = TRUE)
    cell_ids <- .h5ad_index(h5, "obs")
    var_ids  <- NULL
    out <- list()
    for (k in keys) {
        p <- paste0("obsm/", k)
        if (!.h5ad_exists(h5, p) || !identical(.h5ad_encoding(h5[[p]]), "array")) {
            if (verbose) message("[h5adDimReds] skipping '/", p,
                                 "': not a numeric array")
            next
        }
        m <- .h5ad_read_elem(h5[[p]], p)
        if (!is.matrix(m)) next
        dimnames(m) <- list(cell_ids, paste0("Dim.", seq_len(ncol(m))))
        nm <- sub("^X_", "", k)
        method <- if (grepl("pca", nm, ignore.case = TRUE)) "pca"
            else if (grepl("umap", nm, ignore.case = TRUE)) "umap"
            else if (grepl("tsne", nm, ignore.case = TRUE)) "tsne"
            else nm
        misc <- NULL
        if (method == "pca") {
            lk <- paste0("varm/", if (nm == "pca") "PCs" else paste0("PCs_", nm))
            ek <- paste0("uns/", nm, "/variance")
            if (.h5ad_exists(h5, lk)) {
                var_ids <- var_ids %||% .disambiguate_feat_ids(.h5ad_index(h5, "var"))
                L <- .h5ad_read_elem(h5[[lk]], lk)
                dimnames(L) <- list(var_ids, colnames(m)[seq_len(ncol(L))])
                misc$loadings <- L
            }
            if (.h5ad_exists(h5, ek)) misc$eigenvalues <- as.numeric(h5[[ek]][])
        }
        out[[nm]] <- GiottoClass::createDimObj(m, name = nm,
            spat_unit = spat_unit, feat_type = feat_type, method = method,
            reduction = "cells", provenance = spat_unit, misc = misc)
    }
    out
}

#' @name h5adNetworks
#' @title Cell graphs from an AnnData .h5ad
#' @description
#' Finds the graphs described in `uns`: scanpy's `neighbors`, squidpy's
#' `spatial_neighbors`, or any group naming a `connectivities_key` and/or
#' `distances_key`. Each is streamed with [h5adEdgeInput()] into a
#' [parquetEdgeStore-class] of `gsource`.
#'
#' With distances, the graph is the directed kNN they describe, weighted by
#' the connectivity of each pair (0 where the connectivities omit it), with
#' self-loops dropped. With connectivities alone it is undirected. scanpy's
#' `neighbors` is named `kNN.<rep>` (`<rep>` from `params/use_rep`, default
#' `pca`); other keys keep their name. A key containing "spatial" gives a
#' spatial network, and one starting with `sNN` an sNN.
#' @inheritParams h5adCellMeta
#' @param gsource A [gDirSource-class] the edge stores are written to.
#' @param keys character or NULL. Keys under `uns/`; `NULL` (default) scans
#'   all.
#' @param batch_records numeric or NULL. As in [h5adEdgeInput()].
#' @param verbose logical. Report each graph and any that are skipped.
#' @return A named list of `nnNetObj` and `spatialNetworkObj`.
#' @family h5ad readers
#' @export
h5adNetworks <- function(h5ad_path, gsource, keys = NULL,
                         spat_unit = "cell", feat_type = "rna",
                         batch_records = NULL, verbose = TRUE) {
    h5 <- .h5ad_open(h5ad_path)
    graphs <- list()
    if (.h5ad_exists(h5, "uns")) {
        for (k in keys %||% names(h5[["uns"]])) {
            g <- paste0("uns/", k)
            if (!.h5ad_exists(h5, g) || !inherits(h5[[g]], "H5Group")) next
            ck <- if (.h5ad_exists(h5, paste0(g, "/connectivities_key"))) {
                as.character(h5[[paste0(g, "/connectivities_key")]][])
            }
            dk <- if (.h5ad_exists(h5, paste0(g, "/distances_key"))) {
                as.character(h5[[paste0(g, "/distances_key")]][])
            }
            if (is.null(ck) && is.null(dk)) next
            rep <- paste0(g, "/params/use_rep")
            rep <- if (.h5ad_exists(h5, rep)) {
                sub("^X_", "", as.character(h5[[rep]][]))
            } else "pca"
            type <- if (grepl("spatial", k, ignore.case = TRUE)) "spatial"
                else if (startsWith(k, "sNN")) "sNN" else "kNN"
            graphs[[length(graphs) + 1L]] <- list(
                key = k, ck = ck, dk = dk, type = type,
                name = if (k == "neighbors") paste0("kNN.", rep) else k)
        }
    }
    h5$close_all()

    out <- list()
    for (gr in graphs) {
        inp <- tryCatch(
            h5adEdgeInput(h5ad_path, connectivities_key = gr$ck,
                          distances_key = gr$dk, batch_records = batch_records),
            error = function(e) {
                if (verbose) message("[h5adNetworks] skipping uns/", gr$key,
                                     ": ", conditionMessage(e))
                NULL
            })
        if (is.null(inp)) next
        store <- sourceWrite(gsource, inp, store_type = "parquetEdge",
            type = gr$type, directed = inp@params$directed)
        out[[gr$name]] <- if (gr$type == "spatial") {
            GiottoClass::createSpatNetObj(network = store, name = gr$name,
                spat_unit = spat_unit, provenance = spat_unit)
        } else {
            GiottoClass::createNearestNetObj(name = gr$name, network = store,
                nn_type = gr$type, spat_unit = spat_unit, feat_type = feat_type)
        }
        if (verbose) message(sprintf("[h5adNetworks] uns/%s -> %s '%s' (%s edges)",
            gr$key, gr$type, gr$name, format(store@n_edges, big.mark = ",")))
    }
    out
}
