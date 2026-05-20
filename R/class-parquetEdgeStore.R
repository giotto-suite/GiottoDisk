#' @include class-dataStore.R class-parquetStore.R class-fileInputs.R
NULL

# DRAFT: parquetEdgeStore — disk-backed network store.
#
# Mirrors parquetExprStore in shape. Stores network edges in long-format
# parquet with a node-ID sidecar parquet for char<->int translation.
# Subsetting via @ops (no cell_idx slot). Materialization via
# storeRead(output = "igraph") performs match() rank in R on the
# subsetted edge set so the in-mem igraph only ever sees int32 IDs over
# the relevant subgraph.


# CLASS ####

#' @name parquetEdgeStore-class
#' @title Parquet Edge Store (streaming)
#' @description
#' S4 class for **disk-backed networks** stored as long-format Apache
#' Parquet. Carries a node-ID sidecar parquet for char<->int translation
#' so the in-memory class footprint stays bounded regardless of vertex
#' count.
#'
#' Edge schema (locked):
#'
#' | Column     | Type      | Meaning                          |
#' |------------|-----------|----------------------------------|
#' | `from_id`  | int32/64  | node int ID; canonical sort key  |
#' | `to_id`    | int32/64  | node int ID                      |
#' | `weight`   | float32   | edge weight (e.g. sNN Jaccard)   |
#' | `distance` | float32   | PCA / spatial Euclidean distance |
#'
#' On-disk layout:
#'
#' ```
#' <store_root>/                  # @path
#' ├── edges/                     # arrow::open_dataset reads here
#' │   └── *.parquet              # single file today; hive-partitioned later
#' └── nodes/                     # auto-derived sidecar
#'     └── *.parquet
#' ```
#'
#' Both subdirs are read via `arrow::open_dataset`, which handles
#' file-or-dir transparently — extending to `source_id=<uid>` hive
#' partitions for rbind support is a write-time change only.
#'
#' Node sidecar schema (lives at `<store_root>/nodes/`):
#'
#' | Column    | Type      | Required | Meaning                            |
#' |-----------|-----------|----------|------------------------------------|
#' | `node_id` | string    | yes      | original cell barcode / point ID   |
#' | `int_id`  | int32/64  | yes      | dense integer enumeration          |
#' | `*`       | any       | no       | arbitrary writer-supplied metadata |
#'
#' Optional sidecar columns survive round-trips and can be joined into
#' results at read time. Spatial coords (`x_index`, `y_index`) are NOT
#' duplicated here — they live in `parquetGeomStore` and are joined
#' across stores when needed.
#'
#' For undirected networks (sNN, spatial) edges are stored once in
#' canonical form (`from_id <= to_id`). kNN networks are symmetrized
#' at write time and stored canonical. Algorithms requiring directed
#' traversal can union with a swap.
#'
#' @slot n_cells numeric. Vertex count (= nrow nodes sidecar).
#' @slot n_edges numeric. Edge count.
#' @slot nodes parquetStore handle to the node sidecar.
#' @slot type character. One of "kNN", "sNN", "spatial".
#' @slot directed logical. TRUE = edges stored directed as-is;
#'   FALSE = canonical undirected form.
#' @family store types
#' @seealso [parquetEdgeStore()]
NULL

setClass("parquetEdgeStore",
    contains = "queryableStore",
    slots = list(
        n_cells  = "numeric",
        n_edges  = "numeric",
        nodes    = "ANY",
        type     = "character",
        directed = "logical"
    ),
    prototype = list(
        n_cells  = 0,
        n_edges  = 0,
        nodes    = NULL,
        type     = NA_character_,
        directed = FALSE
    )
)


# CONSTRUCTOR ####

#' @name parquetEdgeStore
#' @title Create a Parquet Edge Store handle
#' @description
#' Construct a [parquetEdgeStore-class] handle around an existing store
#' root, or a yet-to-be-written one (the directory and child parquets
#' are materialized later by [storeWrite()]). The node sidecar handle
#' is auto-derived from the store root — caller does not pass it.
#'
#' Typically not called directly — produced by [sourceWrite()] dispatching
#' on an edge input marker ([edgeDTInput()], [igraphInput()],
#' [nnSearchInput()]).
#' @param path character. Path to the store root directory. May not
#'   exist yet; [storeWrite()] creates it. Edges + nodes parquets are
#'   written into `<path>/edges/` and `<path>/nodes/` respectively.
#' @param type character. Network type ("kNN" / "sNN" / "spatial").
#' @param directed logical. Storage convention. Default FALSE.
#' @param n_edges numeric. Edge count. Auto-counted lazily if NA.
#' @param ... additional slots passed to `new()`.
#' @return [parquetEdgeStore-class] object.
#' @family store constructors
#' @export
parquetEdgeStore <- function(
    path     = .dump_tempfile(),
    type     = c("kNN", "sNN", "spatial"),
    directed = FALSE,
    n_edges  = NA_real_,
    ...
) {
    type <- match.arg(type)
    new("parquetEdgeStore",
        path     = path,
        n_edges  = as.numeric(n_edges),
        type     = type,
        directed = as.logical(directed),
        ...
    )
}


# INITIALIZE ####

setMethod("initialize", signature("parquetEdgeStore"),
    function(.Object, ...) {
        .Object <- callNextMethod(.Object, ...)

        # read_fun reads <path>/edges/ — handles file-or-dir
        if (.is_empty_fun(.Object@read_fun)) {
            .Object@read_fun <- function(x, ...) {
                arrow::open_dataset(sources = file.path(x, "edges"), ...)
            }
        }

        # auto-derive nodes sidecar handle at <path>/nodes/.
        # The directory may not exist yet — that's fine; the handle is
        # valid as soon as storeWrite materializes the file. parquetStore
        # initializes lazily, so a non-existent path doesn't error here.
        if (is.null(.Object@nodes) && length(.Object@path) > 0L) {
            .Object@nodes <- parquetStore(
                path = file.path(.Object@path, "nodes")
            )
        }

        # n_cells from sidecar if it materialized
        if (.Object@n_cells == 0 && !is.null(.Object@nodes)) {
            .Object@n_cells <- tryCatch(
                as.numeric(.dplyr_nrow(storeRead(.Object@nodes))),
                error = function(e) 0
            )
        }

        # lazy n_edges
        if (is.na(.Object@n_edges)) {
            .Object@n_edges <- tryCatch(
                as.numeric(.dplyr_nrow(.Object@read_fun(.Object@path))),
                error = function(e) NA_real_
            )
        }
        .Object
    }
)


# INPUT MARKERS ####
#
# Same convention as class-fileInputs.R: edgeInput is a VIRTUAL base;
# concrete markers wrap an in-mem source via @params$data. In-mem
# inputs are a nonstandard escape hatch in the framework but the right
# shape for network construction — networks are typically produced
# in-memory by RANN / dbscan / igraph and never serialize to a 10x-
# style standard layout.


#' @name edgeInput-class
#' @title Base class for edge-table inputs
#' @description Virtual parent of [edgeDTInput-class], [igraphInput-class],
#'   [nnSearchInput-class].
#' @slot from_col,to_col character. Endpoint columns in the source.
#' @slot weight_col,distance_col character or NULL. Optional attr cols.
#' @slot node_meta data.frame or NULL. Optional per-node metadata
#'   joined into the sidecar at write time. Must have a `node_id`
#'   column matching the input's vertex IDs.
#' @family store types
NULL

setClass("edgeInput",
    contains = c("fileStore", "VIRTUAL"),
    slots = list(
        from_col     = "character",
        to_col       = "character",
        weight_col   = "ANY",
        distance_col = "ANY",
        node_meta    = "ANY"
    ),
    prototype = list(
        from_col     = "from",
        to_col       = "to",
        weight_col   = NULL,
        distance_col = NULL,
        node_meta    = NULL
    )
)


# * edgeDTInput ####
setClass("edgeDTInput", contains = "edgeInput")

#' @name edgeDTInput
#' @title Wrap an in-memory data.table as an edge input
#' @param dt data.table with edge rows.
#' @param from_col,to_col character. Column names holding endpoint IDs.
#' @param weight_col,distance_col character or NULL.
#' @param node_meta data.frame or NULL. Optional per-node sidecar
#'   metadata; must have a `node_id` column.
#' @return [edgeDTInput-class] object.
#' @family store constructors
#' @export
edgeDTInput <- function(dt,
                        from_col     = "from",
                        to_col       = "to",
                        weight_col   = NULL,
                        distance_col = NULL,
                        node_meta    = NULL) {
    checkmate::assert_data_frame(dt, min.rows = 1L)
    stopifnot(from_col %in% names(dt), to_col %in% names(dt))
    if (!is.null(node_meta)) {
        checkmate::assert_data_frame(node_meta)
        stopifnot("node_id" %in% names(node_meta))
    }
    new("edgeDTInput",
        path         = .dump_tempfile(),
        params       = list(data = dt),
        from_col     = from_col,
        to_col       = to_col,
        weight_col   = weight_col,
        distance_col = distance_col,
        node_meta    = node_meta
    )
}


# * igraphInput ####
setClass("igraphInput", contains = "edgeInput")

#' @name igraphInput
#' @title Wrap an in-memory igraph as an edge input
#' @description
#' Extracts edges via `igraph::as_data_frame()` at write time; vertex
#' `name` attr (if present) becomes the `node_id` sidecar column; all
#' other vertex attrs become optional sidecar columns.
#' @param g `igraph` object.
#' @param weight_attr,distance_attr character or NULL.
#' @return [igraphInput-class] object.
#' @family store constructors
#' @export
igraphInput <- function(g,
                        weight_attr   = "weight",
                        distance_attr = NULL) {
    checkmate::assert_class(g, "igraph")
    # pull all vertex attrs as a node_meta data.frame
    v_attrs <- igraph::vertex_attr(g)
    if (length(v_attrs) > 0L && !"name" %in% names(v_attrs)) {
        v_attrs$name <- as.character(seq_len(igraph::vcount(g)))
    }
    node_meta <- if (length(v_attrs) > 0L) {
        df <- as.data.frame(v_attrs, stringsAsFactors = FALSE)
        names(df)[names(df) == "name"] <- "node_id"
        df
    } else NULL
    new("igraphInput",
        path         = .dump_tempfile(),
        params       = list(data = g),
        from_col     = "from",
        to_col       = "to",
        weight_col   = weight_attr,
        distance_col = distance_attr,
        node_meta    = node_meta
    )
}


# * nnSearchInput ####
setClass("nnSearchInput", contains = "edgeInput")

#' @name nnSearchInput
#' @title Wrap an RANN / dbscan kNN result as an edge input
#' @description
#' Expects a list with `nn.idx` (n x k integer) and `nn.dists`
#' (n x k float). Each input row expands to k edges at write time.
#' @param nn list with `nn.idx`, `nn.dists`.
#' @param cell_ids character. Vertex labels matching rows of `nn.idx`.
#' @param node_meta data.frame or NULL. Optional per-node sidecar.
#' @return [nnSearchInput-class] object.
#' @family store constructors
#' @export
nnSearchInput <- function(nn, cell_ids, node_meta = NULL) {
    stopifnot(is.list(nn), all(c("nn.idx", "nn.dists") %in% names(nn)))
    if (length(cell_ids) != nrow(nn$nn.idx)) {
        stop("[nnSearchInput] length(cell_ids) must equal nrow(nn$nn.idx)",
             call. = FALSE)
    }
    if (!is.null(node_meta)) {
        checkmate::assert_data_frame(node_meta)
        stopifnot("node_id" %in% names(node_meta))
    }
    new("nnSearchInput",
        path         = .dump_tempfile(),
        params       = list(data = nn, cell_ids = cell_ids),
        from_col     = "from",
        to_col       = "to",
        weight_col   = NULL,
        distance_col = "distance",
        node_meta    = node_meta
    )
}


# STOREWRITE DISPATCH ####
#
# sourceWrite(gsource, inp, store_type = "parquetEdge", ...) routes to
# storeWrite(parquetEdgeStore, inp, ...). Shape:
#   (i) coerce input -> data.table of edges with string IDs + node universe
#   (ii) build node sidecar with dense int enumeration (+ optional meta)
#   (iii) substitute int IDs into edges
#   (iv) canonicalize (sort + dedup if undirected)
#   (v) write edge parquet

#' @rdname storeWrite
setMethod("storeWrite",
    signature(store = "parquetEdgeStore", data = "edgeInput"),
    function(store, data,
             type = c("kNN", "sNN", "spatial"),
             directed = FALSE, ...) {
        type <- match.arg(type)

        # --- (i) coerce input ---------------------------------------
        ext <- .edge_input_to_dt(data)
        edges <- ext$edges       # data.table: <from_col>, <to_col>, [weight, distance]
        nodes <- ext$node_ids    # character vector — full vertex universe

        # --- (ii) build node sidecar -------------------------------
        nodes_dt <- data.table::data.table(
            node_id = as.character(nodes),
            int_id  = seq_along(nodes)
        )
        # join optional writer-supplied metadata onto sidecar
        if (!is.null(data@node_meta)) {
            meta <- data.table::as.data.table(data@node_meta)
            nodes_dt <- merge(nodes_dt, meta, by = "node_id",
                              all.x = TRUE, sort = FALSE)
        }

        # Materialize the directory layout under store@path. Edges +
        # nodes live in their own subdirs so future hive partitions
        # (source_id=<uid>/...) drop in without API change.
        edge_dir <- file.path(store@path, "edges")
        node_dir <- file.path(store@path, "nodes")
        dir.create(edge_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(node_dir, recursive = TRUE, showWarnings = FALSE)
        node_path <- file.path(node_dir, "nodes.parquet")
        edge_path <- file.path(edge_dir, "edges.parquet")
        arrow::write_parquet(nodes_dt, sink = node_path)

        # --- (iii) char -> int via match() (faster than named-vec) -
        edges[, from_id := match(as.character(get(data@from_col)),
                                  nodes_dt$node_id)]
        edges[, to_id   := match(as.character(get(data@to_col)),
                                  nodes_dt$node_id)]

        # --- (iv) canonicalize -------------------------------------
        if (!isTRUE(directed)) {
            swap <- edges$from_id > edges$to_id
            tmp  <- edges$from_id[swap]
            edges$from_id[swap] <- edges$to_id[swap]
            edges$to_id[swap]   <- tmp
            edges <- unique(edges, by = c("from_id", "to_id"))
        }
        data.table::setorder(edges, from_id, to_id)

        # --- (v) write edge parquet --------------------------------
        keep_cols <- c("from_id", "to_id")
        if (!is.null(data@weight_col) && data@weight_col %in% names(edges)) {
            data.table::setnames(edges, data@weight_col, "weight",
                                  skip_absent = TRUE)
            keep_cols <- c(keep_cols, "weight")
        }
        if (!is.null(data@distance_col) && data@distance_col %in% names(edges)) {
            data.table::setnames(edges, data@distance_col, "distance",
                                  skip_absent = TRUE)
            keep_cols <- c(keep_cols, "distance")
        }
        arrow::write_parquet(edges[, ..keep_cols], sink = edge_path)

        # --- assemble store handle ---------------------------------
        # @nodes was already created by initialize() pointing at
        # <path>/nodes/. The file just landed there — handle is valid.
        # Re-initialize the nodes handle so its lazy n_cells / metadata
        # picks up the freshly-written file.
        store@nodes    <- parquetStore(path = node_dir)
        store@n_cells  <- as.numeric(nrow(nodes_dt))
        store@n_edges  <- as.numeric(nrow(edges))
        store@type     <- type
        store@directed <- as.logical(directed)
        # @path stays the store root — read_fun knows to read
        # <path>/edges/. Caller can roundtrip via parquetEdgeStore(path).
        store
    }
)


# helper — pull edges out of any edgeInput subclass
.edge_input_to_dt <- function(inp) {
    if (inherits(inp, "edgeDTInput")) {
        dt <- data.table::as.data.table(inp@params$data)
        node_universe <- if (!is.null(inp@node_meta)) {
            # writer-supplied universe takes precedence (preserves nodes
            # with zero edges, e.g. isolated vertices)
            as.character(inp@node_meta$node_id)
        } else {
            as.character(unique(c(dt[[inp@from_col]], dt[[inp@to_col]])))
        }
        list(edges = dt, node_ids = node_universe)

    } else if (inherits(inp, "igraphInput")) {
        g <- inp@params$data
        dt <- data.table::as.data.table(
            igraph::as_data_frame(g, what = "edges")
        )
        node_universe <- igraph::V(g)$name %||%
            as.character(seq_len(igraph::vcount(g)))
        list(edges = dt, node_ids = node_universe)

    } else if (inherits(inp, "nnSearchInput")) {
        nn <- inp@params$data
        cell_ids <- inp@params$cell_ids
        n  <- nrow(nn$nn.idx); k <- ncol(nn$nn.idx)
        dt <- data.table::data.table(
            from     = rep(cell_ids, each = k),
            to       = cell_ids[as.integer(t(nn$nn.idx))],
            distance = as.numeric(t(nn$nn.dists))
        )
        list(edges = dt, node_ids = as.character(cell_ids))

    } else {
        stop("[.edge_input_to_dt] unhandled input class: ",
             class(inp)[[1L]], call. = FALSE)
    }
}


# STOREREAD — three output modes ####
#
# arrow    : lazy dataset, raw int IDs
# tibble   : collected data.table, character IDs (sidecar-joined)
# igraph   : in-mem igraph; int internals + character V(g)$name
#
# All three apply @ops first (subsetting). For tibble + igraph, the
# rank step uses match() on the subsetted edge set so the in-mem
# materialized object only sees int32 IDs over the relevant subgraph.

#' @rdname storeRead
setMethod("storeRead", signature(store = "parquetEdgeStore"),
    function(store,
             output  = c("arrow", "tibble", "igraph"),
             minimal = TRUE,
             ...) {
        output <- match.arg(output)

        # apply @ops via inherited queryableStore path; pass through
        # as "query" then take over for the materialization modes
        ds <- callNextMethod(store, output = "query", ...)

        if (output == "arrow") return(ds)

        # collect with attr filtering for igraph minimal mode
        if (output == "igraph" && isTRUE(minimal)) {
            keep <- c("from_id", "to_id")
            if ("weight" %in% names(ds)) keep <- c(keep, "weight")
            ds <- dplyr::select(ds, dplyr::all_of(keep))
        }
        edges_dt <- data.table::as.data.table(dplyr::collect(ds))

        if (output == "tibble") {
            return(.edge_relabel_char(edges_dt, store@nodes))
        }

        # igraph
        .edge_to_igraph(edges_dt, store@nodes, directed = store@directed)
    }
)


# Replace int_id columns with character node_id via sidecar join.
# Returns a fresh data.table preserving any extra columns (weight, distance).
.edge_relabel_char <- function(edges_dt, nodes) {
    int_id <- node_id <- NULL  # NSE
    used <- sort(unique(c(edges_dt$from_id, edges_dt$to_id)))
    name_map <- storeRead(nodes) |>
        dplyr::filter(int_id %in% used) |>
        dplyr::select(int_id, node_id) |>
        dplyr::collect() |>
        data.table::as.data.table()
    data.table::setkey(name_map, int_id)
    edges_dt$from_id <- name_map[J(edges_dt$from_id), node_id]
    edges_dt$to_id   <- name_map[J(edges_dt$to_id),   node_id]
    edges_dt
}


# Build an igraph from a subsetted edge table.
# Strategy: rank node ints down to a dense 1..V range over the actually
# referenced vertices (so the in-mem igraph never sees billion-scale
# ids). Set V(g)$name from the sidecar so callers can address vertices
# by their original character barcodes.
.edge_to_igraph <- function(edges_dt, nodes, directed = FALSE) {
    int_id <- node_id <- NULL  # NSE
    used_ints <- sort(unique(c(edges_dt$from_id, edges_dt$to_id)))

    # match() rank — faster than setNames(seq..., used)[as.character(...)]
    edges_dt$from_id <- match(edges_dt$from_id, used_ints)
    edges_dt$to_id   <- match(edges_dt$to_id,   used_ints)

    # vertex names via sidecar join — preserves cross-session stability
    names_dt <- storeRead(nodes) |>
        dplyr::filter(int_id %in% used_ints) |>
        dplyr::select(int_id, node_id) |>
        dplyr::collect() |>
        data.table::as.data.table()
    data.table::setkey(names_dt, int_id)
    v_df <- data.frame(name = names_dt[J(used_ints), node_id],
                       stringsAsFactors = FALSE)

    igraph::graph_from_data_frame(edges_dt,
                                  directed = directed,
                                  vertices = v_df)
}


# SUBSETTING ####
#
# Igraph-shaped ergonomics. All push down as @ops "filter":
#   x[v_set]                          # induced subgraph
#   x[from = v_set]                   # endpoint filter
#   x[from = v_set, to = u_set, negate = TRUE]

setMethod("[", signature(x = "parquetEdgeStore"),
    function(x, i, j, ...,
             from = NULL, to = NULL,
             negate = FALSE, drop = FALSE) {
        from_id <- to_id <- int_id <- node_id <- NULL  # NSE

        as_int <- function(ids) {
            if (is.null(ids)) return(NULL)
            if (is.numeric(ids)) return(as.integer(ids))
            # character -> sidecar lookup
            ns <- storeRead(x@nodes) |>
                dplyr::filter(node_id %in% ids) |>
                dplyr::select(int_id) |>
                dplyr::collect()
            as.integer(ns$int_id)
        }

        # induced subgraph: i = vertex set, both endpoints in set
        if (!missing(i) && is.null(from) && is.null(to)) {
            v_int <- as_int(i)
            expr <- if (isTRUE(negate)) {
                bquote(!(from_id %in% .(v_int) & to_id %in% .(v_int)))
            } else {
                bquote(from_id %in% .(v_int) & to_id %in% .(v_int))
            }
            x@ops <- c(x@ops, list(list(type = "filter", expr = expr)))
            return(x)
        }

        # endpoint filter
        f_int <- as_int(from)
        t_int <- as_int(to)
        clauses <- list()
        if (!is.null(f_int)) clauses <- c(clauses,
            list(bquote(from_id %in% .(f_int))))
        if (!is.null(t_int)) clauses <- c(clauses,
            list(bquote(to_id %in% .(t_int))))
        if (length(clauses) == 0L) return(x)
        expr <- Reduce(function(a, b) bquote(.(a) & .(b)), clauses)
        if (isTRUE(negate)) expr <- bquote(!(.(expr)))
        x@ops <- c(x@ops, list(list(type = "filter", expr = expr)))
        x
    }
)


# DIM / SHOW ####

setMethod("dim", signature(x = "parquetEdgeStore"), function(x) {
    c(as.integer(x@n_edges), 4L)
})

setMethod("nrow", signature(x = "parquetEdgeStore"), function(x) {
    as.integer(x@n_edges)
})

setMethod("show", signature(object = "parquetEdgeStore"), function(object) {
    cat(sprintf("<parquetEdgeStore> type=%s directed=%s\n",
                object@type, object@directed))
    cat(sprintf("  n_cells: %s  n_edges: %s\n",
                format(object@n_cells, big.mark = ","),
                format(object@n_edges, big.mark = ",")))
    cat(sprintf("  path:    %s\n", object@path))
    if (length(object@ops) > 0L) {
        cat(sprintf("  ops:     %d pending\n", length(object@ops)))
    }
})