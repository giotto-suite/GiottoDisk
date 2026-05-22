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
    contains = c("queryableStore", "parquetBase"),
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
             directed = FALSE,
             node_universe = NULL,
             ...) {
        type <- match.arg(type)

        # Coerce input -> (edges DT, node universe). Caller is expected
        # to have produced canonical edges (e.g. .undirected_unique on
        # the Giotto side for sNN/Delaunay) before reaching here.
        ext <- .edge_input_to_dt(data)
        edges <- ext$edges
        nodes <- node_universe %null% ext$node_ids

        # Rename input's endpoint cols to canonical names so the core
        # writer can be column-agnostic across input markers.
        if (data@from_col != "from") {
            data.table::setnames(edges, data@from_col, "from")
        }
        if (data@to_col != "to") {
            data.table::setnames(edges, data@to_col, "to")
        }
        if (!is.null(data@weight_col) && data@weight_col %in% names(edges) &&
            data@weight_col != "weight") {
            data.table::setnames(edges, data@weight_col, "weight")
        }
        if (!is.null(data@distance_col) && data@distance_col %in% names(edges) &&
            data@distance_col != "distance") {
            data.table::setnames(edges, data@distance_col, "distance")
        }

        .edge_storewrite_core(
            store = store, edges = edges, nodes = nodes,
            node_meta = data@node_meta,
            type = type, directed = directed
        )
    }
)


#' @rdname storeWrite
#' @description
#' Direct data.table dispatch: the in-memory caller (e.g. Giotto's
#' `.finalize_network`) hands a pre-canonicalized edge data.table with
#' character `from` / `to` columns plus optional `weight` / `distance` /
#' `shared`. Trusts the caller for canonical form — no swap, no dedup.
setMethod("storeWrite",
    signature(store = "parquetEdgeStore", data = "data.table"),
    function(store, data,
             type = c("kNN", "sNN", "spatial"),
             directed = FALSE,
             node_universe = NULL,
             node_meta = NULL,
             ...) {
        type <- match.arg(type)
        if (!all(c("from", "to") %in% names(data))) {
            stop("[storeWrite(parquetEdgeStore, data.table)] `data` must ",
                 "have `from` and `to` columns.", call. = FALSE)
        }
        nodes <- node_universe %null%
            as.character(unique(c(data$from, data$to)))
        .edge_storewrite_core(
            store = store, edges = data, nodes = nodes,
            node_meta = node_meta,
            type = type, directed = directed
        )
    }
)


# Core write logic shared by edgeInput + data.table storeWrite paths.
# Assumes:
#   - `edges` has columns "from"/"to" with character node IDs, plus any
#     optional payload cols (weight, distance, shared, ...).
#   - `nodes` is the full vertex universe (character).
#   - Caller is responsible for canonicalization. `directed` here is
#     pure metadata — recorded on @directed for readers, never gates
#     processing.
.edge_storewrite_core <- function(store, edges, nodes, node_meta,
        type, directed) {
    from <- to <- NULL  # NSE

    # --- node sidecar build (int auto-promotion at 2^31 boundary) ---
    n <- length(nodes)
    int_ids <- if (n > .Machine$integer.max) {
        # seq_len(n) returns double for long vectors; Arrow writes int64
        # when the schema declares it. No bit64 dependency.
        arrow::Array$create(seq_len(n), type = arrow::int64())
    } else {
        seq_len(n)
    }
    nodes_dt <- data.table::data.table(
        row_index = seq_len(n),  # parquetStore contract; hidden via specialCols
        node_id   = as.character(nodes),
        int_id    = int_ids
    )
    if (!is.null(node_meta)) {
        meta <- data.table::as.data.table(node_meta)
        nodes_dt <- merge(nodes_dt, meta, by = "node_id",
                          all.x = TRUE, sort = FALSE)
        data.table::setorder(nodes_dt, row_index)
    }

    # --- directory layout ------------------------------------------
    edge_dir <- file.path(store@path, "edges")
    node_dir <- file.path(store@path, "nodes")
    dir.create(edge_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(node_dir, recursive = TRUE, showWarnings = FALSE)
    node_path <- file.path(node_dir, "nodes.parquet")
    edge_path <- file.path(edge_dir, "edges.parquet")
    arrow::write_parquet(nodes_dt, sink = node_path)

    # --- substitute char -> int via match() ------------------------
    edges <- data.table::copy(edges)
    edges[, from_id := match(as.character(from), nodes_dt$node_id)]
    edges[, to_id   := match(as.character(to),   nodes_dt$node_id)]
    edges[, c("from", "to") := NULL]
    data.table::setcolorder(edges, c("from_id", "to_id"))
    # Sort once on storage-canonical key.
    data.table::setorder(edges, from_id, to_id)

    # --- write edge parquet ----------------------------------------
    arrow::write_parquet(edges, sink = edge_path)

    # --- assemble store handle -------------------------------------
    store@nodes    <- parquetStore(path = node_dir)
    store@n_cells  <- as.numeric(nrow(nodes_dt))
    store@n_edges  <- as.numeric(nrow(edges))
    store@type     <- type
    store@directed <- as.logical(directed)
    store
}


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

        # Open the edges dataset via inherited queryableStore path,
        # then apply @ops manually (we don't extend parquetStore, so
        # .pbase_storeread_processing isn't in our dispatch chain).
        ds <- callNextMethod(store, output = "query", ...)
        for (op in store@ops) {
            ds <- .do_op(ds, op)
        }

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
    # vertices DF: first col is the vertex identifier (the ranks
    # 1..V that we just substituted into edges_dt). `name` is an
    # attribute igraph sets on V(g)$name automatically.
    v_df <- data.frame(
        id   = seq_along(used_ints),
        name = names_dt[J(used_ints), node_id],
        stringsAsFactors = FALSE
    )

    igraph::graph_from_data_frame(edges_dt,
                                  directed = directed,
                                  vertices = v_df)
}


# SUBSETTING ####
#
# `[` follows igraph's adjacency-style semantics, with one-arg shorthand
# for the common induced-subgraph case:
#
#   x[v_set]                # induced subgraph (both endpoints in v_set)
#   x[v_a, v_b]             # from/to slice (matrix-shape)
#   x[, v]                  # to-endpoint filter
#   x[from = v, to = u]     # explicit named form (same as x[v, u])
#   x[..., negate = TRUE]   # complement of the above
#
# For undirected stores (`@directed = FALSE`) we store edges in canonical
# `from_id <= to_id` form. All asymmetric slices OR both orientations
# so the user-visible behavior matches an in-mem igraph on a symmetric
# graph — the canonical-storage detail doesn't leak.

# helper — translate a character or numeric vertex-id vector to the
# int_id space used in the edges parquet, via sidecar lookup.
.edge_vset_to_int <- function(x, ids) {
    int_id <- node_id <- NULL  # NSE
    if (is.null(ids)) return(NULL)
    if (is.numeric(ids)) return(as.integer(ids))
    ns <- storeRead(x@nodes) |>
        dplyr::filter(node_id %in% ids) |>
        dplyr::select(int_id) |>
        dplyr::collect()
    as.integer(ns$int_id)
}

# helper — induced-subgraph filter (both endpoints in v_int).
# Works for canonical undirected storage without OR because the AND
# requirement catches both orientations regardless.
.edge_induced_op <- function(v_int, negate = FALSE) {
    from_id <- to_id <- NULL  # NSE
    expr <- bquote(from_id %in% .(v_int) & to_id %in% .(v_int))
    if (isTRUE(negate)) expr <- bquote(!(.(expr)))
    list(type = "filter", expr = expr)
}

# helper — directed/undirected adjacency slice (x[v_a, v_b]).
# For undirected canonical storage, also matches the swapped orientation.
.edge_slice_op <- function(v_a_int, v_b_int, directed, negate = FALSE) {
    from_id <- to_id <- NULL  # NSE
    expr <- if (isTRUE(directed)) {
        bquote(from_id %in% .(v_a_int) & to_id %in% .(v_b_int))
    } else {
        bquote(
            (from_id %in% .(v_a_int) & to_id %in% .(v_b_int)) |
            (from_id %in% .(v_b_int) & to_id %in% .(v_a_int))
        )
    }
    if (isTRUE(negate)) expr <- bquote(!(.(expr)))
    list(type = "filter", expr = expr)
}

# helper — single-endpoint filter (from = ... or to = ... alone).
# For undirected, "edge has v as an endpoint" — match either side.
.edge_endpoint_op <- function(v_int, side, directed, negate = FALSE) {
    from_id <- to_id <- NULL  # NSE
    expr <- if (!isTRUE(directed)) {
        bquote(from_id %in% .(v_int) | to_id %in% .(v_int))
    } else if (side == "from") {
        bquote(from_id %in% .(v_int))
    } else {
        bquote(to_id %in% .(v_int))
    }
    if (isTRUE(negate)) expr <- bquote(!(.(expr)))
    list(type = "filter", expr = expr)
}

# --- single-arg forms: x[v_set] = induced subgraph ----------------

setMethod("[",
    signature(x = "parquetEdgeStore", i = "character",
              j = "missing", drop = "missing"),
    function(x, i, j, ..., negate = FALSE, drop) {
        v_int <- .edge_vset_to_int(x, i)
        x@ops <- c(x@ops, list(.edge_induced_op(v_int, negate = negate)))
        x
    }
)

setMethod("[",
    signature(x = "parquetEdgeStore", i = "numeric",
              j = "missing", drop = "missing"),
    function(x, i, j, ..., negate = FALSE, drop) {
        v_int <- .edge_vset_to_int(x, i)
        x@ops <- c(x@ops, list(.edge_induced_op(v_int, negate = negate)))
        x
    }
)

# --- two-arg forms: x[i, j] = from/to slice -----------------------
# Register all four char/num combinations so dispatch is unambiguous.

.parquetedge_slice_method <- function(x, i, j, ..., negate = FALSE, drop) {
    v_a <- .edge_vset_to_int(x, i)
    v_b <- .edge_vset_to_int(x, j)
    x@ops <- c(x@ops, list(.edge_slice_op(
        v_a, v_b, directed = x@directed, negate = negate
    )))
    x
}

setMethod("[",
    signature(x = "parquetEdgeStore", i = "character",
              j = "character", drop = "missing"),
    .parquetedge_slice_method)
setMethod("[",
    signature(x = "parquetEdgeStore", i = "numeric",
              j = "numeric", drop = "missing"),
    .parquetedge_slice_method)
setMethod("[",
    signature(x = "parquetEdgeStore", i = "character",
              j = "numeric", drop = "missing"),
    .parquetedge_slice_method)
setMethod("[",
    signature(x = "parquetEdgeStore", i = "numeric",
              j = "character", drop = "missing"),
    .parquetedge_slice_method)

# --- x[, j] = to-endpoint filter ----------------------------------

.parquetedge_to_method <- function(x, i, j, ..., negate = FALSE, drop) {
    v_int <- .edge_vset_to_int(x, j)
    x@ops <- c(x@ops, list(.edge_endpoint_op(
        v_int, side = "to", directed = x@directed, negate = negate
    )))
    x
}

setMethod("[",
    signature(x = "parquetEdgeStore", i = "missing",
              j = "character", drop = "missing"),
    .parquetedge_to_method)
setMethod("[",
    signature(x = "parquetEdgeStore", i = "missing",
              j = "numeric", drop = "missing"),
    .parquetedge_to_method)

# --- x[from = ..., to = ...] = explicit named form ----------------
# Same observable behavior as x[i, j], routed through named args for
# callers who prefer the explicit shape.

setMethod("[",
    signature(x = "parquetEdgeStore", i = "missing",
              j = "missing", drop = "missing"),
    function(x, i, j, ..., from = NULL, to = NULL,
             negate = FALSE, drop) {
        f_int <- .edge_vset_to_int(x, from)
        t_int <- .edge_vset_to_int(x, to)
        if (is.null(f_int) && is.null(t_int)) return(x)
        if (!is.null(f_int) && !is.null(t_int)) {
            x@ops <- c(x@ops, list(.edge_slice_op(
                f_int, t_int, directed = x@directed, negate = negate
            )))
        } else if (!is.null(f_int)) {
            x@ops <- c(x@ops, list(.edge_endpoint_op(
                f_int, side = "from", directed = x@directed,
                negate = negate
            )))
        } else {
            x@ops <- c(x@ops, list(.edge_endpoint_op(
                t_int, side = "to", directed = x@directed,
                negate = negate
            )))
        }
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