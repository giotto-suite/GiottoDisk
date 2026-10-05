#' @include pkg_imports.R class-parquetEdgeStore.R
NULL

# Motif enrichment on a backed network. Giotto owns the router and the
# igraph-side engine (smotifParam); this is the half that owns the store.
#
# smotifrsParam is a class rather than a knob because the method body differs:
# the edges cross into smotifrs as an Arrow stream from
# storeRead(output = "arrowstream"), so a pending subset is honoured and the
# store is never read into an igraph. Same arrangement as gramEigenPcaParam:
# the substrate selects the engine when "auto" resolves.
#
# Labels arrive named by node ID (Giotto names them by spatIDs()), and every
# path here realigns by name, because the store's node order is not the order
# the caller read the IDs in.


#' @rdname smotifrsParam
#' @exportClass smotifrsParam
setClass("smotifrsParam", contains = "motifParam")


#' @name smotifrsParam
#' @title Streaming motif enrichment parameter
#' @description
#' Motif enrichment on a backed spatial network that streams its edges into
#' \pkg{smotifrs} rather than reading the network into an `igraph`. This is
#' what `method = "auto"` resolves to on a `parquetEdgeStore` when
#' \pkg{smotifrs} can read a stream; it rarely needs to be built by hand.
#'
#' The stream counts over the whole network under a `"label"` or
#' `"conditional"` null. A `"stratified"` null, or a call with `strata` or
#' `anchored_on`, is a smotif feature on an in-memory graph: for those the
#' store is read into an `igraph` and the smotif engine runs instead.
#' @param size motif size: 2, 3 or 4.
#' @param null null model: `"label"`, `"stratified"` or `"conditional"`.
#' @param n_perm number of null draws.
#' @param set_seed,seed_number seed control.
#' @param ... further engine arguments, e.g. `cond_temp`.
#' @returns A `smotifrsParam` object.
#' @seealso [storeRead()] for the `"arrowstream"` output the edges cross as.
#' @examples
#' # p <- smotifrsParam(size = 3L, null = "conditional")
#' # res <- analyzeData(store, p, cell_type = labels_named_by_node_id)
#' @export
smotifrsParam <- function(size = 3L,
                          null = c("label", "stratified", "conditional"),
                          n_perm = 1000L,
                          set_seed = TRUE,
                          seed_number = 1234,
                          ...) {
    # Giotto's factory owns validation of the shared settings
    p <- Giotto::motifParam(
        size = size, null = match.arg(null), n_perm = n_perm,
        set_seed = set_seed, seed_number = seed_number, ...
    )
    methods::new("smotifrsParam", param = p@param)
}


# autoMotifParam on parquetEdgeStore: stream when smotifrs can read one,
# otherwise read the store into an igraph and let Giotto resolve "auto" there.

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetEdgeStore", param = "autoMotifParam"),
    function(x, param, ...) {
        if (.smotifrs_streams()) {
            p <- methods::new("smotifrsParam", param = param@param)
            return(analyzeData(x, p, ...))
        }
        .motif_via_igraph(x, param, ...)
    }
)


#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetEdgeStore", param = "smotifrsParam"),
    function(x, param, cell_type = NULL, strata = NULL, anchored_on = NULL,
             ...) {
        if (is.null(cell_type)) {
            stop("[analyzeData] cell_type labels are required", call. = FALSE)
        }
        if (!is.null(strata) || !is.null(anchored_on) ||
            !param$null %in% c("label", "conditional")) {
            p <- methods::new("smotifParam", param = param@param)
            return(.motif_via_igraph(x, p,
                cell_type = cell_type, strata = strata,
                anchored_on = anchored_on, ...
            ))
        }
        if (!.smotifrs_streams()) {
            stop("[analyzeData] smotifrsParam needs a smotifrs that ",
                "provides motif_enrichment_stream()", call. = FALSE)
        }

        nodes <- .edge_active_nodes(x)
        lab <- .motif_realign(cell_type, nodes$node_id, "cell_type")
        dots <- list(...)
        stream_fn <- getExportedValue("smotifrs", "motif_enrichment_stream")
        stream_fn(
            storeRead(x, output = "arrowstream"),
            int_ids = nodes$int_id,
            cell_type = lab,
            size = param$size,
            n_perm = param$n_perm,
            seed = as.integer(param$seed_number),
            null = param$null,
            cond_temp = dots$cond_temp %null% param$cond_temp %null% 1
        )
    }
)


# helpers ####

.smotifrs_streams <- function() {
    requireNamespace("smotifrs", quietly = TRUE) &&
        "motif_enrichment_stream" %in% getNamespaceExports("smotifrs")
}

# Read the store into an igraph and dispatch there, realigning the named
# labels to the graph's vertex order.
.motif_via_igraph <- function(x, param, cell_type = NULL, strata = NULL, ...) {
    ig <- igraph::as.igraph(x)
    v <- igraph::V(ig)$name
    analyzeData(ig, param,
        cell_type = .motif_realign(cell_type, v, "cell_type"),
        strata = .motif_realign(strata, v, "strata"),
        ...
    )
}

.motif_realign <- function(lab, ids, what) {
    if (is.null(lab)) return(NULL)
    if (is.null(names(lab))) {
        stop(sprintf("[analyzeData] %s must be named by node ID", what),
            call. = FALSE)
    }
    out <- lab[ids]
    if (anyNA(out)) {
        stop(sprintf("[analyzeData] %d network node(s) have no %s label",
            sum(is.na(out)), what), call. = FALSE)
    }
    unname(out)
}

# The vertex set as.igraph() builds for this store -- the recorded selection
# unioned with every node an active edge references -- so the stream and the
# igraph fallback permute labels over the same cells. Ordered by int_id.
.edge_active_nodes <- function(x) {
    int_id <- node_id <- from_id <- to_id <- NULL  # NSE
    edges <- storeRead(x, output = "arrow")
    used <- dplyr::union(
        edges |> dplyr::select(int_id = from_id),
        edges |> dplyr::select(int_id = to_id)
    ) |>
        dplyr::collect() |>
        dplyr::pull(int_id)
    used <- sort(unique(c(as.integer(x@node_idx), as.integer(used))))
    nm <- storeRead(x@nodes) |>
        dplyr::filter(int_id %in% !!used) |>
        dplyr::select(int_id, node_id) |>
        dplyr::collect() |>
        data.table::as.data.table()
    nm <- nm[match(used, nm$int_id)]
    data.table::data.table(
        node_id = as.character(nm$node_id),
        int_id = as.integer(nm$int_id)
    )
}
