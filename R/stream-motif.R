#' @include pkg_imports.R class-parquetEdgeStore.R
NULL

# Motif enrichment on a backed network. Giotto owns the param classes, the
# router and the igraph method; this is the method for the store. Same engine
# (smotif), different carrier: the edges cross as an Arrow stream from
# storeRead(output = "arrowstream"), so a pending subset is honoured and the
# store is never read into an igraph. Whether smotif runs that stream in R or
# hands it to its Rust backend is smotif's decision, not this package's.
#
# Labels arrive named by node ID (Giotto names them by spatIDs()) and are
# realigned by name, because the store's node order is not the order the
# caller read the IDs in.

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetEdgeStore", param = "smotifParam"),
    function(x, param, cell_type = NULL, strata = NULL, anchored_on = NULL,
             ...) {
        if (is.null(cell_type)) {
            stop("[analyzeData] cell_type labels are required", call. = FALSE)
        }
        if (!is.null(strata) || !is.null(anchored_on)) {
            stop("[analyzeData] strata and anchored_on are not supported ",
                "on a backed network yet", call. = FALSE)
        }
        if (isTRUE(param$set_seed)) {
            GiottoUtils::local_seed(seed = param$seed_number)
        }
        package_check("smotif", repository = "github:drieslab/smotif")

        nodes <- .edge_active_nodes(x)
        smotif::motif_enrichment_stream(
            storeRead(x, output = "arrowstream"),
            int_ids = nodes$int_id,
            cell_type = .motif_realign(cell_type, nodes$node_id, "cell_type"),
            size = param$size,
            n_perm = param$n_perm,
            seed = as.integer(param$seed_number),
            null = param$null,
            ...
        )
    }
)


# helpers ####

.motif_realign <- function(lab, ids, what) {
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
# unioned with every node an active edge references -- so the store and the
# igraph method permute labels over the same cells. Ordered by int_id.
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
