#' @include pkg_imports.R class-parquetEdgeStore.R
NULL

# Motif enrichment on a backed network. Giotto owns the param classes, the
# router and the igraph method; this is the method for the store. Same engine
# (smotif), different carrier: the edges cross as an Arrow stream from
# storeRead(output = "arrowstream"), so a pending subset is honoured and the
# store is never read into an igraph. Whether smotif runs that stream in R or
# hands it to its Rust backend is smotif's decision, not this package's.
#
# The node sidecar is this method's concern only. Labels arrive named by node
# ID (Giotto names them by spatIDs()); the sidecar turns them into a label per
# int_id, and only that lookup crosses with the edge stream. The engine takes
# its node set from the edges it receives, so a pending subset decides which
# cells are in the network and the lookup may cover more than that.

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

        lookup <- .edge_vertex_lookup(x, cell_type)
        smotif::motif_enrichment_stream(
            storeRead(x, output = "arrowstream"),
            int_ids = lookup$int_id,
            cell_type = lookup$value,
            size = param$size,
            n_perm = param$n_perm,
            seed = as.integer(param$seed_number),
            null = param$null,
            ...
        )
    }
)

