#' @include class-parquetEdgeStore.R
NULL

# Edge-store building blocks shared by the methods that consume a network
# (label proportions, motifs, and statistics over edges). Each returns either a
# lazy arrow query or a small R lookup; none collects the edge set.
#
# Edges are in int_id space throughout. Per-node values arrive named by node
# ID and cross into that space once, through the node sidecar.

# Both orientations of every stored edge, as one lazy query with columns
# `from_id`, `to_id`, `w`. An undirected store keeps each edge once
# (`from_id <= to_id`), so a per-vertex aggregate grouped on `from_id` needs
# this to see every neighbour. `w` is the stored `weight` column when
# `weight = TRUE`, else 1. Two scans of the same dataset, nothing materialized.
.edge_symmetrise <- function(edges, weight = FALSE) {
    from_id <- to_id <- f <- t <- NULL  # NSE
    w <- if (isTRUE(weight)) rlang::sym("weight") else 1.0
    forward <- edges |>
        dplyr::transmute(from_id = from_id, to_id = to_id, w = !!w)
    # select() renames simultaneously; a transmute would see its own new
    # from_id when computing to_id
    backward <- edges |>
        dplyr::transmute(f = to_id, t = from_id, w = !!w) |>
        dplyr::select(from_id = f, to_id = t, w)
    dplyr::union_all(forward, backward)
}

# One self edge of weight `w` per node in the sidecar, shaped like
# `.edge_symmetrise()` so the two can be unioned.
.edge_self_loops <- function(x, w) {
    int_id <- NULL  # NSE
    storeRead(x@nodes) |>
        dplyr::transmute(from_id = int_id, to_id = int_id,
                         w = !!as.double(w))
}

# Per-node values named by node ID -> one value per int_id, through the node
# sidecar. The sidecar is one row per cell, so it is read whole and matched in
# R rather than filtered in Arrow against a literal set of every name. Nodes
# without a value are left out; what an edge reaching one means is the
# caller's decision.
.edge_vertex_lookup <- function(x, values) {
    int_id <- node_id <- NULL  # NSE
    if (is.null(names(values))) {
        stop("[analyzeData] per-node values must be named by node ID",
            call. = FALSE)
    }
    nm <- storeRead(x@nodes) |>
        dplyr::select(int_id, node_id) |>
        dplyr::collect()
    i <- match(as.character(nm$node_id), names(values))
    keep <- !is.na(i)
    list(
        int_id = as.integer(nm$int_id[keep]),
        value = unname(values[i[keep]])
    )
}

# Attach a `.edge_vertex_lookup()` result to one endpoint of a lazy edge query,
# as column `name`. A left join: edges whose endpoint has no value keep NA.
.edge_join_vertex <- function(edges, lookup, by = "to_id", name = "value") {
    cols <- list(lookup$int_id, lookup$value)
    names(cols) <- c(by, name)
    dplyr::left_join(edges, do.call(arrow::arrow_table, cols), by = by)
}
