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

# Degree of every node with at least one edge in the lazy edge query `edges`
# (as `storeRead(x, output = "arrow")` gives, or narrowed from it):
# data.table(int_id, deg). Each edge counts once at both ends -- on a directed
# store that is in-degree plus out-degree. `weighted = TRUE` sums the `weight`
# column instead of counting. Nodes with no edge are absent; the caller decides
# they are 0.
.edge_degree <- function(edges, weighted = FALSE) {
    from_id <- w <- NULL  # NSE
    .edge_symmetrise(edges, weight = weighted) |>
        dplyr::group_by(from_id) |>
        dplyr::summarise(deg = sum(w, na.rm = TRUE)) |>
        dplyr::collect() |>
        data.table::as.data.table() |>
        data.table::setnames("from_id", "int_id")
}

# The edges of `x` with both ends at cells of the expression store's view, as a
# lazy query with `from_id`, `to_id` and `weight` (the stored weight when
# `weighted = TRUE`, else 1). `map` is `.edge_to_cell_map(x, pe)`.
.edge_view_edges <- function(x, map, weighted = FALSE) {
    from_id <- to_id <- NULL  # NSE
    nodes <- arrow::arrow_table(int_id = as.integer(map$int_id))
    w <- if (isTRUE(weighted)) rlang::sym("weight") else 1.0
    storeRead(x, output = "arrow") |>
        dplyr::transmute(from_id = from_id, to_id = to_id, weight = !!w) |>
        dplyr::semi_join(nodes, by = c(from_id = "int_id")) |>
        dplyr::semi_join(nodes, by = c(to_id = "int_id"))
}

# Network nodes -> cells of an expression store, matched by ID:
# data.table(int_id, pos, row_id[, source_id]). `pos` is the cell's position in
# the store's current view and `row_id` its ON-DISK id (ADR 0003), so a caller
# can key a payload or join edges onto the triplets. A node with no cell in the
# view is absent, as is a cell with no node.
.edge_to_cell_map <- function(x, pe) {
    pos <- NULL  # NSE
    lk <- .edge_vertex_lookup(x,
        stats::setNames(seq_along(pe@cell_ids), pe@cell_ids))
    nodes <- data.table::data.table(int_id = lk$int_id, pos = lk$value)
    out <- merge(nodes, .pe_axis_pos_map(pe, "cell"), by = "pos")
    data.table::setnames(out, "key_id", "row_id")
    data.table::setcolorder(out, c("int_id", "pos", "row_id"))
    data.table::setorder(out, pos)
    out[]
}

# Per-feature sum over the stored edges of `w * x_i * x_j`, where x is the
# value the store's chain produces in the cells at the two ends of an edge and
# `w` the edge's `weight` (`weighted = TRUE`) or 1. Returns a numeric vector in
# the store's feature view order. An undirected store keeps each edge once, so
# this is the sum over undirected edges; an edge with an end outside the view,
# or at a node with no cell, contributes nothing.
#
# The cross-cell term of every edge statistic (Moran's I, Geary's C, join
# counts) -- the per-cell terms are ordinary accumulators. One arrow plan:
# edges -> cell keys of both ends -> triplets at `from` -> triplets at `to`
# with the same feature -> grouped sum. Nothing is collected but the result.
# The plan's hash joins build on the triplets, so peak memory is of the order
# of the stored entries in view.
#
# On a union a feature's on-disk `col_id` is per substore, so the two ends are
# matched on the feature's view position instead. A chain with post ops cannot
# lower into the plan, so it is refused rather than materialized whole.
.pe_edge_products <- function(pe, x, weighted = FALSE) {
    from_id <- to_id <- w <- v1 <- v2 <- fkey <- col_id <- NULL  # NSE
    if (length(pe@post_ops) > 0L) {
        stop("[.pe_edge_products] the store's op chain has R-side (post) ",
             "ops, which cannot run inside the edge join. Materialize it ",
             "with storeWrite() first.", call. = FALSE)
    }
    is_union <- inherits(pe, "unionParquetExprStore")
    cell_key <- c(if (is_union) "source_id", "row_id")

    # node -> cell key, once per end
    cmap <- .edge_to_cell_map(x, pe)
    end_tab <- function(end) {
        d <- as.data.frame(cmap[, c("int_id", cell_key), with = FALSE])
        names(d) <- c(paste0(end, "_id"), paste0(end, "_", cell_key))
        arrow::as_arrow_table(d)
    }

    trip <- storeRead(pe, output = "query")
    if (is_union) {
        fmap <- .pe_axis_pos_map(pe, "feat")
        trip <- trip |>
            dplyr::inner_join(arrow::as_arrow_table(data.frame(
                source_id = as.character(fmap$source_id),
                col_id = as.integer(fmap$key_id),
                fkey = as.integer(fmap$pos))), by = c("source_id", "col_id"))
    } else {
        trip <- dplyr::mutate(trip, fkey = col_id)
    }
    trip_end <- function(end, v) {
        sel <- c(stats::setNames(cell_key, paste0(end, "_", cell_key)),
                 fkey = "fkey")
        sel[[v]] <- "value"
        dplyr::select(trip, !!!rlang::syms(sel))
    }

    w_expr <- if (isTRUE(weighted)) rlang::sym("weight") else 1.0
    res <- storeRead(x, output = "arrow") |>
        dplyr::transmute(from_id = from_id, to_id = to_id, w = !!w_expr) |>
        dplyr::inner_join(end_tab("from"), by = "from_id") |>
        dplyr::inner_join(end_tab("to"), by = "to_id") |>
        dplyr::inner_join(trip_end("from", "v1"),
            by = paste0("from_", cell_key)) |>
        dplyr::inner_join(trip_end("to", "v2"),
            by = c(paste0("to_", cell_key), "fkey")) |>
        dplyr::group_by(fkey) |>
        dplyr::summarise(T = sum(w * v1 * v2, na.rm = TRUE)) |>
        dplyr::collect()

    out <- numeric(as.integer(pe@n_genes))
    if (nrow(res) == 0L) return(out)
    pos <- if (is_union) {
        res$fkey
    } else {
        fmap <- .pe_axis_pos_map(pe, "feat")
        fmap$pos[match(res$fkey, fmap$key_id)]
    }
    keep <- !is.na(pos)
    out[pos[keep]] <- as.numeric(res$T[keep])
    out
}
