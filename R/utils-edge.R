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
# `weighted = TRUE`, else 1), one row per pair of cells. `map` is
# `.edge_to_cell_map(x, pe)`.
#
# A directed store (a kNN network) can hold a mutual pair in both directions;
# those collapse to one canonical row whose weight is the mean of the two, the
# symmetric part of the weights. An undirected store is canonical already.
.edge_view_edges <- function(x, map, weighted = FALSE) {
    from_id <- to_id <- a <- b <- weight <- NULL  # NSE
    nodes <- arrow::arrow_table(int_id = as.integer(map$int_id))
    w <- if (isTRUE(weighted)) rlang::sym("weight") else 1.0
    q <- storeRead(x, output = "arrow") |>
        dplyr::transmute(from_id = from_id, to_id = to_id, weight = !!w) |>
        dplyr::semi_join(nodes, by = c(from_id = "int_id")) |>
        dplyr::semi_join(nodes, by = c(to_id = "int_id"))
    if (!isTRUE(x@directed)) return(q)
    q |>
        dplyr::transmute(
            a = dplyr::if_else(from_id <= to_id, from_id, to_id),
            b = dplyr::if_else(from_id <= to_id, to_id, from_id),
            weight = weight) |>
        dplyr::group_by(a, b) |>
        dplyr::summarise(weight = mean(weight), .groups = "drop") |>
        dplyr::select(from_id = a, to_id = b, weight)
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

# Per-feature sum over the pairs of cells joined by `x` of `w * x_i * x_j`,
# where x is the value the store's chain produces in the two cells and `w` the
# pair's weight (`weighted = TRUE`) or 1 -- each pair once, as
# `.edge_view_edges()` gives them. Returns a numeric vector in the store's
# feature view order. An edge with an end outside the view, or at a node with
# no cell, contributes nothing.
#
# The cross-cell term of every edge statistic (Moran's I, Geary's C, join
# counts) -- the per-cell terms are ordinary accumulators. The pairs are
# collected once as view positions (~16 B a pair); the values are joined a
# cell window at a time: the window's triplets at the `from` end, the
# triplets of the window's neighbours at the `to` end, matched on the
# feature, summed and folded. The joins build on what they read, so the peak
# is the window plus its neighbours rather than the whole view. With one
# window -- the budget covers the view -- the `to` side is the whole view and
# this is a single plan.
#
# On a union a feature's on-disk `col_id` is per substore, so the two ends are
# matched on the feature's view position instead. A chain with post ops cannot
# lower into the plan, so it is refused rather than materialized whole.
.pe_edge_products <- function(pe, x, weighted = FALSE,
    map = .edge_to_cell_map(x, pe)) {
    from_id <- to_id <- weight <- a <- NULL  # NSE
    if (length(pe@post_ops) > 0L) {
        stop("[.pe_edge_products] the store's op chain has R-side (post) ",
             "ops, which cannot run inside the edge join. Materialize it ",
             "with storeWrite() first.", call. = FALSE)
    }
    G <- as.integer(pe@n_genes)
    out <- numeric(G)
    ev <- data.table::as.data.table(dplyr::collect(
        .edge_view_edges(x, map, weighted)))
    if (!nrow(ev)) return(out)
    pairs <- data.table::data.table(
        a = map$pos[match(ev$from_id, map$int_id)],
        b = map$pos[match(ev$to_id, map$int_id)],
        w = as.numeric(ev$weight))

    is_union <- inherits(pe, "unionParquetExprStore")
    cell_key <- c(if (is_union) "source_id", "row_id")
    keys <- .pe_axis_pos_map(pe, "cell")
    if (!is_union) keys[, "source_id" := NA_character_]
    keys <- keys[order(keys$pos)]
    fmap <- .pe_axis_pos_map(pe, "feat")
    ftab <- if (is_union) arrow::as_arrow_table(data.frame(
        source_id = as.character(fmap$source_id),
        col_id = as.integer(fmap$key_id),
        fkey = as.integer(fmap$pos)))

    windows <- .pe_windows(pe, .pe_window_cells(pe, bytes_per_nz = 300, k = 0L))
    single <- length(windows) == 1L
    for (d in windows) {
        lo <- d$offset + d$cs
        hi <- d$offset + d$ce
        pw <- pairs[a >= lo & a <= hi]
        if (!nrow(pw)) next
        nb <- sort(unique(pw$b))
        to_store <- if (single) pe else pe[, nb]
        part <- .pe_edge_join(.pe_window_store(d), to_store, pw, keys,
            cell_key, ftab)
        if (nrow(part) == 0L) next
        fpos <- if (is_union) part$fkey else fmap$pos[match(part$fkey, fmap$key_id)]
        ok <- !is.na(fpos)
        out[fpos[ok]] <- out[fpos[ok]] + as.numeric(part$T[ok])
    }
    out
}

# One window of `.pe_edge_products()`: `pw` (view positions a -> b, weight w)
# as cell keys, joined to the triplets of `from_store` at `a` and of
# `to_store` at `b` with the same feature, summed per feature.
.pe_edge_join <- function(from_store, to_store, pw, keys, cell_key, ftab) {
    w <- v1 <- v2 <- fkey <- col_id <- NULL  # NSE
    ek <- function(end, p) {
        k <- keys[p]
        d <- data.frame(row_id = as.integer(k$key_id))
        if ("source_id" %in% cell_key) d$source_id <- as.character(k$source_id)
        names(d) <- paste0(end, "_", names(d))
        d
    }
    etab <- arrow::as_arrow_table(cbind(ek("from", pw$a), ek("to", pw$b),
        data.frame(w = pw$w)))
    trip <- function(store, end, v) {
        q <- storeRead(store, output = "query")
        q <- if (is.null(ftab)) {
            dplyr::mutate(q, fkey = col_id)
        } else {
            dplyr::inner_join(q, ftab, by = c("source_id", "col_id"))
        }
        sel <- c(stats::setNames(cell_key, paste0(end, "_", cell_key)),
                 fkey = "fkey")
        sel[[v]] <- "value"
        dplyr::select(q, !!!rlang::syms(sel))
    }
    etab |>
        dplyr::inner_join(trip(from_store, "from", "v1"),
            by = paste0("from_", cell_key)) |>
        dplyr::inner_join(trip(to_store, "to", "v2"),
            by = c(paste0("to_", cell_key), "fkey")) |>
        dplyr::group_by(fkey) |>
        dplyr::summarise(T = sum(w * v1 * v2, na.rm = TRUE)) |>
        dplyr::collect()
}

# Group tables of the view's edges for a cell grouping. `code` is an integer
# group (1..K) per view cell position, `map` is `.edge_to_cell_map(x, pe)`.
# Returns
#   a      data.table(int_id, k, l, a): cell int_id of group k has summed edge
#          weight a to cells of group l
#   edges  K x K, sum of a over the cells of group k
#   paths  K x K x K, [k, l, m] = sum over cells of group k of a_l * a_m
.edge_group_tables <- function(x, map, code, K, weighted = FALSE) {
    from_id <- int_id <- l <- w <- k <- a <- i.l <- i.a <- v <- m <- NULL  # NSE
    lk <- list(int_id = map$int_id, value = code[map$pos])
    tab <- .edge_symmetrise(.edge_view_edges(x, map, weighted), weight = TRUE) |>
        .edge_join_vertex(lk, by = "to_id", name = "l") |>
        dplyr::group_by(from_id, l) |>
        dplyr::summarise(a = sum(w, na.rm = TRUE)) |>
        dplyr::collect() |>
        data.table::as.data.table()
    data.table::setnames(tab, "from_id", "int_id")
    tab[, k := lk$value[match(int_id, lk$int_id)]]

    edges <- matrix(0, K, K)
    if (nrow(tab)) {
        e <- tab[, list(a = sum(a)), by = c("k", "l")]
        edges[cbind(e$k, e$l)] <- e$a
    }
    paths <- array(0, c(K, K, K))
    if (nrow(tab)) {
        pp <- tab[tab, on = "int_id", allow.cartesian = TRUE,
            list(k, l, m = i.l, v = a * i.a)][
            , list(v = sum(v)), by = c("k", "l", "m")]
        paths[cbind(pp$k, pp$l, pp$m)] <- pp$v
    }
    list(a = tab[, c("int_id", "k", "l", "a"), with = FALSE],
         edges = edges, paths = paths)
}
