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

# The pairs of cells the edges of `x` join within an expression store's view:
# data.table(a, b, w) of view positions (`a < b`) and weight (the stored
# `weight` when `weighted = TRUE`, else 1), one row per pair. `map` is
# `.edge_to_cell_map(x, pe)`. Edges with an end outside the view are dropped.
#
# A directed store (a kNN network) can hold a mutual pair in both directions;
# it becomes one row at the mean of the two weights, the symmetric part. With
# GiottoKernels the edges stream into its graph_pairs_stream kernel and only
# the pairs return to R. Either way the
# edges are read once, as three numeric columns, and everything edge-only
# downstream (degree, weight totals, group tables, the windows of the edge
# products) is derived from this table in R rather than by another pass.
.edge_view_pairs <- function(x, map, weighted = FALSE) {
    a <- b <- w <- NULL  # NSE
    if (.gd_kernel("graph_pairs_stream")) {
        # the edges stream straight into the kernel and only the pairs come
        # back, so the edge list is never materialised in R
        node_row <- integer(max(c(map$int_id, 0L)))
        node_row[map$int_id] <- map$pos
        p <- GiottoKernels::graph_pairs_stream(
            storeRead(x, output = "arrowstream"), node_row,
            weighted = isTRUE(weighted), dedupe = isTRUE(x@directed))
        return(data.table::data.table(a = p$a, b = p$b, w = p$w))
    }
    cols <- c("from_id", "to_id", if (isTRUE(weighted)) "weight")
    ev <- storeRead(x, output = "arrow") |>
        dplyr::select(dplyr::all_of(cols)) |>
        dplyr::collect()
    pos <- integer(max(c(map$int_id, 0L)))
    pos[map$int_id] <- map$pos
    lookup <- function(id) {
        id <- as.integer(id)
        out <- rep(NA_integer_, length(id))
        ok <- !is.na(id) & id >= 1L & id <= length(pos)
        out[ok] <- pos[id[ok]]
        out[out == 0L] <- NA_integer_
        out
    }
    pa <- lookup(ev$from_id)
    pb <- lookup(ev$to_id)
    pairs <- data.table::data.table(a = pmin(pa, pb), b = pmax(pa, pb),
        w = if (isTRUE(weighted)) as.numeric(ev$weight) else 1)
    rm(ev, pa, pb)
    pairs <- pairs[!is.na(a) & !is.na(b)]
    if (isTRUE(x@directed)) pairs <- pairs[, list(w = mean(w)), by = c("a", "b")]
    pairs[]
}

# Weighted degree of every view cell from `.edge_view_pairs()`, 0 for none.
.pairs_degree <- function(pairs, n_cells) {
    deg <- numeric(n_cells)
    if (!nrow(pairs)) return(deg)
    d <- rowsum(c(pairs$w, pairs$w), c(pairs$a, pairs$b))
    deg[as.integer(rownames(d))] <- d[, 1L]
    deg
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
# `.edge_view_pairs()` gives them. Returns a numeric vector in the store's
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
# Row-range windows need a local stored order, or each window's neighbours
# span the store. When the order is not local and the windows are many
# (`giottodisk.edge_rewrite_min_windows`, default 32), the view is first
# rewritten, out of core, in breadth-first network order into a temporary
# store (.pe_rewrite_cell_order; temp disk about the size of the view's
# values, removed on exit), and the windows run on that. With fewer windows,
# breadth-first windows read by filtered scans are cheaper than the rewrite.
#
# On a union a feature's on-disk `col_id` is per substore, so the two ends are
# matched on the feature's view position instead. A chain with post ops cannot
# lower into the plan, so it is refused rather than materialized whole.
.pe_edge_products <- function(pe, x, weighted = FALSE,
    pairs = .edge_view_pairs(x, .edge_to_cell_map(x, pe), weighted),
    rewritten = FALSE) {
    a <- NULL  # NSE
    if (length(pe@post_ops) > 0L) {
        stop("[.pe_edge_products] the store's op chain has R-side (post) ",
             "ops, which cannot run inside the edge join. Materialize it ",
             "with storeWrite() first.", call. = FALSE)
    }
    G <- as.integer(pe@n_genes)
    out <- numeric(G)
    if (!nrow(pairs)) return(out)
    kern <- .pe_edge_products_kernel(pe, pairs)
    if (!is.null(kern)) return(kern)

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

    n_cells <- as.integer(pe@n_cells)
    budget <- .pe_window_cells(pe, bytes_per_nz = 300, k = 0L)
    # The rewrite is a fixed cost of a few passes; scanning per window costs
    # two scans a window. It pays only when the windows are many (a small
    # machine); with a few large windows the scans are cheaper.
    n_win <- ceiling(n_cells / budget)
    min_win <- getOption("giottodisk.edge_rewrite_min_windows", 32L)
    if (!rewritten && n_win > min_win && !.pairs_local(pairs, budget)) {
        # stored order is not local: rewrite the view once in network order,
        # where the windows are contiguous row ranges with a thin rim (and
        # never again -- a window narrower than a BFS layer keeps failing the
        # locality test, but the order cannot get more local than this)
        ord <- .pairs_bfs_order(pairs, n_cells)
        tmp <- .dump_tempfile()
        on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
        sorted <- .pe_rewrite_cell_order(pe, ord, tmp)
        np <- integer(n_cells)
        np[ord] <- seq_len(n_cells)
        a2 <- np[pairs$a]; b2 <- np[pairs$b]
        moved <- data.table::data.table(a = pmin(a2, b2), b = pmax(a2, b2), w = pairs$w)
        return(.pe_edge_products(sorted, x, pairs = moved, rewritten = TRUE))
    }
    for (win in .pe_edge_windows(pe, budget, n_cells,
        if (!rewritten && !.pairs_local(pairs, budget)) pairs)) {
        pw <- pairs[a %in% win$from]
        if (!nrow(pw)) next
        part <- .pe_edge_join(win$store, if (win$single) pe else pe[, sort(unique(pw$b))],
            pw, keys, cell_key, ftab)
        if (nrow(part) == 0L) next
        fpos <- if (is_union) part$fkey else fmap$pos[match(part$fkey, fmap$key_id)]
        ok <- !is.na(fpos)
        out[fpos[ok]] <- out[fpos[ok]] + as.numeric(part$T[ok])
    }
    out
}

# `.pe_edge_products()` through GiottoKernels' edge_products_stream, when it
# is installed and the store is a single one; NULL means "use the R path".
# The kernel holds the triplets it is given in memory (about 16 B a stored
# value at peak) with random access over every cell, so cell order does not
# matter and nothing is windowed or rewritten. Memory is bounded instead by
# streaming the features in blocks -- their sums are independent -- each block
# sized from the stores' per-feature counts to fit the chunk budget, one
# filtered scan per block.
.pe_edge_products_kernel <- function(pe, pairs) {
    if (inherits(pe, "unionParquetExprStore") ||
        !.gd_kernel("edge_products_stream")) {
        return(NULL)
    }
    cells <- .pe_axis_pos_map(pe, "cell")
    kid <- integer(as.integer(pe@n_cells))
    kid[cells$pos] <- cells$key_id
    feats <- .pe_axis_pos_map(pe, "feat")
    n_rows <- max(kid)
    n_features <- max(feats$key_id)
    w <- if (all(pairs$w == 1)) NULL else pairs$w
    from <- kid[pairs$a]
    to <- kid[pairs$b]

    # feature blocks under the budget, in view order
    nnz <- pe@stats$col_nnz
    per <- if (length(nnz) >= n_features) as.numeric(nnz[feats$key_id]) else
        rep(as.numeric(pe@n_cells), nrow(feats))
    free <- tryCatch(.sc_free_ram(), error = function(e) NA_real_)
    budget <- if (isTRUE(is.finite(free)) && free > 0) {
        .chunk_budget(free, .ram_frac(),
            n_cells = 0, k = 0, p = 0)
    } else 2e9
    pinned <- getOption("giottodisk.edge_kernel_block_nnz")
    # half the budget for the kernel's table (16 B a value at peak), half for
    # the arrow scan feeding it, whose buffers do not shrink with the block
    cap <- if (!is.null(pinned)) as.numeric(pinned) else max(1, budget / 2 / 16)
    block <- cumsum(per) %/% cap
    out <- numeric(as.integer(pe@n_genes))
    for (b in unique(block)) {
        sel <- which(block == b)
        reader <- .gd_kernel_reader(pe[feats$pos[sel], ])
        tot <- GiottoKernels::edge_products_stream(reader,
            from = from, to = to,
            n_rows = n_rows, n_features = n_features, w = w,
            n_threads = .gd_kernel_threads())
        out[feats$pos[sel]] <- tot[feats$key_id[sel]]
        # a finished reader keeps its scan's buffers until it is collected,
        # so release it before the next block's scan starts
        rm(reader)
        invisible(gc(verbose = FALSE))
    }
    out
}

# The cell windows `.pe_edge_products()` walks: list(store, from, single) with
# `from` the window's view positions and `store` the store narrowed to them.
# One window when the budget covers the view, otherwise row ranges
# (`.pe_windows()`, which prune row groups). Given `pairs` (stored order not
# local, too few windows to be worth a rewrite), runs of a breadth-first order
# of the network instead: compact patches whose cells and neighbours are read
# by a filtered scan each.
.pe_edge_windows <- function(pe, budget, n_cells, pairs = NULL) {
    if (budget >= n_cells) {
        return(list(list(store = pe, from = seq_len(n_cells), single = TRUE)))
    }
    if (is.null(pairs)) {
        return(lapply(.pe_windows(pe, budget), function(d) list(
            store = .pe_window_store(d),
            from = (d$offset + d$cs):(d$offset + d$ce), single = FALSE)))
    }
    ord <- .pairs_bfs_order(pairs, n_cells)
    lapply(split(ord, ceiling(seq_along(ord) / budget)), function(w) {
        w <- sort(w)
        list(store = pe[, w], from = w, single = FALSE)
    })
}

# Is the stored cell order local enough for row-range windows of `budget`
# cells? Nearly every pair must join cells less than a window apart: a few far
# pairs already scatter a window's neighbours over the whole store, so its
# neighbour read can no longer prune.
.pairs_local <- function(pairs, budget) {
    mean(abs(pairs$a - pairs$b) < budget) >= 0.99
}

# A breadth-first order of the network's cells (view positions), unreached
# cells last: consecutive runs are compact patches of tissue.
.pairs_bfs_order <- function(pairs, n_cells) {
    g <- igraph::graph_from_edgelist(cbind(pairs$a, pairs$b), directed = FALSE)
    g <- igraph::add_vertices(g, max(0L, n_cells - igraph::vcount(g)))
    ord <- as.integer(igraph::bfs(g, root = 1L, unreachable = TRUE)$order)
    ord <- ord[!is.na(ord)]
    c(ord, setdiff(seq_len(n_cells), ord))
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

# Group tables of the view's pairs for a cell grouping. `code` is an integer
# group (1..K) per view cell position, `pairs` is `.edge_view_pairs()`.
# Returns
#   a      data.table(pos, k, l, a): the cell at view position `pos`, of group
#          k, has summed edge weight a to cells of group l
#   edges  K x K, sum of a over the cells of group k
#   paths  K x K x K, [k, l, m] = sum over cells of group k of a_l * a_m
.edge_group_tables <- function(pairs, code, K) {
    n_cells <- length(code)
    # a_il as a sparse cells x groups matrix (duplicates sum); every table
    # below is a sparse product of it, never a join of cells with themselves
    A <- Matrix::sparseMatrix(
        i = c(pairs$a, pairs$b), j = code[c(pairs$b, pairs$a)],
        x = c(pairs$w, pairs$w), dims = c(n_cells, K))
    Z <- Matrix::sparseMatrix(i = seq_len(n_cells), j = code, x = 1,
        dims = c(n_cells, K))
    edges <- as.matrix(Matrix::crossprod(Z, A))
    paths <- array(0, c(K, K, K))
    for (k in seq_len(K)) {
        sel <- which(code == k)
        if (length(sel)) {
            paths[k, , ] <- as.matrix(Matrix::crossprod(A[sel, , drop = FALSE]))
        }
    }
    s <- Matrix::summary(A)
    tab <- data.table::data.table(pos = s$i, k = code[s$i], l = s$j, a = s$x)
    list(a = tab, edges = unname(edges), paths = paths)
}
