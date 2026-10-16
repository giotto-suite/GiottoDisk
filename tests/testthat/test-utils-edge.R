# Edge-store building blocks (utils-edge.R), checked against the edge list the
# store was written from.

.edge_fixture <- function() {
    dt <- data.table::data.table(
        from = c("a", "a", "b", "c"), to = c("b", "c", "c", "d"),
        weight = c(0.5, 1.5, 2, 3)
    )
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), dt,
        type = "sNN", directed = FALSE)
    list(store = s, dt = dt)
}

.edge_collect_named <- function(q, s) {
    nm <- dplyr::collect(storeRead(s@nodes))
    out <- data.table::as.data.table(dplyr::collect(q))
    out[, from := nm$node_id[match(from_id, nm$int_id)]]
    out[, to := nm$node_id[match(to_id, nm$int_id)]]
    out
}

test_that(".edge_symmetrise returns both orientations of every edge", {
    f <- .edge_fixture()
    got <- .edge_collect_named(
        .edge_symmetrise(storeRead(f$store, output = "arrow")), f$store)
    expect_equal(nrow(got), 2L * nrow(f$dt))
    want <- rbind(f$dt[, .(from, to)], f$dt[, .(from = to, to = from)])
    expect_setequal(paste(got$from, got$to), paste(want$from, want$to))
    expect_true(all(got$w == 1))
})

test_that(".edge_symmetrise carries the stored weight on both orientations", {
    f <- .edge_fixture()
    got <- .edge_collect_named(
        .edge_symmetrise(storeRead(f$store, output = "arrow"), weight = TRUE),
        f$store)
    key <- function(a, b) paste(pmin(a, b), pmax(a, b))
    w_want <- f$dt$weight[match(key(got$from, got$to), key(f$dt$from, f$dt$to))]
    expect_equal(got$w, w_want)
})

test_that(".edge_self_loops gives one self edge per node", {
    f <- .edge_fixture()
    got <- dplyr::collect(.edge_self_loops(f$store, 0.25))
    expect_equal(nrow(got), 4L)
    expect_equal(got$from_id, got$to_id)
    expect_true(all(got$w == 0.25))
})

test_that(".edge_join_vertex attaches the value of the chosen endpoint", {
    f <- .edge_fixture()
    vals <- c(a = 10, b = 20, d = 40)  # c has no value
    lk <- .edge_vertex_lookup(f$store, vals)
    expect_length(lk$int_id, 3L)
    q <- .edge_join_vertex(storeRead(f$store, output = "arrow"), lk,
        by = "to_id", name = "v")
    got <- .edge_collect_named(q, f$store)
    expect_equal(got$v, unname(vals[got$to]))  # NA where `to` is c
})

test_that(".edge_degree matches igraph, counted and weighted", {
    f <- .edge_fixture()
    g <- igraph::graph_from_data_frame(f$dt, directed = FALSE)
    nm <- dplyr::collect(storeRead(f$store@nodes))
    name_of <- function(d) nm$node_id[match(d$int_id, nm$int_id)]

    d <- .edge_degree(f$store)
    expect_equal(d$deg, unname(igraph::degree(g)[name_of(d)]))
    dw <- .edge_degree(f$store, weighted = TRUE)
    expect_equal(dw$deg, unname(igraph::strength(g)[name_of(dw)]))
})

test_that(".edge_to_cell_map matches nodes to cells by ID, in view order", {
    f <- .edge_fixture()
    m <- Matrix::sparseMatrix(i = c(1, 2, 1), j = c(1, 3, 5), x = 1,
        dims = c(2, 5), dimnames = list(c("g1", "g2"),
                                        c("e", "d", "c", "b", "a")))
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    nm <- dplyr::collect(storeRead(f$store@nodes))

    map <- .edge_to_cell_map(f$store, pe)
    # e has no node; every node has a cell
    expect_equal(nrow(map), 4L)
    expect_equal(pe@cell_ids[map$pos], nm$node_id[match(map$int_id, nm$int_id)])
    expect_equal(map$row_id, map$pos)  # an unsliced store: on-disk = view

    # after a slice, `pos` follows the view and `row_id` stays on disk
    sub <- pe[, c(5L, 3L)]
    ms <- .edge_to_cell_map(f$store, sub)
    expect_equal(sub@cell_ids[ms$pos], c("a", "c"))
    expect_equal(ms$row_id, c(5L, 3L))
})

test_that(".edge_to_cell_map carries source_id on a union", {
    f <- .edge_fixture()
    mk <- function(ids) {
        m <- Matrix::sparseMatrix(i = rep(1, length(ids)), j = seq_along(ids),
            x = 1, dims = c(1, length(ids)), dimnames = list("g1", ids))
        storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    }
    u <- unionParquetExprStore(list(mk(c("a", "b")), mk(c("c", "d"))))
    map <- .edge_to_cell_map(f$store, u)
    expect_equal(u@cell_ids[map$pos], c("a", "b", "c", "d"))
    expect_equal(map$row_id, c(1L, 2L, 1L, 2L))
    expect_equal(length(unique(map$source_id)), 2L)
})
