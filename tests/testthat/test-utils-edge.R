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
