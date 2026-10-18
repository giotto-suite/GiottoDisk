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

    d <- .edge_degree(storeRead(f$store, output = "arrow"))
    expect_equal(d$deg, unname(igraph::degree(g)[name_of(d)]))
    dw <- .edge_degree(storeRead(f$store, output = "arrow"), weighted = TRUE)
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

# .pe_edge_products: ground truth is rowSums(X[, from] * X[, to]) over the
# undirected edge list, on the in-memory matrix with the same chain applied.

.ep_fixture <- function(n = 30L, G = 8L, seed = 5L, prefix = "c") {
    set.seed(seed)
    ids <- paste0(prefix, seq_len(n))
    m <- Matrix::rsparsematrix(G, n, density = 0.5,
        rand.x = function(k) as.double(rpois(k, 4L) + 1L))
    dimnames(m) <- list(paste0("g", seq_len(G)), ids)
    from <- rep(seq_len(n), each = 3L)
    to <- ((from - 1L + rep(c(1L, 2L, 5L), n)) %% n) + 1L
    e <- unique(data.table::data.table(a = pmin(from, to), b = pmax(from, to)))
    e[, weight := round(runif(.N, 0.5, 2), 2)]
    dt <- data.table::data.table(from = ids[e$a], to = ids[e$b],
                                 weight = e$weight)
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), dt,
        type = "sNN", directed = FALSE)
    list(m = m, e = e, dt = dt, store = s)
}
.ep_ref <- function(m, from, to, w = rep(1, length(from))) {
    unname(Matrix::rowSums(m[, from, drop = FALSE] *
        Matrix::t(Matrix::t(m[, to, drop = FALSE]) * w)))
}
.ep_pe <- function(m) {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}

test_that(".pe_edge_products matches the in-memory edge sum", {
    f <- .ep_fixture()
    pe <- .ep_pe(f$m)
    expect_equal(.pe_edge_products(pe, f$store), .ep_ref(f$m, f$e$a, f$e$b))
    expect_equal(.pe_edge_products(pe, f$store, weighted = TRUE),
                 .ep_ref(f$m, f$e$a, f$e$b, f$e$weight))
    expect_length(pe@ops, 0L)
})

test_that(".pe_edge_products sees the values the lazy chain produces", {
    f <- .ep_fixture()
    pe <- .ep_pe(f$m) |>
        GiottoClass::processData(Giotto::normParam("library", scalefactor = 100)) |>
        GiottoClass::processData(Giotto::normParam("log", base = 2, offset = 1))
    lib <- Matrix::colSums(f$m)
    x <- f$m %*% Matrix::Diagonal(x = 100 / lib)
    x@x <- log2(x@x + 1)
    dimnames(x) <- dimnames(f$m)
    expect_equal(.pe_edge_products(pe, f$store), .ep_ref(x, f$e$a, f$e$b))
})

test_that(".pe_edge_products follows a sliced view", {
    f <- .ep_fixture()
    pe <- .ep_pe(f$m)
    cells <- c(3:20, 25L)
    genes <- c(2L, 5L, 7L)
    sub <- pe[genes, cells]
    keep <- f$e$a %in% cells & f$e$b %in% cells
    expect_equal(.pe_edge_products(sub, f$store),
                 .ep_ref(f$m[genes, ], f$e$a[keep], f$e$b[keep]))
})

test_that(".pe_edge_products matches features across a union", {
    f <- .ep_fixture()
    m <- f$m
    u <- unionParquetExprStore(list(.ep_pe(m[, 1:12]), .ep_pe(m[, 13:30])))
    # edges within and across the two substores
    expect_true(any(f$e$a <= 12 & f$e$b > 12))
    expect_equal(.pe_edge_products(u, f$store), .ep_ref(m, f$e$a, f$e$b))
    # a gene slice narrows each substore's on-disk index
    expect_equal(.pe_edge_products(u[c(6L, 2L), ], f$store),
                 .ep_ref(m[c(6L, 2L), ], f$e$a, f$e$b))
})

test_that(".pe_edge_products refuses a chain with post ops", {
    f <- .ep_fixture()
    pe <- .pe_push_op(.ep_pe(f$m), list(type = "log", base = 2),
                      phase = "post")
    expect_error(.pe_edge_products(pe, f$store), "post")
})

test_that(".pe_edge_products is the same in many windows as in one", {
    f <- .ep_fixture()
    pe <- .ep_pe(f$m)
    u <- unionParquetExprStore(list(.ep_pe(f$m[, 1:12]), .ep_pe(f$m[, 13:30])))
    ref <- .ep_ref(f$m, f$e$a, f$e$b, f$e$weight)
    withr::with_options(list(giottodisk.chunk_size = 4L), {
        expect_equal(.pe_edge_products(pe, f$store, weighted = TRUE), ref)
        expect_equal(.pe_edge_products(u, f$store, weighted = TRUE), ref)
        cells <- c(3:20, 25L)
        keep <- f$e$a %in% cells & f$e$b %in% cells
        expect_equal(.pe_edge_products(pe[c(2L, 5L), cells], f$store),
                     .ep_ref(f$m[c(2L, 5L), ], f$e$a[keep], f$e$b[keep]))
    })
})

test_that("a store in non-local cell order is rewritten, with the same result", {
    f <- .ep_fixture(n = 60L)
    set.seed(1)
    shuffled <- f$m[, sample(ncol(f$m))]
    pe <- .ep_pe(shuffled)
    ref <- .ep_ref(f$m, f$e$a, f$e$b, f$e$weight)
    withr::with_options(list(giottodisk.chunk_size = 8L,
                             giottodisk.edge_rewrite_min_windows = 2L), {
        pairs <- .edge_view_pairs(f$store, .edge_to_cell_map(f$store, pe), TRUE)
        expect_false(.pairs_local(pairs, 8L))
        dump <- withr::local_tempdir()
        withr::local_options(giottodisk.artifact_dump = dump)
        expect_equal(.pe_edge_products(pe, f$store, weighted = TRUE), ref)
        # the rewritten copy and its spill are removed afterwards
        expect_length(list.files(dump, recursive = TRUE), 0L)
        u <- unionParquetExprStore(list(.ep_pe(shuffled[, 1:25]),
                                        .ep_pe(shuffled[, 26:60])))
        expect_equal(.pe_edge_products(u, f$store, weighted = TRUE), ref)
    })
})

test_that(".pe_rewrite_cell_order writes the view's values in the new order", {
    f <- .ep_fixture(n = 40L)
    m <- f$m
    pe <- .ep_pe(m) |>
        GiottoClass::processData(Giotto::normParam("library", scalefactor = 100))
    ord <- c(seq(2L, 40L, by = 2L), seq(1L, 39L, by = 2L))
    withr::with_options(list(giottodisk.chunk_size = 5L), {
        for (src in list(pe, pe[c(1L, 4L, 6L), ],
                         unionParquetExprStore(list(.ep_pe(m[, 1:17]), .ep_pe(m[, 18:40]))))) {
            out <- .pe_rewrite_cell_order(src, ord, tempfile(fileext = ".parquet"))
            got <- as.matrix(storeRead(out, output = "dgcmatrix",
                max_rows = Inf, max_cols = Inf))
            want <- as.matrix(storeRead(src, output = "dgcmatrix",
                max_rows = Inf, max_cols = Inf))[, ord]
            expect_equal(got, want, ignore_attr = TRUE)
            expect_equal(out@cell_ids, src@cell_ids[ord])
            expect_length(out@ops, 0L)
            rid <- dplyr::collect(storeRead(out, output = "query"))$row_id
            expect_false(is.unsorted(rid))
        }
    })
})

test_that("with few windows a non-local store uses breadth-first windows", {
    f <- .ep_fixture(n = 60L)
    set.seed(2)
    pe <- .ep_pe(f$m[, sample(ncol(f$m))])
    withr::with_options(list(giottodisk.chunk_size = 8L,
                             giottodisk.edge_rewrite_min_windows = 100L), {
        dump <- withr::local_tempdir()
        withr::local_options(giottodisk.artifact_dump = dump)
        expect_equal(.pe_edge_products(pe, f$store, weighted = TRUE),
                     .ep_ref(f$m, f$e$a, f$e$b, f$e$weight))
        expect_length(list.files(dump, recursive = TRUE), 0L)
    })
})

test_that("the GiottoKernels path matches the R path, in feature blocks", {
    skip_if_not(.gd_kernel("edge_products_stream"))
    f <- .ep_fixture(n = 60L)
    set.seed(3)
    pe <- .ep_pe(f$m[, sample(ncol(f$m))]) |>
        GiottoClass::processData(Giotto::normParam("library", scalefactor = 100))
    r_path <- function(x, ...) withr::with_options(
        list(giottodisk.use_kernels = FALSE), .pe_edge_products(x, f$store, ...))
    for (blk in list(NULL, 20, 1)) {
        withr::with_options(list(giottodisk.edge_kernel_block_nnz = blk), {
            expect_equal(.pe_edge_products(pe, f$store, weighted = TRUE),
                         r_path(pe, weighted = TRUE))
            expect_equal(.pe_edge_products(pe[c(7L, 2L, 5L), ], f$store),
                         r_path(pe[c(7L, 2L, 5L), ]))
        })
    }
})
