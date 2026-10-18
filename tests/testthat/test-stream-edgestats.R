# analyzeData(parquetExprBase, edgeStatsParam) against the in-memory method in
# GiottoClass on the same matrix and edge list: the whole contract must agree.

.es_fixture <- function(n = 40L, G = 6L, seed = 8L, prefix = "c") {
    set.seed(seed)
    ids <- paste0(prefix, seq_len(n))
    m <- Matrix::rsparsematrix(G, n, density = 0.5,
        rand.x = function(k) as.double(rpois(k, 4L) + 1L))
    dimnames(m) <- list(paste0("g", seq_len(G)), ids)
    from <- rep(seq_len(n), each = 3L)
    to <- ((from - 1L + rep(c(1L, 3L, 7L), n)) %% n) + 1L
    e <- unique(data.table::data.table(a = pmin(from, to), b = pmax(from, to)))
    net <- data.frame(from = ids[e$a], to = ids[e$b],
                      weight = round(runif(nrow(e), 0.5, 2), 3))
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
        data.table::as.data.table(net), type = "sNN", directed = FALSE)
    list(m = m, net = net, store = s)
}
.es_pe <- function(m) {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}
.es_expect_contract <- function(got, ref) {
    expect_equal(got$feats, ref$feats)
    for (nm in setdiff(names(ref), "feats")) {
        expect_equal(got[[nm]], unname(ref[[nm]]), label = nm)
    }
    expect_equal(attr(got, "graph"), attr(ref, "graph"))
}

test_that("the store method reproduces the in-memory contract", {
    f <- .es_fixture()
    pe <- .es_pe(f$m)
    for (wt in c(FALSE, TRUE)) {
        p <- GiottoClass::edgeStatsParam(weighted = wt)
        .es_expect_contract(
            GiottoClass::analyzeData(pe, p, network = f$store),
            GiottoClass::analyzeData(f$m, p, network = f$net))
    }
    expect_length(pe@ops, 0L)
})

test_that("only edges with both ends in the view count", {
    f <- .es_fixture()
    pe <- .es_pe(f$m)
    cells <- c(2:25, 31L)
    p <- GiottoClass::edgeStatsParam()
    .es_expect_contract(
        GiottoClass::analyzeData(pe[c(1L, 4L, 5L), cells], p, network = f$store),
        GiottoClass::analyzeData(f$m[c(1L, 4L, 5L), cells], p, network = f$net))
})

test_that("a union and a normalized chain give the in-memory contract", {
    f <- .es_fixture()
    u <- unionParquetExprStore(list(.es_pe(f$m[, 1:15]), .es_pe(f$m[, 16:40])))
    p <- GiottoClass::edgeStatsParam(weighted = TRUE)
    .es_expect_contract(GiottoClass::analyzeData(u, p, network = f$store),
                        GiottoClass::analyzeData(f$m, p, network = f$net))

    pe <- .es_pe(f$m) |>
        GiottoClass::processData(Giotto::normParam("library", scalefactor = 100)) |>
        GiottoClass::processData(Giotto::normParam("log", base = 2, offset = 1))
    x <- f$m %*% Matrix::Diagonal(x = 100 / Matrix::colSums(f$m))
    x@x <- log2(x@x + 1)
    dimnames(x) <- dimnames(f$m)
    .es_expect_contract(
        GiottoClass::analyzeData(pe, GiottoClass::edgeStatsParam(), network = f$store),
        GiottoClass::analyzeData(x, GiottoClass::edgeStatsParam(), network = f$net))
})

test_that("the multiply path for many degree classes agrees with the grouped one", {
    f <- .es_fixture()
    pe <- .es_pe(f$m)
    deg <- runif(ncol(f$m), 0, 5)
    expect_equal(.pe_degree_sums(pe, deg, max_classes = 0L),
                 .pe_degree_sums(pe, deg))
})

test_that("the grouped attribute reproduces the in-memory one", {
    f <- .es_fixture()
    ids <- colnames(f$m)
    grp <- stats::setNames(rep(c("lo", "mid", "hi"), length.out = length(ids)), ids)
    u <- unionParquetExprStore(list(.es_pe(f$m[, 1:15]), .es_pe(f$m[, 16:40])))
    for (wt in c(FALSE, TRUE)) {
        p <- GiottoClass::edgeStatsParam(weighted = wt)
        ref <- attr(GiottoClass::analyzeData(f$m, p, network = f$net,
            groups = grp), "grouped")
        for (x in list(.es_pe(f$m), u)) {
            got <- attr(GiottoClass::analyzeData(x, p, network = f$store,
                groups = rev(grp)), "grouped")
            expect_equal(got, ref)
        }
    }
    sub_cells <- c(3:30)
    p <- GiottoClass::edgeStatsParam()
    expect_equal(
        attr(GiottoClass::analyzeData(.es_pe(f$m)[, sub_cells], p,
            network = f$store, groups = grp), "grouped"),
        attr(GiottoClass::analyzeData(f$m[, sub_cells], p, network = f$net,
            groups = grp), "grouped"))
})

test_that("a directed store's mutual pairs count once, at their mean weight", {
    f <- .es_fixture()
    both <- data.table::rbindlist(list(f$net,
        data.frame(from = f$net$to, to = f$net$from, weight = f$net$weight * 3)))
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), both,
        type = "kNN", directed = TRUE)
    mean_w <- f$net
    mean_w$weight <- f$net$weight * 2
    pe <- .es_pe(f$m)
    for (wt in c(FALSE, TRUE)) {
        p <- GiottoClass::edgeStatsParam(weighted = wt)
        .es_expect_contract(GiottoClass::analyzeData(pe, p, network = s),
            GiottoClass::analyzeData(f$m, p, network = if (wt) mean_w else f$net))
    }
})
