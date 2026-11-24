# Giotto's neighborhoodPCAParam on an expression store must equal the
# in-memory embedding.

.npca_parts <- function(G = 25L, N = 150L, seed = 3L, envir = parent.frame()) {
    withr::local_options(giotto.check_valid = FALSE, giotto.verbose = FALSE,
                         giotto.no_python_warn = TRUE, .local_envir = envir)
    set.seed(seed)
    x <- runif(N, 0, 10); y <- runif(N, 0, 10)
    lam <- cbind(matrix(1.5, N, G - 4), x, y, 5 * (x < 4), 5 * (y < 3))
    m <- t(matrix(rpois(N * G, lam), N))
    m[G - 5L, ] <- 0                                  # a feature without values
    dimnames(m) <- list(paste0("g", seq_len(G)), paste0("c", seq_len(N)))
    locs <- data.table::data.table(cell_ID = colnames(m), sdimx = x, sdimy = y)
    dir <- file.path(withr::local_tempdir(.local_envir = envir), "proj")
    g <- GiottoClass::createGiottoObject(expression = m, spatial_locs = locs,
        backend = dir)
    g <- Giotto::normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE,
        verbose = FALSE)
    store <- GiottoClass::getExpression(g, values = "normalized",
        output = "exprObj")[]
    mem <- as.matrix(storeRead(store, output = "dgcmatrix", max_rows = Inf,
        max_cols = Inf))
    xy <- as.matrix(locs[, .(sdimx, sdimy)])
    rownames(xy) <- locs$cell_ID
    list(store = store, mem = mem[, colnames(m)], locs = xy)
}

# the same kNN spatial network in memory and as an edge store
.npca_net <- function(locs, k) {
    p <- function(out) GiottoClass::kNNNetworkParam(k = k, output = out)
    list(mem = GiottoClass::createNetwork(locs, p("data.table")),
        disk = GiottoClass::createNetwork(locs, p("parquet")))
}

.npca_same <- function(got, ref) {
    expect_equal(dim(got$u), dim(ref$u))
    expect_identical(rownames(got$u), rownames(ref$u))
    r <- abs(diag(stats::cor(got$u, ref$u[rownames(got$u), ])))
    expect_gt(min(r), 1 - 1e-8)
    expect_equal(got$eigenvalues, ref$eigenvalues, tolerance = 1e-8)
}

test_that("neighborhoodPCAParam on a store equals the in-memory result", {
    skip_if_not(methods::existsMethod("reduceData",
        c("allMatrix", "neighborhoodPCAParam")))
    o <- .npca_parts()
    net <- .npca_net(o$locs, 6L)
    for (agf in c(FALSE, TRUE)) {
        p <- Giotto::neighborhoodPCAParam(lambda = 0.8, k_geom = 6L,
            npcs = 5L, compute_agf = agf)
        ref <- GiottoClass::reduceData(o$mem, p, spatial_locs = o$locs,
            network = net$mem)
        # one window, then windows small enough to need the local reorder
        for (cs in list(NULL, 12L)) {
            got <- withr::with_options(list(giottodisk.chunk_size = cs),
                GiottoClass::reduceData(o$store, p, spatial_locs = o$locs,
                    network = net$disk))
            .npca_same(got, ref)
        }
    }
})

test_that("a feature subset of the store gives the subset's embedding", {
    skip_if_not(methods::existsMethod("reduceData",
        c("allMatrix", "neighborhoodPCAParam")))
    o <- .npca_parts()
    feats <- rownames(o$mem)[c(1:4, 19:25)]
    p <- Giotto::neighborhoodPCAParam(lambda = 0.2, k_geom = 6L, npcs = 4L,
        compute_agf = TRUE)
    n6 <- .npca_net(o$locs, 6L)
    n9 <- .npca_net(o$locs, 9L)
    ref <- GiottoClass::reduceData(o$mem[feats, ], p, spatial_locs = o$locs,
        network = n6$mem, agf_network = n9$mem)
    got <- GiottoClass::reduceData(o$store[feats, ], p, spatial_locs = o$locs,
        network = n6$disk, agf_network = n9$disk)
    .npca_same(got, ref)
})

test_that("runNeighborhoodPCA on a backed object equals the in-memory one", {
    skip_if_not(methods::existsMethod("reduceData",
        c("allMatrix", "neighborhoodPCAParam")))
    withr::local_options(giotto.check_valid = FALSE, giotto.verbose = FALSE,
        giotto.no_python_warn = TRUE)
    set.seed(5)
    N <- 120L
    m <- matrix(rpois(20L * N, 3), 20L, dimnames = list(paste0("g", 1:20),
        paste0("c", seq_len(N))))
    locs <- data.table::data.table(cell_ID = colnames(m),
        sdimx = runif(N, 0, 10), sdimy = runif(N, 0, 10))
    mk <- function(...) {
        g <- GiottoClass::createGiottoObject(expression = m,
            spatial_locs = locs, ...)
        g <- Giotto::normalizeGiotto(g, scale_feats = FALSE,
            scale_cells = FALSE, verbose = FALSE)
        Giotto::runNeighborhoodPCA(g, lambda = 0.8, k_geom = 6L, npcs = 4L,
            compute_agf = TRUE)
    }
    dir <- file.path(withr::local_tempdir(), "proj")
    emb <- function(g) GiottoClass::getDimReduction(g,
        reduction_method = "neighborhood_pca", name = "neighborhood_pca",
        output = "matrix")
    a <- emb(mk())
    b <- emb(withr::with_options(list(giottodisk.chunk_size = 16L),
        mk(backend = dir)))
    expect_setequal(rownames(b), rownames(a))
    r <- abs(diag(stats::cor(a, b[rownames(a), ])))
    expect_gt(min(r), 1 - 1e-8)
})

test_that("an undirected network on a store equals the in-memory result", {
    skip_if_not(methods::existsMethod("reduceData",
        c("allMatrix", "neighborhoodPCAParam")))
    o <- .npca_parts()
    p <- Giotto::neighborhoodPCAParam(lambda = 0.8, npcs = 4L,
        compute_agf = TRUE)
    dp <- function(out) GiottoClass::delaunayNetworkParam(output = out)
    ref <- GiottoClass::reduceData(o$mem, p, spatial_locs = o$locs,
        network = GiottoClass::createNetwork(o$locs, dp("igraph")))
    disk <- GiottoClass::createNetwork(o$locs, dp("parquet"))
    expect_false(disk@directed)
    got <- withr::with_options(list(giottodisk.chunk_size = 12L),
        GiottoClass::reduceData(o$store, p, spatial_locs = o$locs,
            network = disk))
    .npca_same(got, ref)
})
