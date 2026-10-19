# spatialSmoothParam / featCorParam on expression stores, and Giotto's
# spatialCoexpression() on a backed object: each must equal the in-memory
# result.

.ss_objects <- function(G = 25L, N = 120L, seed = 7L, envir = parent.frame()) {
    withr::local_options(giotto.check_valid = FALSE, giotto.verbose = FALSE,
                         giotto.no_python_warn = TRUE, .local_envir = envir)
    set.seed(seed)
    x <- runif(N, 0, 10); y <- runif(N, 0, 10)
    lam <- cbind(matrix(1.5, N, G - 4), x, y, 5 * (x < 4), 5 * (y < 3))
    m <- t(matrix(rpois(N * G, lam), N))
    dimnames(m) <- list(paste0("g", seq_len(G)), paste0("c", seq_len(N)))
    locs <- data.table::data.table(cell_ID = colnames(m), sdimx = x, sdimy = y)
    dir <- file.path(withr::local_tempdir(.local_envir = envir), "proj")
    mk <- function(...) {
        g <- GiottoClass::createGiottoObject(expression = m,
            spatial_locs = locs, ...)
        g <- Giotto::normalizeGiotto(g, scale_feats = FALSE,
            scale_cells = FALSE, verbose = FALSE)
        suppressMessages(GiottoClass::createSpatialNetwork(g,
            method = "Delaunay", verbose = FALSE))
    }
    list(mem = mk(), backed = mk(backend = dir))
}

.ss_parts <- function(g) {
    list(
        x = GiottoClass::getExpression(g, values = "normalized",
            output = "exprObj")[],
        net = GiottoClass::getSpatialNetwork(g, name = "Delaunay_network",
            output = "spatialNetworkObj", verbose = FALSE)[]
    )
}

.ss_mat <- function(store) {
    m <- storeRead(store, output = "dgcmatrix", max_rows = Inf,
        max_cols = Inf)
    as.matrix(m)
}

test_that("smoothing a store equals the in-memory result", {
    skip_if_not_installed("Giotto")
    o <- .ss_objects()
    mem <- .ss_parts(o$mem)
    dsk <- .ss_parts(o$backed)
    feats <- rownames(mem$x)[c(2, 5, 21:25)]
    for (p in list(GiottoClass::spatialSmoothParam(smoothing = 0.6),
                   GiottoClass::spatialSmoothParam(),
                   GiottoClass::spatialSmoothParam(smoothing = 0.5, hops = 2L))) {
        ref <- as.matrix(GiottoClass::processData(mem$x[feats, ], p,
            network = mem$net))
        # one window, then windows small enough to need the network-order
        # rewrite of the (randomly ordered) cells
        for (cs in list(NULL, 12L)) {
            got <- withr::with_options(list(giottodisk.chunk_size = cs),
                GiottoClass::processData(dsk$x[feats, ], p,
                    network = dsk$net))
            expect_s4_class(got, "parquetExprStore")
            gm <- .ss_mat(got)
            expect_identical(colnames(gm), colnames(ref))
            expect_equal(gm[rownames(ref), ], ref, tolerance = 1e-12)
        }
    }
})

test_that("feature correlation of a store equals cor()", {
    skip_if_not_installed("Giotto")
    o <- .ss_objects()
    mem <- .ss_parts(o$mem)
    dsk <- .ss_parts(o$backed)
    got <- GiottoClass::analyzeData(dsk$x, GiottoClass::featCorParam())
    ref <- GiottoClass::analyzeData(mem$x, GiottoClass::featCorParam())
    expect_equal(got[rownames(ref), colnames(ref)], ref, tolerance = 1e-10)
    expect_error(GiottoClass::analyzeData(dsk$x,
        GiottoClass::featCorParam("spearman")), "pearson")
})

test_that("spatialCoexpression on a backed object equals in memory", {
    skip_if_not_installed("Giotto")
    skip_if_not("spatialCoexpression" %in% getNamespaceExports("Giotto"))
    o <- .ss_objects()
    feats <- paste0("g", c(1:6, 21:25))
    args <- list(feats = feats, n_modules = 3L, return_gobject = FALSE)
    mem <- do.call(Giotto::spatialCoexpression, c(list(o$mem), args))
    dsk <- do.call(Giotto::spatialCoexpression, c(list(o$backed), args))
    expect_equal(dsk$cor[feats, feats], mem$cor[feats, feats],
        tolerance = 1e-10)
    expect_identical(dsk$modules[feats], mem$modules[feats])

    g <- Giotto::spatialCoexpression(o$backed, feats = feats, n_modules = 3L)
    sc <- GiottoClass::getSpatialEnrichment(g, name = "spatial_module",
        output = "data.table")
    gm <- Giotto::spatialCoexpression(o$mem, feats = feats, n_modules = 3L)
    sm <- GiottoClass::getSpatialEnrichment(gm, name = "spatial_module",
        output = "data.table")
    expect_equal(sc[match(sm$cell_ID, cell_ID), !"cell_ID"],
        sm[, !"cell_ID"], tolerance = 1e-10)
})

test_that("correlation of smoothed values without storing them equals memory", {
    skip_if_not_installed("Giotto")
    o <- .ss_objects()
    mem <- .ss_parts(o$mem)
    dsk <- .ss_parts(o$backed)
    feats <- rownames(mem$x)[c(1:4, 21:25)]
    p <- GiottoClass::featCorParam(
        smooth = GiottoClass::spatialSmoothParam(smoothing = 0.7))
    ref <- GiottoClass::analyzeData(mem$x[feats, ], p, network = mem$net)
    for (cs in list(NULL, 12L)) {
        got <- withr::with_options(list(giottodisk.chunk_size = cs),
            GiottoClass::analyzeData(dsk$x[feats, ], p, network = dsk$net))
        expect_equal(got[feats, feats], ref[feats, feats], tolerance = 1e-10)
    }
})

test_that("crossprod of a weight matrix with a store equals memory", {
    skip_if_not_installed("Giotto")
    o <- .ss_objects()
    mem <- .ss_parts(o$mem)
    dsk <- .ss_parts(o$backed)
    W <- matrix(runif(nrow(mem$x) * 3), ncol = 3,
        dimnames = list(NULL, c("a", "b", "c")))
    ref <- as.matrix(Matrix::crossprod(W, mem$x))
    for (cs in list(NULL, 12L)) {
        got <- withr::with_options(list(giottodisk.chunk_size = cs),
            Matrix::crossprod(W, dsk$x))
        expect_equal(got[, colnames(ref)], ref, tolerance = 1e-12)
    }
})
