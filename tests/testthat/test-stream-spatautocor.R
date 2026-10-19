# Giotto's spatAutoCorParam on a disk-backed giotto object: the composition
# (binarize -> edgeStatsParam -> statistic) dispatches to the store methods
# and must equal the in-memory object's result.

.sac_objects <- function(G = 30L, N = 150L, seed = 4L, envir = parent.frame()) {
    withr::local_options(giotto.check_valid = FALSE, giotto.verbose = FALSE,
                         giotto.no_python_warn = TRUE, .local_envir = envir)
    set.seed(seed)
    x <- runif(N, 0, 10); y <- runif(N, 0, 10)
    lam <- cbind(matrix(2, N, G - 3), x, y, 6 * (x < 4))
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

test_that("the backed object runs on the store methods", {
    skip_if_not_installed("Giotto")
    o <- .sac_objects()
    expr <- GiottoClass::getExpression(o$backed, values = "normalized",
                                       output = "exprObj")[]
    net <- GiottoClass::getSpatialNetwork(o$backed, name = "Delaunay_network",
        output = "spatialNetworkObj", verbose = FALSE)[]
    expect_s4_class(expr, "parquetExprBase")
    expect_s4_class(net, "parquetEdgeStore")
})

test_that("moran and geary match the in-memory object", {
    skip_if_not_installed("Giotto")
    o <- .sac_objects()
    for (p in list(Giotto::spatAutoCorParam("moran"),
                   Giotto::spatAutoCorParam("geary"))) {
        mem <- GiottoClass::analyzeData(o$mem, p)
        got <- GiottoClass::analyzeData(o$backed, p)
        got <- got[match(mem$feats, got$feats)]
        expect_equal(got$statistic, mem$statistic, label = p$method)
        expect_equal(got$p.value, mem$p.value, label = p$method)
    }
    expect_true(all(c("g28", "g29", "g30") %in%
        head(GiottoClass::analyzeData(o$backed, Giotto::spatAutoCorParam("moran"))$feats, 3)))
})

test_that("binspect on continuous values agrees up to the binning", {
    skip_if_not_installed("Giotto")
    o <- .sac_objects()
    # normalized values are continuous: values sharing a bin become ties, so
    # a split or a rank cut can land one bin over
    for (rule in c("kmeans", "rank")) {
        p <- Giotto::spatAutoCorParam("binspect", bin_method = rule)
        mem <- GiottoClass::analyzeData(o$mem, p)
        got <- GiottoClass::analyzeData(o$backed, p)
        got <- got[match(mem$feats, got$feats)]
        expect_gte(mean(got$statistic == mem$statistic), 0.9)
        expect_equal(got$statistic, mem$statistic, tolerance = 0.1)
    }
})

test_that("binSpect fast writes binSpect.pval on a backed object", {
    skip_if_not_installed("Giotto")
    o <- .sac_objects()
    g <- Giotto::binSpect(o$backed, implementation = "fast", return_gobject = TRUE)
    fm <- GiottoClass::getFeatureMetadata(g, output = "data.table")
    ref <- Giotto::binSpect(o$mem, implementation = "fast")
    expect_equal(fm$binSpect.pval[match(ref$feats, fm$feat_ID)], ref$adj.p.value)
})

test_that("the covariate null matches the in-memory object", {
    skip_if_not_installed("Giotto")
    o <- .sac_objects()
    p <- Giotto::spatAutoCorParam("moran", null = "covariate", n_groups = 4L)
    mem <- GiottoClass::analyzeData(o$mem, p)
    got <- GiottoClass::analyzeData(o$backed, p)
    got <- got[match(mem$feats, got$feats)]
    expect_equal(got$statistic, mem$statistic)
    expect_equal(got$p.value, mem$p.value)
})
