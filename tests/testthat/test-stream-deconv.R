# Streaming spatialDWLS.
#
# The bar is parity with the in-memory method: the same cells with the same
# proportions. The solves are identical; only the cluster means are summed
# window by window, so the comparison allows floating-point slack and no more.

# Three cell types with disjoint markers over a background panel, spots in
# three clusters of different compositions. `envir` is the caller's frame, as
# in test-stream-enrich.R, so the project directory outlives this helper.
.deconv_pe_fixture <- function(N = 150L, seed = 5L, envir = parent.frame()) {
    withr::local_options(giotto.check_valid = FALSE, giotto.verbose = FALSE,
                         giotto.no_python_warn = TRUE, .local_envir = envir)
    set.seed(seed)
    G <- 120L
    genes <- paste0("g", seq_len(G))
    types <- c("tA", "tB", "tC")
    mu <- matrix(1, G, 3, dimnames = list(genes, types))
    mu[1:15, "tA"] <- 30
    mu[16:30, "tB"] <- 30
    mu[31:45, "tC"] <- 30

    clus <- rep(c("k1", "k2", "k3"), length.out = N)
    comp <- vapply(clus, function(k) {
        p <- switch(k,
            k1 = c(0.7, 0.3, 0), k2 = c(0.2, 0.5, 0.3), k3 = c(0, 0.4, 0.6))
        stats::rmultinom(1, 8, p)
    }, numeric(3))
    m <- matrix(stats::rpois(G * N, mu %*% comp), G, N,
        dimnames = list(genes, paste0("c", seq_len(N))))
    m[, 7] <- 0  # one spot with no signal at all

    locs <- data.table::data.table(
        cell_ID = colnames(m),
        sdimx = stats::runif(N, 0, 10), sdimy = stats::runif(N, 0, 10)
    )
    dir <- file.path(withr::local_tempdir(.local_envir = envir), "proj")
    make <- function(...) {
        g <- GiottoClass::createGiottoObject(
            expression = m, spatial_locs = locs, ...)
        # c7 is all-zero on purpose; normalisation warns about it
        g <- suppressWarnings(Giotto::normalizeGiotto(g, scale_feats = FALSE,
            scale_cells = FALSE, verbose = FALSE))
        GiottoClass::addCellMetadata(g,
            new_metadata = data.frame(cell_ID = colnames(m), clus = clus),
            by_column = TRUE, column_cell_ID = "cell_ID")
    }
    mem <- make()
    norm <- GiottoClass::getExpression(mem, values = "normalized",
        output = "matrix")
    # the signature: mean reverse-logged expression of the pure profiles
    sig <- 2^log2(mu / colSums(mu) * 1e4 * 0.1 + 1) - 1
    list(
        backed = make(backend = dir),
        mem = mem,
        sig = sig[1:45, ],
        clusters = stats::setNames(clus, colnames(m))
    )
}

.as_prop <- function(dt) {
    m <- as.matrix(dt[, setdiff(names(dt), "cell_ID"), with = FALSE])
    rownames(m) <- dt$cell_ID
    m[sort(rownames(m)), sort(colnames(m)), drop = FALSE]
}


test_that("streamed DWLS matches the in-memory method over several windows", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("spatialDWLS")
    f <- .deconv_pe_fixture()
    withr::local_options(giottodisk.chunk_size = 40L)

    pe <- GiottoClass::getExpression(f$backed, values = "normalized",
        output = "exprObj")[]
    expect_s4_class(pe, "parquetExprBase")
    sub <- pe[intersect(rownames(f$sig), pe@feat_ids), ]
    expect_gt(length(.pe_windows(sub,
        .pe_window_cells(sub, bytes_per_nz = 48, k = 0L))), 1L)
    mat <- GiottoClass::getExpression(f$mem, values = "normalized",
        output = "matrix")
    p <- Giotto::deconvParam("DWLS", n_cell = 5)

    streamed <- analyzeData(pe, p, sign_matrix = f$sig,
        cell_clusters = f$clusters)
    in_mem <- analyzeData(mat, p, sign_matrix = f$sig,
        cell_clusters = f$clusters)

    expect_identical(sort(streamed$cell_ID), sort(in_mem$cell_ID))
    expect_equal(.as_prop(streamed), .as_prop(in_mem), tolerance = 1e-10)

    # the no-signal spot is unsolvable on both backends
    expect_true(all(is.na(.as_prop(streamed)["c7", ])))
})


test_that("runDWLSDeconv() runs end to end on a backed object", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("spatialDWLS")
    f <- .deconv_pe_fixture()
    withr::local_options(giottodisk.chunk_size = 40L)

    res_b <- Giotto::runDWLSDeconv(f$backed, sign_matrix = f$sig,
        cluster_column = "clus", n_cell = 5, return_gobject = FALSE,
        na_as_zero = TRUE)
    res_m <- Giotto::runDWLSDeconv(f$mem, sign_matrix = f$sig,
        cluster_column = "clus", n_cell = 5, return_gobject = FALSE,
        na_as_zero = TRUE)
    expect_s4_class(res_b, "spatEnrObj")
    expect_equal(.as_prop(res_b[]), .as_prop(res_m[]), tolerance = 1e-10)
    expect_false(anyNA(res_b[]))
})


test_that("unnamed clusters follow the store's cell order", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("spatialDWLS")
    f <- .deconv_pe_fixture()
    pe <- GiottoClass::getExpression(f$backed, values = "normalized",
        output = "exprObj")[]
    p <- Giotto::deconvParam("DWLS", n_cell = 5)

    named <- analyzeData(pe, p, sign_matrix = f$sig,
        cell_clusters = f$clusters)
    unnamed <- analyzeData(pe, p, sign_matrix = f$sig,
        cell_clusters = unname(f$clusters[pe@cell_ids]))
    expect_equal(unnamed, named)

    expect_error(
        analyzeData(pe, p, sign_matrix = f$sig,
            cell_clusters = unname(f$clusters)[-1]),
        "one value per cell"
    )
})
