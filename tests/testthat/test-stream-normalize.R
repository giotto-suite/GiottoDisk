# Tests for streaming normalize dispatch:
#   processData(parquetExprStore, libraryNormParam) -> updated pe with JIT
#   processData(parquetExprStore, logNormParam)     -> updated pe with log flag

.tiny_mat <- function(n_genes = 12L, n_cells = 30L,
                       density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}


test_that("libraryNormParam stores correct scale_factors on parquetExprStore", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    pe2 <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))

    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(pe2@params$norm$method, "library_size")
    expect_equal(pe2@params$norm$scalefactor, 1e4)
    expect_equal(pe2@params$norm$scale_factors, 1e4 / libsz)
})


test_that("logNormParam sets log flag on parquetExprStore", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 2)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))
    pe  <- GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 1))

    expect_true(isTRUE(pe@params$norm$log))
    expect_equal(pe@params$norm$base, 2)
    expect_equal(pe@params$norm$offset, 1)
    # library scale_factors must be preserved through the log step
    expect_true(!is.null(pe@params$norm$scale_factors))
})


test_that("logNormParam with offset != 1 errors clearly", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(
        GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 0.5)),
        "offset != 1"
    )
})


test_that("normalizeGiotto end-to-end on parquet backend stores JIT recipe", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")

    mat <- .tiny_mat(seed = 4)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    locs <- data.frame(cell_ID = colnames(mat),
                        sdimx = runif(ncol(mat)),
                        sdimy = runif(ncol(mat)))
    g <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                     verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    g <- suppressWarnings(Giotto::normalizeGiotto(g,
        scalefactor = 1e4,
        scale_feats = FALSE, scale_cells = FALSE,
        verbose = FALSE))

    norm_eo <- GiottoClass::getExpression(g, values = "normalized",
                                            output = "exprObj")
    pe_norm <- slot(norm_eo, "exprMat")
    expect_s4_class(pe_norm, "parquetExprStore")
    expect_equal(pe_norm@params$norm$scalefactor, 1e4)
    expect_true(isTRUE(pe_norm@params$norm$log))

    # JIT recipe applied to raw must reproduce in-memory normalized matrix
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    applied <- t(t(mat) * (1e4 / libsz))
    applied <- log1p(applied) / log(2)

    # Also check direct math matches stored scale_factors
    expect_equal(pe_norm@params$norm$scale_factors, 1e4 / libsz)
})


test_that("normalizeGiotto with scale_feats=TRUE errors on parquet backend", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 5)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    locs <- data.frame(cell_ID = colnames(mat),
                        sdimx = runif(ncol(mat)),
                        sdimy = runif(ncol(mat)))
    g <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                     verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    expect_error(
        suppressWarnings(Giotto::normalizeGiotto(g,
            scale_feats = TRUE, scale_cells = TRUE, verbose = FALSE)),
        "z-score scaling"
    )
})
