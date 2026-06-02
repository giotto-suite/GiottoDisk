# Tests for streaming normalize dispatch:
#   processData(parquetExprStore, libraryNormParam) -> appends norm_libsize_log
#       op (log = FALSE) to pe@ops
#   processData(parquetExprStore, logNormParam)     -> flips the same op's
#       log flag to TRUE (libsize+log fuse into one op)

.tiny_mat <- function(n_genes = 12L, n_cells = 30L,
                       density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}


test_that("libraryNormParam appends norm_libsize_log op to @ops", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    pe2 <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))

    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_length(pe2@ops, 1L)
    expect_equal(pe2@ops[[1]]$type, "norm_libsize_log")
    expect_equal(pe2@ops[[1]]$scalef$scalef, 1e4 / libsz)
    expect_equal(pe2@ops[[1]]$scalef$orig_row_id, seq_along(libsz))
    expect_false(pe2@ops[[1]]$log)
    expect_equal(pe2@ops[[1]]$base, 2)
})


test_that("logNormParam fuses log flag onto existing libsize op", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 2)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))
    pe  <- GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 1))

    # Still one fused op, not two appended.
    expect_length(pe@ops, 1L)
    expect_equal(pe@ops[[1]]$type, "norm_libsize_log")
    expect_true(pe@ops[[1]]$log)
    expect_equal(pe@ops[[1]]$base, 2)
    # scale_factors preserved through the fuse.
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(pe@ops[[1]]$scalef$scalef, 1e4 / libsz)
})


test_that("logNormParam with offset != 1 errors clearly", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))
    expect_error(
        GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 0.5)),
        "offset != 1"
    )
})


test_that("logNormParam without prior libsize op errors", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 6)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(
        GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 1)),
        "no.*library-size normalization op"
    )
})


test_that("normalizeGiotto end-to-end on parquet backend builds fused @ops", {
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
    expect_length(pe_norm@ops, 1L)
    expect_equal(pe_norm@ops[[1]]$type, "norm_libsize_log")
    expect_true(pe_norm@ops[[1]]$log)

    # Direct math matches stored scale_factors
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(pe_norm@ops[[1]]$scalef$scalef, 1e4 / libsz)
})


test_that("cell subset on parquetExprStore slices scalef table", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 7)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    expect_equal(nrow(pe@ops[[1]]$scalef), ncol(mat))

    pe2 <- pe[, 1:5]
    expect_equal(nrow(pe2@ops[[1]]$scalef), 5L)
    # source_id stays constant for single store
    expect_setequal(unique(pe2@ops[[1]]$scalef$source_id), pe@uid)
    expect_setequal(pe2@ops[[1]]$scalef$orig_row_id, 1:5)

    # v_norm for surviving cells must match pre-subset values
    v_full <- data.table::as.data.table(storeRead(pe,  output = "tibble"))
    v_sub  <- data.table::as.data.table(storeRead(pe2, output = "tibble"))
    common <- v_full[row_id %in% 1:5]
    data.table::setorder(common, row_id, col_id)
    data.table::setorder(v_sub, row_id, col_id)
    expect_equal(common$v_norm, v_sub$v_norm)
})


test_that("gene subset on parquetExprStore preserves scalef table", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 8)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    before <- nrow(pe@ops[[1]]$scalef)

    pe2 <- pe[1:5, ]  # 5 of the n_genes
    expect_equal(nrow(pe2@ops[[1]]$scalef), before)
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
