# Tests for storeRead(output = "dgcmatrix") on parquetExprStore and
# unionParquetExprStore.
#
# * Matches raw matrix when no @ops queued
# * Matches normalized matrix when norm @ops queued (single + union)
# * Asymmetric dimension guard: at least one axis must be <= its cap
# * Bioconductor convention dimnames (rownames = feat_ids, colnames = cell_ids)
# * Sparsity preserved (uses v_norm if present, else raw value)

.dgc_mat <- function(n_genes = 20L, n_cells = 30L, density = 0.4, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", sprintf("%03d", seq_len(n_genes)))
    colnames(m) <- paste0("c", sprintf("%03d", seq_len(n_cells)))
    m
}


test_that("storeRead(dgcmatrix) on parquetExprStore matches raw matrix", {
    m <- .dgc_mat(seed = 11)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    out <- storeRead(pe, output = "dgcmatrix",
        max_rows = 100, max_cols = 100)
    expect_s4_class(out, "dgCMatrix")
    expect_equal(dim(out), dim(m))
    expect_identical(rownames(out), rownames(m))
    expect_identical(colnames(out), colnames(m))
    expect_equal(as.matrix(out), as.matrix(m))
})


test_that("storeRead(dgcmatrix) with norm @ops matches in-memory libnorm", {
    skip_if_not_installed("Giotto")
    m <- .dgc_mat(seed = 12)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    pe <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))

    out <- storeRead(pe, output = "dgcmatrix",
        max_rows = 100, max_cols = 100)
    libsz <- as.numeric(Matrix::colSums(m)); libsz[libsz == 0] <- 1
    expected <- t(t(as.matrix(m)) * (1e4 / libsz))
    expect_equal(as.matrix(out), expected, tolerance = 1e-10)
})


test_that("storeRead(dgcmatrix) with norm + log fuse applies log1p / log(base)", {
    skip_if_not_installed("Giotto")
    m <- .dgc_mat(seed = 13)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    pe <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    pe <- GiottoClass::processData(pe,
        Giotto::normParam("log", base = 2, offset = 1))

    out <- storeRead(pe, output = "dgcmatrix",
        max_rows = 100, max_cols = 100)
    libsz <- as.numeric(Matrix::colSums(m)); libsz[libsz == 0] <- 1
    expected <- log1p(t(t(as.matrix(m)) * (1e4 / libsz))) / log(2)
    expect_equal(as.matrix(out), expected, tolerance = 1e-10)
})


test_that("asymmetric guard: both axes above caps errors; one below passes", {
    big <- Matrix::rsparsematrix(200, 200, density = 0.1,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(big) <- paste0("g", 1:200)
    colnames(big) <- paste0("c", 1:200)
    pe_big <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), big)

    # 200 × 200 with default 100/100 caps — should error.
    expect_error(
        storeRead(pe_big, output = "dgcmatrix"),
        "would materialize"
    )
    # 200 × 50 (col below cap) — should pass.
    out_col <- storeRead(pe_big[, 1:50], output = "dgcmatrix")
    expect_equal(dim(out_col), c(200L, 50L))
    # 50 × 200 (row below cap) — should pass.
    out_row <- storeRead(pe_big[1:50, ], output = "dgcmatrix")
    expect_equal(dim(out_row), c(50L, 200L))
})


test_that("override via max_rows/max_cols and options", {
    big <- Matrix::rsparsematrix(150, 150, density = 0.1,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(big) <- paste0("g", 1:150); colnames(big) <- paste0("c", 1:150)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), big)

    # Default 100x100 caps -> error.
    expect_error(storeRead(pe, output = "dgcmatrix"), "would materialize")
    # Call-arg override -> passes.
    out <- storeRead(pe, output = "dgcmatrix",
        max_rows = 200, max_cols = 200)
    expect_equal(dim(out), c(150L, 150L))
    # Option override -> passes.
    GiottoUtils::gwith_options(list(
        giottodisk.dgc_max_rows = 200L,
        giottodisk.dgc_max_cols = 200L
    ), {
        out2 <- storeRead(pe, output = "dgcmatrix")
        expect_equal(dim(out2), c(150L, 150L))
    })
})


test_that("dgcmatrix output on unionParquetExprStore spans substores", {
    skip_if_not_installed("Giotto")
    mat1 <- .dgc_mat(n_genes = 15L, n_cells = 10L, seed = 21)
    mat2 <- .dgc_mat(n_genes = 15L, n_cells = 8L, seed = 22)
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat1)
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat2)
    # Make cell_ids globally unique
    pe1@cell_ids <- paste0("a_", pe1@cell_ids)
    pe2@cell_ids <- paste0("b_", pe2@cell_ids)
    u <- unionParquetExprStore(list(pe1, pe2))

    # raw union dgcmatrix
    out <- storeRead(u, output = "dgcmatrix",
        max_rows = 100, max_cols = 100)
    expect_equal(dim(out), c(15L, 18L))
    # union cell_ids in order
    expect_identical(colnames(out), c(pe1@cell_ids, pe2@cell_ids))
    # value contents = column-binding the two matrices (after raw rename)
    expected <- cbind(as.matrix(mat1), as.matrix(mat2))
    colnames(expected) <- c(pe1@cell_ids, pe2@cell_ids)
    rownames(expected) <- rownames(mat1)
    expect_equal(as.matrix(out), expected)

    # normalized union dgcmatrix
    u_n <- GiottoClass::processData(u,
        Giotto::normParam("library", scalefactor = 1e4))
    out_n <- storeRead(u_n, output = "dgcmatrix",
        max_rows = 100, max_cols = 100)
    libsz <- as.numeric(Matrix::colSums(expected))
    libsz[libsz == 0] <- 1
    expected_n <- t(t(expected) * (1e4 / libsz))
    expect_equal(as.matrix(out_n), expected_n, tolerance = 1e-10)
})
