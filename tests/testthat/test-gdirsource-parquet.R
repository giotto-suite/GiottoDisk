# Tests for gDirSource integration of parquetExpr (Phase 3 Layer 3)

.tiny_mat <- function(n_genes = 12L, n_cells = 30L, density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}


test_that("storeCreate('parquetExpr') builds a parquetExprStore", {
    s <- storeCreate(type = "parquetExpr")
    expect_s4_class(s, "parquetExprStore")
})


test_that("storeCreate accepts the long-form alias 'parquetexprstore'", {
    s <- storeCreate(type = "parquetexprstore")
    expect_s4_class(s, "parquetExprStore")
})


test_that("sourceWrite(gDirSource, dgCMatrix, store_type = 'parquetExpr') works", {
    proj_dir <- file.path(tempdir(),
                            paste0("gd_pq_", as.integer(Sys.time())))
    src <- gDirSource(proj_dir)

    mat <- .tiny_mat()
    written <- sourceWrite(src, mat, store_type = "parquetExpr")

    expect_s4_class(written, "parquetExprStore")
    expect_equal(nrow(written), nrow(mat))
    expect_equal(ncol(written), ncol(mat))

    # Re-open via storeRead
    ds <- storeRead(written)
    expect_true(inherits(ds, "Dataset"))

    df <- dplyr::collect(ds)
    expect_equal(nrow(df), Matrix::nnzero(mat))

    # Element-wise round-trip
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(written), ncol(written))
    )
    expect_equal(unname(as.matrix(mat)), unname(as.matrix(rt)))

    # Manifest file exists
    expect_true(file.exists(file.path(proj_dir, "giottodir.json")))

    unlink(proj_dir, recursive = TRUE)
})


test_that("parquetExpr artifact path has the expected layout", {
    proj_dir <- file.path(tempdir(),
                            paste0("gd_layout_", as.integer(Sys.time())))
    src <- gDirSource(proj_dir)
    mat <- .tiny_mat(seed = 11)
    written <- sourceWrite(src, mat, store_type = "parquetExpr")

    # Path is under proj_dir/artifacts/<uid>/  (use normalized real paths
    # because macOS resolves /var -> /private/var)
    artifacts_real <- normalizePath(file.path(proj_dir, "artifacts"),
                                       mustWork = FALSE)
    written_real   <- normalizePath(written@path, mustWork = FALSE)
    expect_true(startsWith(written_real, artifacts_real))
    expect_true(grepl(written@uid, written@path, fixed = TRUE))

    unlink(proj_dir, recursive = TRUE)
})
