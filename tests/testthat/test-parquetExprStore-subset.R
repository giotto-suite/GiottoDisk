# Tests for parquetExprStore `[` subset method (Step 2.3-bis)

.tiny_mat <- function(n_genes = 12L, n_cells = 30L, density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}


test_that("`[` with integer indices narrows feat_ids and cell_ids", {
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    sub <- pe[c(2L, 5L, 7L), c(1L, 3L, 10L, 20L)]
    expect_s4_class(sub, "parquetExprStore")
    expect_equal(sub@n_genes, 3)
    expect_equal(sub@n_cells, 4)
    expect_equal(sub@feat_ids, paste0("g", c(2, 5, 7)))
    expect_equal(sub@cell_ids, paste0("c", c(1, 3, 10, 20)))
    expect_equal(sub@gene_idx, c(2L, 5L, 7L))
    expect_equal(sub@cell_idx, c(1L, 3L, 10L, 20L))
})


test_that("`[` with character IDs subsets correctly", {
    mat <- .tiny_mat(seed = 2)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    sub <- pe[c("g3", "g8"), c("c5", "c10", "c25")]
    expect_equal(sub@n_genes, 2)
    expect_equal(sub@n_cells, 3)
    expect_equal(sub@feat_ids, c("g3", "g8"))
    expect_equal(sub@cell_ids, c("c5", "c10", "c25"))
    expect_equal(sub@gene_idx, c(3L, 8L))
    expect_equal(sub@cell_idx, c(5L, 10L, 25L))
})


test_that("`[` with logical indices works", {
    mat <- .tiny_mat(n_genes = 5L, n_cells = 10L, seed = 3)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    keep_g <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
    keep_c <- rep(c(TRUE, FALSE), 5)

    sub <- pe[keep_g, keep_c]
    expect_equal(sub@n_genes, 3)
    expect_equal(sub@n_cells, 5)
    expect_equal(sub@gene_idx, c(1L, 3L, 5L))
    expect_equal(sub@cell_idx, c(1L, 3L, 5L, 7L, 9L))
})


test_that("`[` chains correctly (subset of subset)", {
    mat <- .tiny_mat(n_genes = 10L, n_cells = 20L, seed = 4)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    s1 <- pe[3L:8L, 5L:15L]               # 6 genes x 11 cells
    s2 <- s1[c(1L, 4L, 6L), c(2L, 5L, 9L)] # 3 x 3
    # Original-parquet positions
    expect_equal(s2@gene_idx, c(3L, 6L, 8L))
    expect_equal(s2@cell_idx, c(6L, 9L, 13L))
    expect_equal(s2@feat_ids, rownames(mat)[c(3, 6, 8)])
    expect_equal(s2@cell_ids, colnames(mat)[c(6, 9, 13)])
})


test_that("dim / dimnames respect the subset", {
    mat <- .tiny_mat(seed = 5)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    sub <- pe[1L:4L, 1L:7L]
    expect_equal(dim(sub), c(4, 7))
    expect_equal(nrow(sub), 4)
    expect_equal(ncol(sub), 7)
    expect_equal(rownames(sub), paste0("g", 1:4))
    expect_equal(colnames(sub), paste0("c", 1:7))
})


test_that("storeRead on a subset filters via Arrow (lazy, no rewrite)", {
    mat <- .tiny_mat(seed = 6)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    sub <- pe[c(2L, 4L, 9L), c(3L, 7L, 12L)]

    df <- dplyr::collect(storeRead(sub))
    # All returned row_ids should be in the subset's cell_idx
    expect_true(all(df$row_id %in% sub@cell_idx))
    expect_true(all(df$col_id %in% sub@gene_idx))

    # Element-wise check: reconstruct the subset matrix and compare
    rt <- Matrix::sparseMatrix(
        i = match(df$col_id, sub@gene_idx),
        j = match(df$row_id, sub@cell_idx),
        x = df$value,
        dims = c(nrow(sub), ncol(sub))
    )
    rownames(rt) <- sub@feat_ids
    colnames(rt) <- sub@cell_ids
    ref <- mat[c(2, 4, 9), c(3, 7, 12)]
    expect_equal(unname(as.matrix(rt)), unname(as.matrix(ref)))
})


test_that("`[` errors on invalid character IDs", {
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(pe[c("g1", "missing_gene"), ], "unknown IDs")
    expect_error(pe[, c("c1", "missing_cell")], "unknown IDs")
})


test_that("`[` errors on logical of wrong length", {
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(pe[c(TRUE, FALSE), ], "logical row")
})


test_that("filterGiotto applies on parquet backend (Step 2.3-bis end-to-end)", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(n_genes = 30L, n_cells = 80L, density = 0.4, seed = 99)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    locs <- data.frame(cell_ID = colnames(mat),
                        sdimx = runif(ncol(mat)),
                        sdimy = runif(ncol(mat)))
    g <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                     verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    # Reference: filterGiotto on the in-memory dgCMatrix path
    g_mem <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                         verbose = FALSE)
    g_mem_f <- suppressWarnings(Giotto::filterGiotto(g_mem,
        expression_threshold   = 1,
        feat_det_in_min_cells  = 5,
        min_det_feats_per_cell = 5,
        verbose = FALSE))

    # Streaming path
    g_pq_f <- suppressWarnings(Giotto::filterGiotto(g,
        expression_threshold   = 1,
        feat_det_in_min_cells  = 5,
        min_det_feats_per_cell = 5,
        verbose = FALSE))

    # Both should keep the same cell_IDs and feat_IDs
    expect_setequal(
        colnames(GiottoClass::getExpression(g_pq_f)),
        colnames(GiottoClass::getExpression(g_mem_f))
    )
    expect_setequal(
        rownames(GiottoClass::getExpression(g_pq_f)),
        rownames(GiottoClass::getExpression(g_mem_f))
    )

    # Backend after filter is still parquetExprStore (lazy subset)
    em <- GiottoClass::getExpression(g_pq_f)
    expect_s4_class(slot(em, "exprMat"), "parquetExprStore")
})
