# Tests for streaming QC: addStreamStatistics on a parquetExprStore-backed
# giotto object. The result must be bit-for-bit identical to Giotto's
# in-memory addStatistics on the same expression values.

.tiny_mat <- function(n_genes = 12L, n_cells = 30L,
                       density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}


test_that("addStreamStatistics matches reference per-cell + per-gene stats", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")

    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # Direct reference: replicate Giotto's addStatistics formulas by hand
    # so the test does not depend on Giotto::addStatistics internals
    # (which have a sys.call-3 lookup that misbehaves under testthat).
    n_genes <- nrow(mat)
    n_cells <- ncol(mat)
    detect  <- mat > 0
    ref_cell_total <- as.numeric(Matrix::colSums(mat))
    ref_cell_nfeat <- as.integer(Matrix::colSums(detect))
    ref_feat_total <- as.numeric(Matrix::rowSums(mat))
    ref_feat_ncell <- as.integer(Matrix::rowSums(detect))

    # Parquet-backed giotto + streaming pass
    g_pq <- Giotto::createGiottoObject(expression = mat, verbose = FALSE)
    eo   <- new("exprObj", name = "raw", exprMat = pe,
                spat_unit = "cell", feat_type = "rna")
    g_pq <- GiottoClass::setExpression(g_pq, x = eo, name = "raw",
                                        verbose = FALSE)
    g_pq <- addStreamStatistics(g_pq, verbose = FALSE)

    cm <- GiottoClass::getCellMetadata(g_pq,    output = "data.table")
    fm <- GiottoClass::getFeatureMetadata(g_pq, output = "data.table")
    cm <- cm[match(colnames(mat), cm$cell_ID), ]
    fm <- fm[match(rownames(mat), fm$feat_ID), ]

    expect_equal(cm$total_expr, ref_cell_total)
    expect_equal(cm$nr_feats,   ref_cell_nfeat)
    expect_equal(cm$perc_feats, ref_cell_nfeat / n_genes * 100)

    expect_equal(fm$total_expr, ref_feat_total)
    expect_equal(fm$nr_cells,   ref_feat_ncell)
    expect_equal(fm$perc_cells, ref_feat_ncell / n_cells * 100)
    expect_equal(fm$mean_expr,  ref_feat_total / n_cells)
    expect_equal(
        fm$mean_expr_det,
        ifelse(ref_feat_ncell > 0L, ref_feat_total / ref_feat_ncell, 0)
    )
})


test_that("addStreamStatistics errors clearly on non-parquet backend", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 2)
    g   <- Giotto::createGiottoObject(expression = mat, verbose = FALSE)
    expect_error(
        addStreamStatistics(g),
        "requires a parquetExprStore"
    )
})


test_that("addStreamStatistics rejects non-giotto input", {
    expect_error(
        addStreamStatistics(list()),
        "must be a `giotto` object"
    )
})
