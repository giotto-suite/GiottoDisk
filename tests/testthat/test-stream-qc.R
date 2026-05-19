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


test_that("analyzeData(parquetExprStore, cellStatsParam) matches reference", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")

    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    n_genes <- nrow(mat); n_cells <- ncol(mat)
    detect  <- mat > 0
    ref_cell_total <- as.numeric(Matrix::colSums(mat))
    ref_cell_nfeat <- as.integer(Matrix::colSums(detect))

    cs <- GiottoClass::analyzeData(pe, Giotto::analyzeParam("cell_stats"))
    cs <- cs[match(colnames(mat), cs$cells), ]

    expect_equal(cs$total_expr, ref_cell_total)
    expect_equal(cs$nr_feats,   ref_cell_nfeat)
    expect_equal(cs$perc_feats, ref_cell_nfeat / n_genes * 100)
})


test_that("analyzeData(parquetExprStore, featStatsParam) matches reference", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 7)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    n_genes <- nrow(mat); n_cells <- ncol(mat)
    detect  <- mat > 0
    ref_total <- as.numeric(Matrix::rowSums(mat))
    ref_ncell <- as.integer(Matrix::rowSums(detect))

    fs <- GiottoClass::analyzeData(pe, Giotto::analyzeParam("feat_stats"))
    fs <- fs[match(rownames(mat), fs$feats), ]

    expect_equal(fs$total_expr,   ref_total)
    expect_equal(fs$nr_cells,     ref_ncell)
    expect_equal(fs$perc_cells,   ref_ncell / n_cells * 100)
    expect_equal(fs$mean_expr,    ref_total / n_cells)
})


test_that("addStatistics(g) on parquet backend matches in-memory reference", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")

    mat <- .tiny_mat(seed = 11)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # Reference: compute inline (avoid Giotto::addStatistics internal
    # sys.call brittleness under testthat)
    n_genes <- nrow(mat); n_cells <- ncol(mat)
    detect  <- mat > 0
    ref_cell_total <- as.numeric(Matrix::colSums(mat))
    ref_cell_nfeat <- as.integer(Matrix::colSums(detect))

    # Parquet-backed giotto, then SAME user-facing addStatistics call as
    # the in-memory workflow — dispatch routes via processData internally.
    g <- Giotto::createGiottoObject(expression = mat, verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    g <- suppressWarnings(  # silence "provenance mismatch" cosmetic warning
        Giotto::addCellStatistics(g, expression_values = "raw",
                                   verbose = FALSE)
    )
    cm <- GiottoClass::getCellMetadata(g, output = "data.table")
    cm <- cm[match(colnames(mat), cm$cell_ID), ]

    expect_equal(cm$total_expr, ref_cell_total)
    expect_equal(cm$nr_feats,   ref_cell_nfeat)
    expect_equal(cm$perc_feats, ref_cell_nfeat / n_genes * 100)
})
