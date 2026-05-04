# Tests for createGiottoFromParquet (Phase 3 Layer 1)

.tiny_mat <- function(n_genes = 12L, n_cells = 30L, density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}

# Helper: write a tiny 10x triple to a directory
.write_mtx_triple <- function(mat, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(colnames(mat), file.path(out_dir, "barcodes.tsv"))
    feats <- data.frame(
        ensg   = paste0("ENSG", sprintf("%011d", seq_len(nrow(mat)))),
        symbol = rownames(mat),
        type   = "Gene Expression"
    )
    write.table(feats, file.path(out_dir, "features.tsv"),
        sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    Matrix::writeMM(mat, file.path(out_dir, "matrix.mtx"))
    out_dir
}


test_that("createGiottoFromParquet errors when neither path nor mtx_dir supplied", {
    expect_error(createGiottoFromParquet(),
                 "must be supplied")
})


test_that("createGiottoFromParquet errors when both path and mtx_dir supplied", {
    expect_error(
        createGiottoFromParquet(parquet_path = tempfile(), mtx_dir = tempdir()),
        "not both"
    )
})


test_that("createGiottoFromParquet from existing parquet builds a giotto", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
                      mat)
    g <- createGiottoFromParquet(
        parquet_path = pe@path,
        cell_ids     = colnames(mat),
        feat_ids     = rownames(mat),
        spatial_locs = NULL,
        chunk_size   = 1000L,    # skip sc_recommend_chunk
        verbose      = FALSE
    )
    expect_s4_class(g, "giotto")
    expr <- GiottoClass::getExpression(g, output = "exprObj")
    expect_s4_class(slot(expr, "exprMat"), "parquetExprStore")
    expect_equal(nrow(expr), nrow(mat))
    expect_equal(ncol(expr), ncol(mat))
})


test_that("createGiottoFromParquet from mtx_dir streams without dgCMatrix", {
    skip_if_not_installed("Giotto")
    mat     <- .tiny_mat(n_genes = 8, n_cells = 15, density = 0.6, seed = 7)
    src_dir <- .write_mtx_triple(mat, file.path(tempdir(),
                                                  "tiny_10x_layer1"))
    out_pq  <- file.path(tempdir(), "tiny_10x_layer1_out.parquet")

    g <- createGiottoFromParquet(
        mtx_dir      = src_dir,
        output_path  = out_pq,
        chunk_size   = 1000L,
        verbose      = FALSE
    )
    expect_s4_class(g, "giotto")
    pe <- slot(GiottoClass::getExpression(g, output = "exprObj"), "exprMat")
    expect_s4_class(pe, "parquetExprStore")
    expect_equal(nrow(pe), nrow(mat))
    expect_equal(ncol(pe), ncol(mat))
    expect_equal(pe@feat_ids, rownames(mat))
    expect_equal(pe@cell_ids, colnames(mat))

    # Element-wise round-trip
    df <- dplyr::collect(storeRead(pe))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    expect_equal(unname(as.matrix(mat)), unname(as.matrix(rt)))

    unlink(src_dir, recursive = TRUE)
    unlink(out_pq,  recursive = TRUE)
})


test_that("createGiottoFromParquet requires output_path when mtx_dir supplied", {
    expect_error(
        createGiottoFromParquet(mtx_dir = tempdir()),
        "output_path` is required"
    )
})


test_that("createGiottoFromParquet requires cell_ids/feat_ids when parquet_path supplied", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
                      mat)
    expect_error(
        createGiottoFromParquet(parquet_path = pe@path),
        "cell_ids` and `feat_ids` are"
    )
})


test_that("sc_recommend_chunk returns a sensible integer", {
    rec <- sc_recommend_chunk(n_cells = 10000, n_genes = 1000, density = 0.1,
                                verbose = FALSE)
    expect_true(is.integer(rec))
    expect_true(rec >= 10000)   # floor
})
