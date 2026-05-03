# Tests for parquetExprStore: storage class for streaming-friendly
# long-format expression matrices.

# helper: build a tiny labeled dgCMatrix
.tiny_mat <- function(n_genes = 5L, n_cells = 8L, density = 0.4, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("gene", seq_len(n_genes))
    colnames(m) <- paste0("cell", seq_len(n_cells))
    m
}

# helper: write a tiny MatrixMarket triple to a directory in 10x layout
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


# ---- Class basics ----------------------------------------------------------

test_that("parquetExprStore constructor builds a valid object", {
    pe <- parquetExprStore(
        path     = tempfile(fileext = ".parquet"),
        cell_ids = c("c1", "c2", "c3"),
        feat_ids = c("g1", "g2")
    )
    expect_s4_class(pe, "parquetExprStore")
    expect_s4_class(pe, "fileStore")
    expect_equal(pe@n_cells, 3)
    expect_equal(pe@n_genes, 2)
    expect_equal(nrow(pe), 2)        # genes
    expect_equal(ncol(pe), 3)        # cells
    expect_equal(dim(pe), c(2, 3))
    expect_equal(rownames(pe), c("g1", "g2"))
    expect_equal(colnames(pe), c("c1", "c2", "c3"))
})

test_that("constructor validates id length vs n_cells / n_genes", {
    expect_error(
        parquetExprStore(cell_ids = c("c1", "c2"), n_cells = 5),
        "length\\(cell_ids\\)"
    )
    expect_error(
        parquetExprStore(feat_ids = c("g1"), n_genes = 4),
        "length\\(feat_ids\\)"
    )
})


# ---- storeWrite (memoryMatrix) round-trip ---------------------------------

test_that("storeWrite from dgCMatrix is lossless", {
    mat <- .tiny_mat()
    pe  <- parquetExprStore(path = tempfile(fileext = ".parquet"))
    pe  <- storeWrite(pe, mat)

    expect_equal(pe@n_genes, nrow(mat))
    expect_equal(pe@n_cells, ncol(mat))
    expect_equal(pe@feat_ids, rownames(mat))
    expect_equal(pe@cell_ids, colnames(mat))

    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_setequal(names(df), c("row_id", "col_id", "value"))

    # Reconstruct and compare element-wise
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    rownames(rt) <- pe@feat_ids
    colnames(rt) <- pe@cell_ids
    expect_equal(as.matrix(mat), as.matrix(rt))
})

test_that("Parquet payload has only row_id, col_id, value (sorted by row_id)", {
    mat <- .tiny_mat(seed = 2)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    df  <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_equal(sort(names(df)), c("col_id", "row_id", "value"))
    expect_true(!is.unsorted(df$row_id))
})


# ---- storeRead returns a lazy Arrow Dataset ------------------------------

test_that("storeRead returns an Arrow Dataset (lazy, not materialized)", {
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    ds  <- storeRead(pe)
    expect_true(inherits(ds, "Dataset"))
    # Can be filtered without collecting
    sub <- ds |> dplyr::filter(row_id == 1L) |> dplyr::collect()
    expect_true(all(sub$row_id == 1L))
})


# ---- Streaming mtx_to_parquetExprStore -----------------------------------

test_that("mtx_to_parquetExprStore on a synthetic 10x triple is lossless", {
    mat     <- .tiny_mat(n_genes = 6, n_cells = 10, density = 0.5, seed = 3)
    src_dir <- .write_mtx_triple(mat, file.path(tempdir(), "tiny_10x"))
    out     <- file.path(tempdir(), "tiny_10x_out.parquet")

    pe <- mtx_to_parquetExprStore(
        mtx_path      = file.path(src_dir, "matrix.mtx"),
        barcodes_path = file.path(src_dir, "barcodes.tsv"),
        features_path = file.path(src_dir, "features.tsv"),
        output_path   = out,
        feature_id_col = 2L,
        overwrite      = TRUE
    )
    expect_s4_class(pe, "parquetExprStore")
    expect_equal(pe@n_genes, nrow(mat))
    expect_equal(pe@n_cells, ncol(mat))
    expect_equal(pe@cell_ids, colnames(mat))
    expect_equal(pe@feat_ids, rownames(mat))

    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_equal(nrow(df), Matrix::nnzero(mat))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    expect_equal(unname(as.matrix(mat)), unname(as.matrix(rt)))

    # Cleanup
    unlink(src_dir, recursive = TRUE)
    unlink(out, recursive = TRUE)
})

test_that("mtx_to_parquetExprStore writes a directory when nnz > batch_lines", {
    mat     <- .tiny_mat(n_genes = 4, n_cells = 20, density = 0.6, seed = 4)
    src_dir <- .write_mtx_triple(mat, file.path(tempdir(), "tiny_10x_b"))
    out     <- file.path(tempdir(), "tiny_10x_b_out")

    # Force directory mode by setting a tiny batch_lines
    pe <- mtx_to_parquetExprStore(
        mtx_path       = file.path(src_dir, "matrix.mtx"),
        barcodes_path  = file.path(src_dir, "barcodes.tsv"),
        features_path  = file.path(src_dir, "features.tsv"),
        output_path    = out,
        feature_id_col = 2L,
        batch_lines    = 5L,                       # forces multi-batch
        overwrite      = TRUE
    )
    expect_true(dir.exists(out))
    chunks <- list.files(out, pattern = "\\.parquet$", full.names = TRUE)
    expect_gt(length(chunks), 1L)

    # Aggregate read still works (Arrow handles directory transparently)
    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_equal(nrow(df), Matrix::nnzero(mat))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    expect_equal(unname(as.matrix(mat)), unname(as.matrix(rt)))

    unlink(src_dir, recursive = TRUE)
    unlink(out,     recursive = TRUE)
})

test_that("mtx_to_parquetExprStore rejects mismatched barcode / feature counts", {
    mat     <- .tiny_mat(n_genes = 4, n_cells = 6, seed = 5)
    src_dir <- .write_mtx_triple(mat, file.path(tempdir(), "tiny_10x_mm"))
    # corrupt barcodes file: too few lines
    writeLines(colnames(mat)[1:3], file.path(src_dir, "barcodes.tsv"))
    out <- file.path(tempdir(), "tiny_10x_mm_out.parquet")

    expect_error(
        mtx_to_parquetExprStore(
            mtx_path      = file.path(src_dir, "matrix.mtx"),
            barcodes_path = file.path(src_dir, "barcodes.tsv"),
            features_path = file.path(src_dir, "features.tsv"),
            output_path   = out,
            overwrite     = TRUE
        ),
        "number of barcodes"
    )
    unlink(src_dir, recursive = TRUE)
})


# ---- exprObj + giotto integration -----------------------------------------

test_that("parquetExprStore can be embedded in an exprObj@exprMat slot", {
    skip_if_not_installed("GiottoClass")
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    # Construct directly via new() — createExprObj would route through
    # .evaluate_expr_matrix which doesn't yet recognize parquetExprStore.
    eo  <- new("exprObj", name = "raw", exprMat = pe,
               spat_unit = "cell", feat_type = "rna")
    expect_s4_class(eo, "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
})

test_that("parquetExprStore swaps into a giotto object via setExpression", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")
    mat <- .tiny_mat(n_genes = 8, n_cells = 12, density = 0.5, seed = 9)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # Build a giotto skeleton with the in-memory matrix, then swap pe in.
    g  <- Giotto::createGiottoObject(expression = mat, verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    # @exprMat is now the parquetExprStore; structural shape matches.
    em <- GiottoClass::getExpression(g)
    expect_s4_class(slot(em, "exprMat"), "parquetExprStore")
    expect_equal(nrow(em), nrow(mat))
    expect_equal(ncol(em), ncol(mat))
    expect_equal(rownames(em), rownames(mat))
    expect_equal(colnames(em), colnames(mat))
})
