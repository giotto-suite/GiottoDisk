# Every writer of an expression store must lay it out cell-major: each part
# file covers its own range of cells (helper-cell-major.R). Each test forces
# several parts, checks the layout, and checks the values against the source
# matrix. The Stereo-seq GEF readers are covered in test-stereoseq-gef.R.

.cm_mat <- function(n_genes = 6L, n_cells = 40L, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.5,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("gene", seq_len(n_genes))
    colnames(m) <- paste0("cell", seq_len(n_cells))
    m
}

# store -> genes x cells matrix, columns in store order
.cm_as_matrix <- function(pe) {
    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    Matrix::sparseMatrix(i = df$col_id, j = df$row_id, x = df$value,
                         dims = c(nrow(pe), ncol(pe)))
}

.expect_same_values <- function(pe, mat) {
    expect_equal(unname(as.matrix(.cm_as_matrix(pe))),
                 unname(as.matrix(mat[pe@feat_ids, pe@cell_ids])))
}


test_that("storeWrite of an in-memory matrix is cell-major", {
    mat <- .cm_mat(seed = 2L)
    out <- tempfile("cm_mem_")
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    pe <- storeWrite(parquetExprStore(path = out), mat)
    .expect_cell_major(out)
    .expect_same_values(pe, mat)
})

test_that("mtxInput writes cell-major parts", {
    mat <- .cm_mat(seed = 3L)
    src <- tempfile("cm_mtx_src_")
    out <- tempfile("cm_mtx_")
    on.exit(unlink(c(src, out), recursive = TRUE), add = TRUE)
    dir.create(src)
    writeLines(colnames(mat), file.path(src, "barcodes.tsv"))
    utils::write.table(data.frame(id = paste0("ENSG", seq_len(nrow(mat))),
                                  name = rownames(mat), type = "Gene Expression"),
                       file.path(src, "features.tsv"), sep = "\t", quote = FALSE,
                       col.names = FALSE, row.names = FALSE)
    Matrix::writeMM(mat, file.path(src, "matrix.mtx"))

    inp <- mtxInput(mtx_path = file.path(src, "matrix.mtx"), batch_lines = 17L)
    pe <- storeWrite(parquetExprStore(path = out), inp)
    # batches are cut by line count, so a cell can span two neighbouring parts
    expect_gt(nrow(.expect_cell_major(out, strict = FALSE)), 1L)
    .expect_same_values(pe, mat)
})

test_that("tenxH5Input writes cell-major parts", {
    skip_if_not_installed("hdf5r")
    mat <- .cm_mat(seed = 4L)
    h5p <- tempfile(fileext = ".h5")
    out <- tempfile("cm_h5_")
    on.exit(unlink(c(h5p, out), recursive = TRUE), add = TRUE)
    csc <- methods::as(mat, "CsparseMatrix")          # genes x cells, CSC by cell
    h5 <- hdf5r::H5File$new(h5p, mode = "w")
    g <- h5$create_group("matrix")
    g[["data"]] <- as.integer(csc@x)
    g[["indices"]] <- as.integer(csc@i)
    g[["indptr"]] <- as.integer(csc@p)
    g[["shape"]] <- as.integer(dim(mat))
    g[["barcodes"]] <- colnames(mat)
    f <- g$create_group("features")
    f[["id"]] <- paste0("ENSG", seq_len(nrow(mat)))
    f[["name"]] <- rownames(mat)
    f[["feature_type"]] <- rep("Gene Expression", nrow(mat))
    h5$close_all()

    pe <- storeWrite(parquetExprStore(path = out), tenxH5Input(h5p, batch_cells = 9L))
    expect_gt(nrow(.expect_cell_major(out)), 1L)
    .expect_same_values(pe, mat)
})

test_that("tenxZarrInput writes cell-major parts in cellblock mode", {
    skip_if_no_zarr_deps()
    fx <- make_zarr_fixture()
    out <- tempfile("cm_zarr_")
    on.exit(unlink(c(fx$dir, out), recursive = TRUE), add = TRUE)

    inp <- tenxZarrInput(fx$paths$cell_feature_matrix, mode = "cellblock",
                         cells_per_block = 5L)
    pe <- storeWrite(parquetExprStore(path = out), inp)
    expect_gt(nrow(.expect_cell_major(out)), 1L)
    m <- fx$truth$cfm
    expect_equal(sum(.cm_as_matrix(pe)), sum(m))
})

test_that("csvWideInput writes cell-major parts", {
    mat <- .cm_mat(seed = 5L)
    csv <- tempfile(fileext = ".csv")
    out <- tempfile("cm_csv_")
    on.exit(unlink(c(csv, out), recursive = TRUE), add = TRUE)
    wide <- data.frame(cell_ID = colnames(mat), t(as.matrix(mat)), check.names = FALSE)
    utils::write.csv(wide, csv, row.names = FALSE, quote = FALSE)

    pe <- storeWrite(parquetExprStore(path = out), csvWideInput(csv, batch_rows = 7L))
    expect_gt(nrow(.expect_cell_major(out)), 1L)
    .expect_same_values(pe, mat)
})

test_that("a store baked through its op chain is written cell-major", {
    mat <- .cm_mat(n_cells = 60L, seed = 6L)
    src <- tempfile("cm_bake_src_")
    out <- tempfile("cm_bake_")
    on.exit(unlink(c(src, out), recursive = TRUE), add = TRUE)
    pe <- storeWrite(parquetExprStore(path = src), mat)
    # a post-phase op forces the windowed bake rather than a lazy rewrite
    pe_post <- .pe_push_op(pe, list(type = "log", base = 2), phase = "post")
    withr::local_options(list(giottodisk.chunk_size = 11L))
    baked <- storeWrite(parquetExprStore(path = out), pe_post)
    expect_gt(nrow(.expect_cell_major(out)), 1L)
    .expect_same_values(baked, log2(mat + 1))
})
