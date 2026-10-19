# storeWrite(parquetExprStore, <parquetExprStore | union>) writes cell-major:
# every file sorted by (row_id, col_id), files covering disjoint cell ranges.
# AGENTS.md states the layout and the Gram kernel requires it (adr/0018).

.layout_mat <- function(n_genes = 60L, n_cells = 500L, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    dimnames(m) <- list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
    m
}

.write_pe <- function(x) {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), x)
}

# one data.frame per written file, in file order
.files_of <- function(pe) {
    lapply(arrow::open_dataset(pe@path)$files,
        function(f) as.data.frame(arrow::read_parquet(f)))
}

# The shared check (helper-cell-major.R) asserts files cover disjoint cell
# ranges; the kernels also need each file sorted within, so check that too.
.expect_sorted_cell_major <- function(pe) {
    .expect_cell_major(pe@path)
    for (d in Filter(nrow, .files_of(pe))) {
        expect_identical(order(d$row_id, d$col_id), seq_len(nrow(d)))
        expect_false(anyDuplicated(d[c("row_id", "col_id")]) > 0L)
    }
}

.as_dense <- function(pe) {
    as.matrix(storeRead(pe, output = "dgcmatrix", max_rows = Inf, max_cols = Inf))
}

test_that("a subset written in caller order comes out cell-major and exact", {
    m <- .layout_mat()
    set.seed(2L)
    v <- .write_pe(m)[sample(60L, 25L), sample(500L, 300L)]
    out <- .write_pe(v)
    .expect_sorted_cell_major(out)
    M <- .as_dense(out)
    expect_equal(M, as.matrix(m[rownames(M), colnames(M)]))
})

test_that("several windows give the same store as one", {
    m <- .layout_mat()
    v <- .write_pe(m)[1:40, ]
    one <- .write_pe(v)
    old <- options(giottodisk.chunk_size = 37L)
    on.exit(options(old), add = TRUE)
    many <- .write_pe(v)
    expect_gt(length(arrow::open_dataset(many@path)$files), 1L)
    .expect_sorted_cell_major(many)
    expect_identical(.as_dense(many), .as_dense(one))
})

test_that("a union writes cell-major across its substores", {
    m <- .layout_mat()
    u <- unionParquetExprStore(list(.write_pe(m[, 1:200]), .write_pe(m[, 201:500])))
    old <- options(giottodisk.chunk_size = 90L)
    on.exit(options(old), add = TRUE)
    out <- .write_pe(u[5:30, ])
    .expect_sorted_cell_major(out)
    M <- .as_dense(out)
    expect_equal(M, as.matrix(m[rownames(M), colnames(M)]))
})

test_that("the written layout is the same on every run", {
    m <- .layout_mat()
    v <- .write_pe(m)[sample(60L, 30L), ]
    a <- .files_of(.write_pe(v))
    b <- .files_of(.write_pe(v))
    expect_identical(a, b)
})

test_that("a view with no stored values still writes a readable store", {
    m <- .layout_mat()
    m[, 1:10] <- 0
    m <- Matrix::drop0(m)
    out <- .write_pe(.write_pe(m)[, 1:10])
    expect_equal(sum(.as_dense(out)), 0)
})
