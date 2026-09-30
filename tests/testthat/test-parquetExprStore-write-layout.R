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

# Fanned out (adr/0019), the windows are written by workers without GiottoDisk
# and renamed into place; the store must be the one the serial loop writes.
.write_pe_with <- function(x, workers) {
    old <- options(giottodisk.par_workers = workers)
    on.exit(options(old), add = TRUE)
    .write_pe(x)
}

.expect_fanout_matches_serial <- function(v) {
    old <- options(giottodisk.chunk_size = 37L)
    on.exit(options(old), add = TRUE)
    # otherwise the write falls back to serial and matches trivially
    wins <- .pe_windows(v, 37L)
    expect_gt(length(wins), 1L)
    for (d in wins) expect_false(is.null(.pe_lower_read(.pe_window_store(d))))
    old_w <- options(giottodisk.par_workers = 2L)
    expect_identical(.isolated_workers(), 2L)
    options(old_w)
    serial <- .write_pe_with(v, 1L)
    fanned <- .write_pe_with(v, 2L)
    expect_gt(length(arrow::open_dataset(serial@path)$files), 1L)
    expect_identical(basename(arrow::open_dataset(fanned@path)$files),
                     basename(arrow::open_dataset(serial@path)$files))
    expect_identical(.files_of(fanned), .files_of(serial))
    .expect_sorted_cell_major(fanned)
    invisible(fanned)
}

# a multiply payload keyed by on-disk id plus a log, so both the payload join
# and a plain transform run in the worker
.with_ops <- function(pe) {
    set.seed(3L)
    pe@ops <- list(
        list(type = "multiply", axis = "cell",
             factors = stats::setNames(list(stats::runif(pe@n_cells)), pe@uid)),
        list(type = "log", base = 2))
    pe
}

test_that("a fanned-out write gives the serial store", {
    skip_on_cran()
    skip_if_not_installed("mirai")
    skip_if_not_installed("codetools")
    m <- .layout_mat()
    set.seed(4L)
    v <- .with_ops(.write_pe(m))[sample(60L, 25L), 20:480]
    .expect_fanout_matches_serial(v)
})

test_that("a fanned-out union write gives the serial store", {
    skip_on_cran()
    skip_if_not_installed("mirai")
    skip_if_not_installed("codetools")
    m <- .layout_mat()
    u <- unionParquetExprStore(list(.write_pe(m[, 1:200]), .write_pe(m[, 201:500])))
    .expect_fanout_matches_serial(u[5:30, ])
})

test_that("empty windows keep the serial file numbering when fanned out", {
    skip_on_cran()
    skip_if_not_installed("mirai")
    skip_if_not_installed("codetools")
    m <- .layout_mat()
    m[, 1:120] <- 0
    m <- Matrix::drop0(m)
    .expect_fanout_matches_serial(.write_pe(m)[, 1:300])
})

test_that("a store with its own reader is written serially, unchanged", {
    m <- .layout_mat()
    v <- .write_pe(m)[1:40, ]
    v@read_fun <- function(x, ...) arrow::open_dataset(x, ...)
    expect_null(.pe_lower_read(v))
    old <- options(giottodisk.chunk_size = 37L)
    on.exit(options(old), add = TRUE)
    expect_identical(.files_of(.write_pe_with(v, 2L)), .files_of(.write_pe_with(v, 1L)))
})
