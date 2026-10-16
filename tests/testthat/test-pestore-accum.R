# The `max` accumulator, `by_value` grouping keys, and the per-feature value
# histogram built on them. Ground truth is the same quantity computed on the
# in-memory dgCMatrix.

.acc_mat <- function(n_genes = 15L, n_cells = 40L, seed = 1L, prefix = "c") {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.4,
        rand.x = function(n) as.double(rpois(n, 6L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0(prefix, seq_len(n_cells))
    m
}
.acc_store <- function(m) {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}
# per-row max over stored entries (0 for a row with none)
.row_max <- function(m) {
    m <- methods::as(m, "TsparseMatrix")
    out <- numeric(nrow(m))
    if (length(m@x)) {
        r <- tapply(m@x, m@i + 1L, max)
        out[as.integer(names(r))] <- r
    }
    out
}
# in-memory reference for .pe_value_hist()
.ref_hist <- function(m, n_bins) {
    mx <- .row_max(m)
    t <- methods::as(m, "TsparseMatrix")
    keep <- t@x > 0
    g <- t@i[keep] + 1L
    x <- t@x[keep] * (n_bins / mx[g])
    b <- pmin(ceiling(x) - 1, n_bins - 1L)
    d <- data.table::data.table(pos = g, bin = as.integer(b))[
        , list(n = as.numeric(.N)), by = c("pos", "bin")]
    data.table::setorder(d, pos, bin)
    d
}
.sorted <- function(d) {
    d <- data.table::copy(d)
    data.table::setorder(d, pos, bin)
    d
}

test_that("max matches the in-memory per-feature and per-cell max", {
    m  <- .acc_mat()
    pe <- .acc_store(m)
    expect_equal(.stream_expr_accum(pe, stats = "max")$max, .row_max(m))
    expect_equal(.stream_expr_accum(pe, axis = "cell", stats = "max")$max,
                 .row_max(Matrix::t(m)))
})

test_that("max folds by max across substores and R-side windows", {
    m1 <- .acc_mat(n_cells = 9L, seed = 2L, prefix = "a")
    m2 <- .acc_mat(n_cells = 6L, seed = 3L, prefix = "b")
    u <- unionParquetExprStore(list(.acc_store(m1), .acc_store(m2)))
    expect_equal(.stream_expr_accum(u, stats = "max")$max,
                 .row_max(cbind(m1, m2)))

    pe <- .acc_store(m1)
    logged <- .pe_push_op(pe, list(type = "log", base = 2), phase = "post")
    withr::with_options(list(giottodisk.chunk_size = 4L), {
        got <- .stream_expr_accum(logged, stats = "max")$max
    })
    expect_equal(got, log2(.row_max(m1) + 1))
})

test_that("by_value keys group the aggregate by a function of value", {
    m  <- .acc_mat()
    pe <- .acc_store(m)
    got <- .pe_accum_raw(pe, stats = "nnz",
        by_value = list(hi = rlang::expr(as.integer(value > 6))))
    t <- methods::as(m, "TsparseMatrix")
    ref <- data.table::data.table(pos = t@i + 1L, hi = as.integer(t@x > 6))[
        , list(nnz = as.numeric(.N)), by = c("pos", "hi")]
    data.table::setorder(ref, pos, hi)
    data.table::setorder(got, pos, hi)
    expect_equal(got$nnz, ref$nnz)
    expect_equal(as.integer(got$hi), ref$hi)
})

test_that("the value histogram matches the in-memory binning", {
    m  <- .acc_mat()
    pe <- .acc_store(m)
    h <- .pe_value_hist(pe, n_bins = 8L)
    expect_equal(.sorted(h$counts), .ref_hist(m, 8L))
    expect_equal(h$max, .row_max(m))
    expect_equal(h$n_zero, unname(ncol(m) - Matrix::rowSums(m > 0)))
    expect_length(pe@ops, 0L)
})

test_that("the value histogram is the same over windows, unions and post ops", {
    m1 <- .acc_mat(n_cells = 9L, seed = 2L, prefix = "a")
    m2 <- .acc_mat(n_cells = 6L, seed = 3L, prefix = "b")
    u <- unionParquetExprStore(list(.acc_store(m1), .acc_store(m2)))
    withr::with_options(list(giottodisk.chunk_size = 4L), {
        h <- .pe_value_hist(u, n_bins = 5L)
    })
    expect_equal(.sorted(h$counts), .ref_hist(cbind(m1, m2), 5L))

    pe <- .acc_store(m1)
    logged <- .pe_push_op(pe, list(type = "log", base = 2), phase = "post")
    lm <- m1
    lm@x <- log2(lm@x + 1)
    expect_equal(.sorted(.pe_value_hist(logged, n_bins = 6L)$counts),
                 .ref_hist(lm, 6L))
})
