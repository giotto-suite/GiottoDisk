# Compare group methods on parquetExprBase: `x >= t` queues a `compare`
# indicator record, so the margin methods count through the accumulator. The
# ground truth throughout is the same comparison on the in-memory dgCMatrix.

.cmp_mat <- function(n_genes = 12L, n_cells = 30L, seed = 1L,
                     prefix = "c") {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.5,
        rand.x = function(n) as.double(rpois(n, 3L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0(prefix, seq_len(n_cells))
    m
}
.cmp_store <- function(m) {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}
.cmp_vals <- function(v) unname(as.numeric(v))

test_that("every sparse-safe comparison matches the in-memory margins", {
    m  <- .cmp_mat()
    pe <- .cmp_store(m)
    cases <- list(
        list(">=", 3), list(">", 2), list("<", 0), list("<=", -1),
        list("==", 4), list("!=", 0)
    )
    for (cs in cases) {
        op <- cs[[1L]]; t <- cs[[2L]]
        got <- do.call(op, list(pe, t))
        ref <- do.call(op, list(m, t))
        # guard against a vacuous pass: the chosen thresholds must split
        if (!op %in% c("<", "<=")) expect_gt(sum(ref), 0)
        expect_equal(.cmp_vals(rowSums(got)), .cmp_vals(Matrix::rowSums(ref)),
                     label = paste("rowSums x", op, t))
        expect_equal(.cmp_vals(colSums(got)), .cmp_vals(Matrix::colSums(ref)),
                     label = paste("colSums x", op, t))
    }
})

test_that("the comparison is lazy and reads back as a 0/1 matrix", {
    m  <- .cmp_mat()
    pe <- .cmp_store(m)
    got <- pe >= 3
    expect_s4_class(got, "parquetExprStore")
    expect_equal(got@ops[[length(got@ops)]],
                 list(type = "compare", op = ">=", axis = "all", e2 = 3))

    # a per-axis record is reserved, not executable yet
    feat <- pe
    feat@ops <- list(list(type = "compare", op = ">=", axis = "feat",
                          e2 = 3))
    expect_error(rowSums(feat), "only `axis = \"all\"`")

    M <- storeRead(got, output = "dgcmatrix", max_rows = Inf, max_cols = Inf)
    expect_equal(as.matrix(M), as.matrix(1 * (m >= 3)),
                 ignore_attr = TRUE)
})

test_that("a scalar on the left is flipped onto the store", {
    m  <- .cmp_mat()
    pe <- .cmp_store(m)
    expect_equal((2 < pe)@ops, (pe > 2)@ops)
    expect_equal((0 > pe)@ops, (pe < 0)@ops)
    expect_equal(.cmp_vals(rowSums(4 == pe)),
                 .cmp_vals(Matrix::rowSums(m == 4)))
})

test_that("a comparison TRUE at zero is refused, not densified", {
    pe <- .cmp_store(.cmp_mat())
    for (f in list(function(x) x >= 0, function(x) x < 1,
                   function(x) x == 0, function(x) x != 2,
                   function(x) 0 <= x)) {
        expect_error(f(pe), "TRUE for every unstored zero")
    }
    expect_error(pe >= c(1, 2), "single non-NA number")
    expect_error(pe >= NA_real_, "single non-NA number")
    expect_error(pe >= pe, "two expression stores")
})

test_that("the comparison sees the values the chain produces", {
    m  <- .cmp_mat()
    pe <- .cmp_store(m)

    # lazy phase: a queued multiply is applied before the threshold
    halved <- .pe_push_op(pe,
        list(type = "multiply", axis = "all", factors = 0.5), phase = "lazy")
    expect_equal(.cmp_vals(rowSums(halved >= 2)),
                 .cmp_vals(Matrix::rowSums(m * 0.5 >= 2)))

    # post phase: after a post op the record has to follow it R-side
    logged <- .pe_push_op(pe, list(type = "log", base = 2), phase = "post")
    got <- logged > 1.5
    expect_length(got@ops, 0L)
    expect_equal(got@post_ops[[2L]]$type, "compare")
    ref <- log2(as.matrix(m) + 1) > 1.5
    expect_equal(.cmp_vals(rowSums(got)), .cmp_vals(rowSums(ref)))
    expect_equal(.cmp_vals(colSums(got)), .cmp_vals(colSums(ref)))
})

test_that("comparisons compose with subsetting, transpose and unions", {
    m  <- .cmp_mat()
    pe <- .cmp_store(m)

    # narrowing before or after the comparison is the same view
    ref <- Matrix::rowSums(m[3:8, 5:20] >= 3)
    expect_equal(.cmp_vals(rowSums(pe[3:8, 5:20] >= 3)), .cmp_vals(ref))
    expect_equal(.cmp_vals(rowSums((pe >= 3)[3:8, 5:20])), .cmp_vals(ref))

    # transposed: the logical row margin is the cell margin
    expect_equal(.cmp_vals(rowSums(t(pe) >= 3)),
                 .cmp_vals(Matrix::colSums(m >= 3)))

    m1 <- .cmp_mat(n_cells = 7L, seed = 2L, prefix = "a")
    m2 <- .cmp_mat(n_cells = 5L, seed = 3L, prefix = "b")
    u  <- unionParquetExprStore(list(.cmp_store(m1), .cmp_store(m2)))
    mm <- cbind(m1, m2)
    expect_equal(.cmp_vals(rowSums(u >= 3)),
                 .cmp_vals(Matrix::rowSums(mm >= 3)))
    expect_equal(.cmp_vals(colSums(u >= 3)),
                 .cmp_vals(Matrix::colSums(mm >= 3)))
})

test_that("the duckdb carrier returns the same indicator", {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("dbplyr")
    m  <- .cmp_mat()
    pe <- .cmp_store(m) >= 3
    via_arrow <- data.table::as.data.table(
        dplyr::collect(storeRead(pe, output = "query")))
    via_duck <- data.table::as.data.table(
        dplyr::collect(storeRead(pe, output = "duckdb")))
    data.table::setorder(via_arrow, row_id, col_id)
    data.table::setorder(via_duck, row_id, col_id)
    expect_equal(nrow(via_duck), sum(m >= 3))
    expect_equal(via_duck$value, via_arrow$value)
    expect_true(all(via_duck$value == 1))
})

test_that("the flex margins reach the store through the comparison", {
    m  <- .cmp_mat()
    pe <- .cmp_store(m)
    expect_equal(.cmp_vals(GiottoClass::rowSums_flex(pe >= 1)),
                 .cmp_vals(Matrix::rowSums(m >= 1)))
    expect_equal(.cmp_vals(GiottoClass::colSums_flex(pe >= 1)),
                 .cmp_vals(Matrix::colSums(m >= 1)))
})
