# hnswKNN: approximate kNN search returning a dbscan-compatible shape

.nn_mat <- function(n = 500L, d = 10L, seed = 1L) {
    set.seed(seed)
    matrix(stats::rnorm(n * d), nrow = n, ncol = d)
}

test_that("hnswKNN returns the same shape as dbscan::kNN", {
    skip_if_not_installed("RcppHNSW")
    skip_if_not_installed("dbscan")
    m <- .nn_mat()
    k <- 15L

    hn <- hnswKNN(m, k = k, n_threads = 2L)
    ex <- dbscan::kNN(m, k = k, sort = TRUE)

    expect_identical(class(hn), class(ex))
    expect_true(all(c("id", "dist", "k", "sort", "metric") %in% names(hn)))
    expect_identical(dim(hn$id), dim(ex$id))
    expect_identical(dim(hn$dist), dim(ex$dist))
    expect_identical(hn$k, k)
    expect_type(hn$id, "integer")
})

test_that("hnswKNN excludes self and returns sorted distances", {
    skip_if_not_installed("RcppHNSW")
    m <- .nn_mat()
    hn <- hnswKNN(m, k = 10L, n_threads = 2L)

    expect_false(any(hn$id == seq_len(nrow(m))))
    expect_true(all(apply(hn$dist, 1L, function(r) !is.unsorted(r))))
})

test_that("hnswKNN drops the correct entry when self is not column 1", {
    # Duplicated coordinates put a zero-distance twin alongside the self-hit,
    # so the self can land at any column. A column-major mask extract
    # misaligns rows here while looking correct when self is always first.
    skip_if_not_installed("RcppHNSW")
    m <- .nn_mat(n = 100L, d = 5L)
    md <- rbind(m, m) # every row has an exact duplicate

    hn <- hnswKNN(md, k = 5L, n_threads = 2L)

    expect_false(any(hn$id == seq_len(nrow(md))))
    expect_identical(dim(hn$id), c(nrow(md), 5L))
    expect_false(anyNA(hn$id))
})

test_that("hnswKNN recall is high against exact search", {
    skip_if_not_installed("RcppHNSW")
    skip_if_not_installed("dbscan")
    m <- .nn_mat(n = 1000L, d = 15L)
    k <- 20L

    hn <- hnswKNN(m, k = k, n_threads = 2L)
    ex <- dbscan::kNN(m, k = k, sort = TRUE)
    recall <- mean(vapply(seq_len(nrow(m)), function(i) {
        length(intersect(hn$id[i, ], ex$id[i, ])) / k
    }, numeric(1L)))

    expect_gt(recall, 0.95)
})

test_that("dbscan::sNN consumes an hnswKNN result", {
    skip_if_not_installed("RcppHNSW")
    skip_if_not_installed("dbscan")
    m <- .nn_mat()
    hn <- hnswKNN(m, k = 15L, n_threads = 2L)

    snn <- dbscan::sNN(x = hn, k = 15L, kt = NULL)
    expect_true(all(c("shared", "id", "dist") %in% names(snn)))
    expect_identical(dim(snn$shared), dim(hn$id))
})

test_that("hnswKNN rejects k >= nrow(x)", {
    skip_if_not_installed("RcppHNSW")
    m <- .nn_mat(n = 20L)
    expect_error(hnswKNN(m, k = 20L), "must be less than")
})
