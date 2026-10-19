# processData(parquetExprBase, thresholdParam("binarize")) against Giotto's
# in-memory method on the same matrix.

.bz_mat <- function(G = 10L, N = 50L, seed = 6L, prefix = "c") {
    set.seed(seed)
    m <- Matrix::rsparsematrix(G, N, density = 0.5,
        rand.x = function(k) as.double(rpois(k, 6L) + 1L))
    dimnames(m) <- list(paste0("g", seq_len(G)), paste0(prefix, seq_len(N)))
    m
}
.bz_pe <- function(m) {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}
.bz_read <- function(p) {
    as.matrix(storeRead(p, output = "dgcmatrix", max_rows = Inf, max_cols = Inf))
}
.bz_ref <- function(m, p) as.matrix(GiottoClass::processData(m, p))

test_that("numeric thresholds binarize lazily, one or per feature", {
    m <- .bz_mat()
    pe <- .bz_pe(m)
    for (t in list(3, seq(2, 8, length.out = nrow(m)))) {
        p <- Giotto::thresholdParam("binarize", threshold = t)
        got <- GiottoClass::processData(pe, p)
        expect_length(got@ops, 1L)
        expect_equal(.bz_read(got), .bz_ref(m, p), ignore_attr = TRUE)
    }
    expect_length(pe@ops, 0L)
})

test_that("rules on integer counts match the exact in-memory split", {
    m <- .bz_mat()
    pe <- .bz_pe(m)
    for (p in list(Giotto::thresholdParam("binarize", threshold = "kmeans"),
                   Giotto::thresholdParam("binarize", threshold = "rank",
                                          percentage_rank = 20))) {
        expect_equal(.bz_read(GiottoClass::processData(pe, p)), .bz_ref(m, p),
                     ignore_attr = TRUE)
    }
})

test_that("rules apply per feature on unions and transposed stores", {
    m <- .bz_mat()
    p <- Giotto::thresholdParam("binarize", threshold = "kmeans")
    u <- unionParquetExprStore(list(.bz_pe(m[, 1:20]), .bz_pe(m[, 21:50])))
    expect_equal(.bz_read(GiottoClass::processData(u, p)), .bz_ref(m, p),
                 ignore_attr = TRUE)
    tp <- GiottoClass::processData(t(.bz_pe(m)), p)
    expect_equal(.bz_read(tp), t(.bz_ref(m, p)), ignore_attr = TRUE)
})

test_that("rules on continuous values agree with the exact split almost everywhere", {
    m <- .bz_mat(G = 20L, N = 200L)
    m@x <- log1p(m@x * runif(length(m@x), 0.5, 1.5))
    pe <- .bz_pe(m)
    p <- Giotto::thresholdParam("binarize", threshold = "kmeans")
    expect_gt(mean(.bz_read(GiottoClass::processData(pe, p)) == .bz_ref(m, p)),
              0.99)
})
