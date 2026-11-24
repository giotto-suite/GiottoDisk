# zscoreScaleParam / defaultScaleParam on parquetExprBase match Giotto's
# in-memory scaling, read as a bounded slice

.scale_fixture <- function(seed = 4L) {
    mat <- .rand_counts(n_genes = 12L, n_cells = 30L, seed = seed)
    mat[5L, ] <- 0          # a feature with nothing stored: NaN, as in memory
    mat <- Matrix::drop0(mat)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    norm <- list(Giotto::normParam("library", scalefactor = 1e3),
        Giotto::normParam("log", base = 2, offset = 1))
    list(mat = processData(mat, norm), pe = processData(pe, norm))
}

.read_all <- function(pe) .bare(storeRead(pe, output = "dgcmatrix",
    max_rows = Inf, max_cols = Inf))

# Values and shape only (base::scale() adds attributes)
.bare <- function(m) {
    m <- as.matrix(m)
    attributes(m) <- list(dim = dim(m))
    m
}

test_that("zscoreScaleParam matches the in-memory result on either margin", {
    skip_if_not_installed("Giotto")
    fx <- .scale_fixture()
    for (margin in 1:2) {
        p <- Giotto::scaleParam("zscore", MARGIN = margin)
        pe <- processData(fx$pe, p)
        types <- vapply(c(pe@ops, pe@post_ops), `[[`, "", "type")
        expect_identical(utils::tail(types, 2L), c("multiply", "add"))
        expect_identical(pe@post_ops[[length(pe@post_ops)]]$type, "add")
        ref <- .bare(processData(fx$mat, p))
        expect_equal(.read_all(pe), ref, , label = paste("MARGIN", margin))
    }
})

test_that("scaling without centring stays sparse and matches base::scale()", {
    skip_if_not_installed("Giotto")
    fx <- .scale_fixture()
    # A feature with nothing stored stays 0 (in memory: NaN), being sparse
    keep <- c(1:4, 6:12)
    fx$mat <- fx$mat[keep, ]
    fx$pe <- fx$pe[keep, ]
    pe <- processData(fx$pe, Giotto::scaleParam("zscore", center = FALSE, MARGIN = 1))
    expect_false("add" %in% vapply(c(pe@ops, pe@post_ops), `[[`, "", "type"))
    ref <- .bare(t(scale(t(as.matrix(fx$mat)), center = FALSE)))
    expect_equal(.read_all(pe), ref)
    # Sparse, so streaming passes still read it
    expect_equal(.stream_expr_accum(pe, axis = "feat", stats = "sum")$sum,
        rowSums(ref))
})

test_that("defaultScaleParam matches in-memory in both orders", {
    skip_if_not_installed("Giotto")
    fx <- .scale_fixture(seed = 7L)
    keep <- c(1:4, 6:12)           # in memory, a NaN feature spreads to every cell
    fx$mat <- fx$mat[keep, ]
    fx$pe <- fx$pe[keep, ]
    for (ord in c("first_feats", "first_cells")) {
        p <- Giotto::scaleParam("default", scale_order = ord, verbose = FALSE)
        pe <- processData(fx$pe, p)
        expect_equal(.read_all(pe), .bare(processData(fx$mat, p)),
            , tolerance = 1e-10, label = ord)
    }
})

test_that("the statistics are frozen: a slice reads the values of the whole", {
    skip_if_not_installed("Giotto")
    fx <- .scale_fixture()
    p <- Giotto::scaleParam("zscore", MARGIN = 1)
    pe <- processData(fx$pe, p)
    ref <- .bare(processData(fx$mat, p))
    expect_equal(.bare(storeRead(pe[2:4, 6:20], output = "dgcmatrix")), ref[2:4, 6:20])
    # Scaling a slice computes its own statistics
    sub <- processData(fx$pe[2:4, 6:20], p)
    expect_equal(.bare(storeRead(sub, output = "dgcmatrix")),
        .bare(processData(fx$mat[2:4, 6:20], p)))
})

test_that("a centred store is read only as a bounded slice", {
    skip_if_not_installed("Giotto")
    fx <- .scale_fixture()
    pe <- processData(fx$pe, Giotto::scaleParam("zscore", MARGIN = 1))
    expect_error(storeRead(pe, output = "tibble"), "bounded slice")
    expect_error(.stream_expr_accum(pe, axis = "feat"), "bounded slice")
    expect_error(storeRead(pe, output = "dgcmatrix", max_rows = Inf, max_cols = Inf),
        NA)
    withr::local_options(list(giottodisk.dense_max = 50))
    expect_error(storeRead(pe[1:2, 1:30], output = "dgcmatrix"), "dense_max")
    expect_error(storeRead(pe[1:2, 1:20], output = "dgcmatrix"), NA)
    # Two centring steps in a row cannot be streamed
    pe2 <- processData(pe, Giotto::scaleParam("zscore", MARGIN = 2))
    expect_error(processData(pe2, Giotto::scaleParam("zscore", MARGIN = 1)),
        "storeWrite")
})

test_that("processExpression scales an exprObj on disk through the list method", {
    skip_if_not_installed("Giotto")
    skip_if_not(methods::existsMethod(GiottoClass::processData, c("ANY", "list")))
    fx <- .scale_fixture()
    ex <- GiottoClass::createExprObj(fx$pe, name = "normalized", spat_unit = "cell",
        feat_type = "rna")
    out <- processData(ex, list(Giotto::scaleParam("zscore", MARGIN = 1)), name = "scaled")
    expect_s4_class(out[], "parquetExprStore")
    expect_equal(.read_all(out[]), .bare(processData(fx$mat,
        Giotto::scaleParam("zscore", MARGIN = 1))))
})
