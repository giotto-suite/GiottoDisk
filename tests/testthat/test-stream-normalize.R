# Tests for streaming normalize dispatch:
#   processData(parquetExprStore, libraryNormParam) -> appends norm_libsize
#       op (log = FALSE) to pe@post_ops
#   processData(parquetExprStore, logNormParam)     -> flips the same op's
#       log flag to TRUE (libsize+log fuse into one op)

.tiny_mat <- function(n_genes = 12L, n_cells = 30L,
                       density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}


test_that("libraryNormParam appends norm_libsize op to @post_ops", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    pe2 <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))

    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_length(pe2@post_ops, 1L)
    expect_equal(pe2@post_ops[[1]]$type, "norm_libsize")
    expect_equal(pe2@post_ops[[1]]$scalef$scalef, 1e4 / libsz)
    expect_equal(pe2@post_ops[[1]]$scalef$orig_row_id, seq_along(libsz))
    # scaling carries no log state -- that is a separate op
    expect_null(pe2@post_ops[[1]]$log)
})


test_that("logNormParam appends an independent log op", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 2)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))
    pe  <- GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 1))

    # Two records in chain order, not one fused record.
    expect_length(pe@post_ops, 2L)
    expect_equal(vapply(pe@post_ops, function(o) o$type, character(1L)),
                 c("norm_libsize", "log"))
    expect_equal(pe@post_ops[[2]]$base, 2)

    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(pe@post_ops[[1]]$scalef$scalef, 1e4 / libsz)
})


# log-only on raw counts used to be rejected outright: the log flag lived on
# the libsize record, so there was nowhere to put it without one. As its own
# op it stands alone.
test_that("logNormParam works with no library normalization", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 3)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 1))

    expect_length(pe@post_ops, 1L)
    expect_equal(pe@post_ops[[1]]$type, "log")

    got <- storeRead(pe, output = "dgcmatrix")
    expect_equal(as.matrix(got[rownames(mat), colnames(mat)]),
                 as.matrix(log1p(mat) / log(2)), tolerance = 1e-12)
})


test_that("logNormParam with offset != 1 errors clearly", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))
    expect_error(
        GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 0.5)),
        "offset != 1"
    )
})


test_that("libraryNormParam re-run replaces scale factors, keeps the log op", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 6)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("log", base = 2, offset = 1))
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 5e3))

    # Second library norm replaces the first rather than stacking a second
    # scaling, and leaves the log op in place after it.
    expect_equal(vapply(pe@post_ops, function(o) o$type, character(1L)),
                 c("norm_libsize", "log"))
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(pe@post_ops[[1]]$scalef$scalef, 5e3 / libsz)
})


test_that("normalizeGiotto end-to-end on parquet backend builds the op chain", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")

    mat <- .tiny_mat(seed = 4)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    locs <- data.frame(cell_ID = colnames(mat),
                        sdimx = runif(ncol(mat)),
                        sdimy = runif(ncol(mat)))
    g <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                     verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    g <- suppressWarnings(Giotto::normalizeGiotto(g,
        scalefactor = 1e4,
        scale_feats = FALSE, scale_cells = FALSE,
        verbose = FALSE))

    norm_eo <- GiottoClass::getExpression(g, values = "normalized",
                                            output = "exprObj")
    pe_norm <- slot(norm_eo, "exprMat")
    expect_s4_class(pe_norm, "parquetExprStore")
    expect_equal(vapply(pe_norm@post_ops, function(o) o$type, character(1L)),
                 c("norm_libsize", "log"))

    # Direct math matches stored scale_factors
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(pe_norm@post_ops[[1]]$scalef$scalef, 1e4 / libsz)
})


test_that("cell subset on parquetExprStore slices scalef table", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 7)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    expect_equal(nrow(pe@post_ops[[1]]$scalef), ncol(mat))

    pe2 <- pe[, 1:5]
    expect_equal(nrow(pe2@post_ops[[1]]$scalef), 5L)
    # source_id stays constant for single store
    expect_setequal(unique(pe2@post_ops[[1]]$scalef$source_id), pe@uid)
    expect_setequal(pe2@post_ops[[1]]$scalef$orig_row_id, 1:5)

    # Normalized value for surviving cells must match pre-subset values
    # (value is mutated in place by @post_ops during tibble collect).
    v_full <- data.table::as.data.table(storeRead(pe,  output = "tibble"))
    v_sub  <- data.table::as.data.table(storeRead(pe2, output = "tibble"))
    common <- v_full[row_id %in% 1:5]
    data.table::setorder(common, row_id, col_id)
    data.table::setorder(v_sub, row_id, col_id)
    expect_equal(common$value, v_sub$value)
})


test_that("gene subset on parquetExprStore preserves scalef table", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 8)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    before <- nrow(pe@post_ops[[1]]$scalef)

    pe2 <- pe[1:5, ]  # 5 of the n_genes
    expect_equal(nrow(pe2@post_ops[[1]]$scalef), before)
})


test_that("normalizeGiotto with scale_feats=TRUE errors on parquet backend", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 5)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    locs <- data.frame(cell_ID = colnames(mat),
                        sdimx = runif(ncol(mat)),
                        sdimy = runif(ncol(mat)))
    g <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                     verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    expect_error(
        suppressWarnings(Giotto::normalizeGiotto(g,
            scale_feats = TRUE, scale_cells = TRUE, verbose = FALSE)),
        "z-score scaling"
    )
})


# @post_ops is a chain applied in order, so the two norm ops need not be
# adjacent, need not both be present, and need not appear in a fixed order --
# the chain itself supplies the sequencing. Both executors must honour that:
# the R-side one inside storeRead's materializing modes, and the arrow
# lowering used by the pushed-down stats aggregate.
#
# The reversed case is reachable through the public verbs (log then library
# appends the scaling after the log) and is genuinely distinguishing, since
# log1p(x)/log(2) * s != log1p(x * s)/log(2).
test_that("norm and log ops compose in any order on both executors", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 4)
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    sf <- 1e4 / libsz
    dm <- as.matrix(mat)

    mk <- function() storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    LN <- function(x) GiottoClass::processData(x,
        Giotto::normParam("library", scalefactor = 1e4))
    LG <- function(x) GiottoClass::processData(x,
        Giotto::normParam("log", base = 2, offset = 1))

    cases <- list(
        list(pe = LG(mk()),         types = "log",
             ref = log1p(dm) / log(2)),
        list(pe = LN(mk()),         types = "norm_libsize",
             ref = t(t(dm) * sf)),
        list(pe = LG(LN(mk())),     types = c("norm_libsize", "log"),
             ref = log1p(t(t(dm) * sf)) / log(2)),
        list(pe = LN(LG(mk())),     types = c("log", "norm_libsize"),
             ref = t(t(log1p(dm) / log(2)) * sf)),
        list(pe = LG(LN(LG(mk()))), types = c("log", "norm_libsize", "log"),
             ref = log1p(t(t(log1p(dm) / log(2)) * sf)) / log(2))
    )

    for (cs in cases) {
        expect_equal(
            vapply(cs$pe@post_ops, function(o) o$type, character(1L)),
            cs$types)

        # R-side executor, via storeRead materialization
        got <- as.matrix(storeRead(cs$pe, output = "dgcmatrix")[
            rownames(mat), colnames(mat)])
        expect_equal(got, cs$ref, tolerance = 1e-12, ignore_attr = TRUE)

        # arrow lowering, via the pushed-down per-gene aggregate
        st <- suppressWarnings(GiottoDisk:::.stream_gene_stats(cs$pe))
        expect_equal(st$mean_expr[match(rownames(mat), st$feats)],
                     unname(rowMeans(cs$ref)), tolerance = 1e-12)
    }
})
