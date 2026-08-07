# Tests for streaming normalize dispatch:
#   processData(parquetExprStore, libraryNormParam) -> appends multiply
#       op (log = FALSE) to pe@ops
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


test_that("libraryNormParam appends a multiply op to @post_ops", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    pe2 <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))

    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_length(pe2@ops, 1L)
    expect_equal(pe2@ops[[1]]$type, "multiply")
    expect_equal(pe2@ops[[1]]$axis, "cell")
    expect_equal(unname(pe2@ops[[1]]$factors[[pe@uid]]), unname(1e4 / libsz))
    # scaling carries no log state -- that is a separate op
    expect_null(pe2@ops[[1]]$log)
})


test_that("logNormParam appends an independent log op", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 2)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("library", scalefactor = 1e4))
    pe  <- GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 1))

    # Two records in chain order, not one fused record.
    expect_length(pe@ops, 2L)
    expect_equal(vapply(pe@ops, function(o) o$type, character(1L)),
                 c("multiply", "log"))
    expect_equal(pe@ops[[2]]$base, 2)

    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(unname(pe@ops[[1]]$factors[[pe@uid]]),
                 unname(1e4 / libsz))
})


# log-only on raw counts used to be rejected outright: the log flag lived on
# the libsize record, so there was nowhere to put it without one. As its own
# op it stands alone.
test_that("logNormParam works with no library normalization", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 3)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe, Giotto::normParam("log", base = 2, offset = 1))

    expect_length(pe@ops, 1L)
    expect_equal(pe@ops[[1]]$type, "log")

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


test_that("libraryNormParam re-run appends rather than editing its own record", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 6)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("log", base = 2, offset = 1))
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 5e3))

    # Records are positional: the second library norm appends at the end of
    # the chain rather than reaching back to rewrite the first, which would
    # be rewriting it *underneath* the intervening log.
    expect_equal(vapply(pe@ops, function(o) o$type, character(1L)),
                 c("multiply", "log", "multiply"))

    # The first still holds its original factors ...
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(unname(pe@ops[[1]]$factors[[pe@uid]]),
                 unname(1e4 / libsz))

    # ... and the last normalizes whatever the chain produced before it, so
    # the result's columns land on the new scalefactor.
    got <- as.matrix(storeRead(pe, output = "dgcmatrix"))
    expect_equal(unname(colSums(got)), rep(5e3, ncol(mat)), tolerance = 1e-8)
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
    expect_equal(vapply(pe_norm@ops, function(o) o$type, character(1L)),
                 c("multiply", "log"))

    # Direct math matches stored scale_factors
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    expect_equal(unname(pe_norm@ops[[1]]$factors[[pe_norm@uid]]),
                 unname(1e4 / libsz))
})


test_that("cell subset on parquetExprStore slices scalef table", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 7)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    expect_length(pe@ops[[1]]$factors[[pe@uid]], ncol(mat))

    pe2 <- pe[, 1:5]
    # The payload is a vector indexed by ON-DISK row_id, so narrowing the view
    # neither invalidates nor misaligns it -- there is nothing to slice.
    expect_identical(pe2@ops[[1]]$factors, pe@ops[[1]]$factors)

    # Normalized value for surviving cells must match pre-subset values
    # (value is mutated in place by @post_ops during tibble collect).
    v_full <- data.table::as.data.table(storeRead(pe,  output = "tibble"))
    v_sub  <- data.table::as.data.table(storeRead(pe2, output = "tibble"))
    common <- v_full[row_id %in% 1:5]
    data.table::setorder(common, row_id, col_id)
    data.table::setorder(v_sub, row_id, col_id)
    expect_equal(common$value, v_sub$value)
})


test_that("gene subset on parquetExprStore preserves the factor payload", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 8)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe  <- GiottoClass::processData(pe,
        Giotto::normParam("library", scalefactor = 1e4))
    before <- pe@ops[[1]]$factors

    pe2 <- pe[1:5, ]  # 5 of the n_genes
    # A cell-axis payload is untouched by a gene subset -- which is what makes
    # normalize-then-subset keep full-library scaling.
    expect_identical(pe2@ops[[1]]$factors, before)
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
    dm  <- as.matrix(mat)

    # Scale factors are taken from the values as the chain leaves them at the
    # norm's own position -- not from the raw counts. So a log ahead of the
    # norm changes them, and each case needs its own.
    sfac <- function(x) {
        cs <- colSums(x); cs[cs == 0] <- 1; 1e4 / cs
    }
    lg2 <- function(x) log1p(x) / log(2)

    mk <- function() storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    LN <- function(x) GiottoClass::processData(x,
        Giotto::normParam("library", scalefactor = 1e4))
    LG <- function(x) GiottoClass::processData(x,
        Giotto::normParam("log", base = 2, offset = 1))

    cases <- list(
        list(pe = LG(mk()),         types = "log",
             ref = log1p(dm) / log(2)),
        list(pe = LN(mk()),         types = "multiply",
             ref = t(t(dm) * sfac(dm))),
        list(pe = LG(LN(mk())),     types = c("multiply", "log"),
             ref = lg2(t(t(dm) * sfac(dm)))),
        # log first: the norm scales the LOGGED values, so its factors come
        # from colSums of those, not of the counts.
        list(pe = LN(LG(mk())),     types = c("log", "multiply"),
             ref = t(t(lg2(dm)) * sfac(lg2(dm)))),
        list(pe = LG(LN(LG(mk()))), types = c("log", "multiply", "log"),
             ref = lg2(t(t(lg2(dm)) * sfac(lg2(dm)))))
    )

    for (cs in cases) {
        expect_equal(
            vapply(cs$pe@ops, function(o) o$type, character(1L)),
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


# The `multiply` op has two executors that cannot share code: arrow cannot
# reach an R vector from inside a plan, so the lowering joins a factor table
# while the R-side path indexes the vector directly. That index is ~4x faster
# than an in-R join and is why the second implementation exists -- but two
# implementations of one op drift. They have already disagreed once, on rows
# the payload does not cover (index yields 0 or NA by position; the join
# yields NA for every miss).
#
# So assert they agree with EACH OTHER on the same store, not merely that each
# agrees with a reference separately.
test_that("multiply: arrow lowering and R-side executor agree", {
    skip_if_not_installed("Giotto")

    # Same record, moved between phases: @ops takes the arrow lowering,
    # @post_ops takes the R-side executor.
    as_lazy <- function(p) { p@ops <- c(p@ops, p@post_ops); p@post_ops <- list(); p }
    dense   <- function(p) as.matrix(storeRead(p, output = "dgcmatrix"))

    mk_store <- function(seed, n_cells, prefix) {
        set.seed(seed)
        m <- Matrix::rsparsematrix(9, n_cells, density = 0.7,
            rand.x = function(n) as.double(rpois(n, 5L) + 1L))
        rownames(m) <- paste0("g", seq_len(9))
        colnames(m) <- paste0(prefix, seq_len(n_cells))
        storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    }
    norm <- function(p) GiottoClass::processData(p,
        Giotto::normParam("library", scalefactor = 6e3))

    pe <- mk_store(21L, 10L, "a")
    pe2 <- mk_store(22L, 7L, "b")

    cases <- list(
        "single store"          = norm(pe),
        "gene subset"           = norm(pe)[c(2L, 5L, 8L), ],
        "cell subset"           = norm(pe)[, c(1L, 4L, 9L)],
        "both axes"             = norm(pe)[c(1L, 3L), c(2L, 5L, 10L)],
        # multi-source: the R executor takes its split-by-source branch and
        # the lowering takes a composite (source_id, id) join
        "union"                 = norm(unionParquetExprStore(list(pe, pe2))),
        "union cross-substore"  = norm(unionParquetExprStore(list(pe, pe2)))[
                                      , c(2L, 6L, 12L, 16L)]
    )

    for (nm in names(cases)) {
        p <- cases[[nm]]
        expect_equal(dense(p), dense(as_lazy(p)),
                     tolerance = 1e-12, info = nm)
    }

    # feat axis and the scalar case are only reachable by hand today, but both
    # executors claim to serve them
    pf <- pe
    pf@post_ops <- list(list(type = "multiply", axis = "feat",
        factors = stats::setNames(list(seq_len(9) / 3), pe@uid)))
    expect_equal(dense(pf), dense(as_lazy(pf)), tolerance = 1e-12)

    pa <- pe
    pa@post_ops <- list(list(type = "multiply", axis = "all", factors = 2.5))
    expect_equal(dense(pa), dense(as_lazy(pa)), tolerance = 1e-12)

    # Rows the payload does not cover are where these two most easily diverge:
    # the join yields NA for any miss, while a positional index yields whatever
    # sits at that slot. The NA-filled payload makes both give NA -- the older
    # zero-filled table gave 0 on the R side. Nothing produces an uncovered row
    # today, so this only holds while someone keeps checking it.
    pg <- pe
    pg@post_ops <- list(list(type = "multiply", axis = "cell",
        factors = stats::setNames(
            list(c(2, 2, NA, 2, 2, 2)), pe@uid)))   # gap at 3, short of 7:10
    got_r <- dense(pg)
    got_a <- dense(as_lazy(pg))
    expect_equal(got_r, got_a, tolerance = 1e-12)
    # Only STORED entries become NA -- structural zeros in an uncovered column
    # stay 0, since a sparse matrix has no slot for them to be NA in.
    expect_true(anyNA(got_r[, 7:10]))                  # uncovered ids
    expect_false(anyNA(got_r[, c(1L, 2L, 4L, 5L)]))    # covered ids
})
