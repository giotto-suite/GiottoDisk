# Streaming PAGE enrichment, and the two methods that are refused.
#
# The bar is parity with the in-memory implementation: the same cells in the
# same order with the same scores. Scores agree to floating point rather than
# bit-for-bit, because the per-cell sd is recovered from raw moments and sums
# in a different order than stats::sd. Everything downstream divides by that
# sd, so the tolerance is stated rather than left to expect_equal's default.

# `envir` is the caller's frame: withr's local_* would otherwise unwind when
# this helper returns, deleting the project directory before a single test has
# read from it.
.enrich_pe_fixture <- function(G = 300L, N = 200L, seed = 21L,
                               envir = parent.frame()) {
    withr::local_options(giotto.check_valid = FALSE, giotto.verbose = FALSE,
                         giotto.no_python_warn = TRUE, .local_envir = envir)
    set.seed(seed)
    m <- matrix(rpois(G * N, 2), G, N,
        dimnames = list(paste0("g", seq_len(G)), paste0("c", seq_len(N))))
    locs <- data.table::data.table(
        cell_ID = colnames(m),
        sdimx = runif(N, 0, 10), sdimy = runif(N, 0, 10)
    )
    dir <- file.path(withr::local_tempdir(.local_envir = envir), "proj")
    list(
        m = m,
        backed = GiottoClass::createGiottoObject(
            expression = m, spatial_locs = locs, backend = dir),
        mem = GiottoClass::createGiottoObject(
            expression = m, spatial_locs = locs)
    )
}

.enrich_sign_matrix <- function(genes, sizes = c(60L, 40L, 80L, 3L),
                                seed = 3L) {
    set.seed(seed)
    types <- paste0("t", LETTERS[seq_along(sizes)])
    sm <- matrix(0L, length(genes), length(sizes),
                 dimnames = list(genes, types))
    for (j in seq_along(sizes)) sm[sample(length(genes), sizes[j]), j] <- 1L
    sm
}

# the accumulators the streaming statistic is built out of ------------------

test_that("the expm1 op inverts log and stays bit-exact", {
    skip_if_not_installed("Giotto")
    f <- .enrich_pe_fixture(G = 40L, N = 60L)
    pe <- GiottoClass::getExpression(f$backed, values = "raw",
                                     output = "exprObj")[]

    for (base in c(2, exp(1), 10)) {
        src <- .pe_push_op(pe, list(type = "expm1", base = base))
        got <- .stream_expr_accum(src, axis = "feat", stats = "sum")$sum
        # bit-exact, not merely close: `m` feeds a z-score's denominator, and
        # the op is written as base^v - 1 for exactly this reason
        expect_identical(unname(got), unname(rowSums(base^f$m - 1)),
                         info = paste("base", base))
    }

    # the op is also honoured on the R-side (post-op) path
    dt <- data.table::data.table(value = c(0, 1, 2.5))
    expect_equal(.op_transform_expm1(dt, list(base = 2))$value,
                 2^c(0, 1, 2.5) - 1)
})

test_that("a per-feature multiply payload reaches the right rows", {
    skip_if_not_installed("Giotto")
    f <- .enrich_pe_fixture(G = 40L, N = 60L)
    pe <- GiottoClass::getExpression(f$backed, values = "raw",
                                     output = "exprObj")[]

    set.seed(7)
    w <- runif(nrow(f$m), 0.5, 2)
    weighted <- .pe_push_op(pe, list(
        type = "multiply", axis = "feat",
        factors = .pe_feat_factor_payload(pe, w)
    ))
    got <- .stream_expr_accum(weighted, axis = "cell", stats = "sum")$sum
    expect_equal(unname(got), unname(colSums(f$m * w)))
})


# parity -------------------------------------------------------------------

test_that("streaming PAGE matches the in-memory method", {
    skip_if_not_installed("Giotto")
    f <- .enrich_pe_fixture()
    sm <- .enrich_sign_matrix(rownames(f$m))
    pe <- GiottoClass::getExpression(f$backed, values = "raw",
                                     output = "exprObj")[]
    expect_s4_class(pe, "parquetExprStore")
    mm <- as.matrix(GiottoClass::getExpression(f$mem, values = "raw",
                                               output = "matrix"))

    configs <- list(
        default = list(),
        zscore = list(output_enrichment = "zscore"),
        no_reverse = list(reverse_log_scale = FALSE),
        natural_log = list(logbase = exp(1)),
        strict_overlap = list(min_overlap_genes = 50)
    )
    for (nm in names(configs)) {
        p <- do.call(Giotto::enrichParam,
                     c(list("PAGE", verbose = FALSE), configs[[nm]]))
        streamed <- GiottoClass::analyzeData(pe, p, sign_matrix = sm)
        memory <- data.table::as.data.table(
            GiottoClass::analyzeData(mm, p, sign_matrix = sm)
        )
        data.table::setcolorder(memory, names(streamed))

        # same cells, same order, same cell types kept
        expect_identical(streamed$cell_ID, memory$cell_ID, info = nm)
        expect_identical(names(streamed), names(memory), info = nm)
        for (cl in setdiff(names(streamed), "cell_ID")) {
            expect_equal(streamed[[cl]], memory[[cl]],
                         tolerance = 1e-10, info = paste(nm, cl))
        }
    }

    # tD has 3 markers and is dropped at the default threshold of 5; tB has 40
    # and survives it but not min_overlap_genes = 50. Both assertions above
    # would pass vacuously if the drop never happened.
    d <- GiottoClass::analyzeData(pe, Giotto::enrichParam("PAGE",
        verbose = FALSE), sign_matrix = sm)
    expect_false("tD" %in% names(d))
    expect_true("tB" %in% names(d))
    s <- GiottoClass::analyzeData(pe, Giotto::enrichParam("PAGE",
        min_overlap_genes = 50, verbose = FALSE), sign_matrix = sm)
    expect_false("tB" %in% names(s))
})

test_that("runPAGEEnrich runs end to end on a backed object", {
    skip_if_not_installed("Giotto")
    f <- .enrich_pe_fixture()
    sm <- .enrich_sign_matrix(rownames(f$m))

    # This is the call that used to fail outright: the wrapper densified the
    # store before handing it on.
    backed <- Giotto::runPAGEEnrich(f$backed, sign_matrix = sm,
        expression_values = "raw", return_gobject = FALSE,
        verbose = FALSE)$matrix[]
    memory <- Giotto::runPAGEEnrich(f$mem, sign_matrix = sm,
        expression_values = "raw", return_gobject = FALSE,
        verbose = FALSE)$matrix[]
    backed <- data.table::as.data.table(backed)
    memory <- data.table::as.data.table(memory)
    data.table::setcolorder(memory, names(backed))

    expect_identical(backed$cell_ID, memory$cell_ID)
    for (cl in setdiff(names(backed), "cell_ID")) {
        expect_equal(backed[[cl]], memory[[cl]], tolerance = 1e-10, info = cl)
    }

    g2 <- Giotto::runPAGEEnrich(f$backed, sign_matrix = sm,
        expression_values = "raw", verbose = FALSE)
    expect_true("PAGE" %in% GiottoClass::list_spatial_enrichments_names(
        g2, spat_unit = "cell", feat_type = "rna"))
})


# what is refused ----------------------------------------------------------

test_that("the unsupported enrichment paths say why", {
    skip_if_not_installed("Giotto")
    f <- .enrich_pe_fixture(G = 150L, N = 60L)
    sm <- .enrich_sign_matrix(rownames(f$m), sizes = c(30L, 30L))
    pe <- GiottoClass::getExpression(f$backed, values = "raw",
                                     output = "exprObj")[]

    expect_error(
        GiottoClass::analyzeData(pe,
            Giotto::enrichParam("PAGE", p_value = TRUE), sign_matrix = sm),
        "p_value = TRUE` is not supported"
    )
    expect_error(
        GiottoClass::analyzeData(pe, Giotto::enrichParam("rank"),
                                 sign_matrix = sm),
        "ranks each gene across every cell"
    )
    expect_error(
        GiottoClass::analyzeData(pe,
            Giotto::enrichParam("hypergeometric"), sign_matrix = sm),
        "not supported on a disk-backed"
    )
})
