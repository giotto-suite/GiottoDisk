# Streaming PAGE and hypergeometric enrichment, and rank, which is refused.
#
# The bar is parity with the in-memory implementation: the same cells in the
# same order with the same scores. Scores agree to floating point rather than
# bit-for-bit, because the reference level is summed window by window, so the
# tolerance is stated rather than left to expect_equal's default.

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
        factors = .pe_axis_factor_payload(pe, w, "feat")
    ))
    got <- .stream_expr_accum(weighted, axis = "cell", stats = "sum")$sum
    expect_equal(unname(got), unname(colSums(f$m * w)))
})


# parity -------------------------------------------------------------------

test_that("streaming enrichment matches the in-memory method", {
    skip_if_not_installed("Giotto")
    # several windows, so cells scored apart must agree with all at once
    withr::local_options(giottodisk.chunk_size = 37L)
    f <- .enrich_pe_fixture()
    sm <- .enrich_sign_matrix(rownames(f$m))
    pe <- GiottoClass::getExpression(f$backed, values = "raw",
                                     output = "exprObj")[]
    expect_s4_class(pe, "parquetExprStore")
    mm <- as.matrix(GiottoClass::getExpression(f$mem, values = "raw",
                                               output = "matrix"))

    page <- function(...) Giotto::enrichParam("PAGE", verbose = FALSE, ...)
    hyper <- function(...) Giotto::enrichParam("hypergeometric", ...)
    configs <- list(
        default = page(),
        zscore = page(output_enrichment = "zscore"),
        p_value = page(p_value = TRUE),
        depletion = page(p_value = TRUE, include_depletion = TRUE),
        no_reverse = page(reverse_log_scale = FALSE),
        natural_log = page(logbase = exp(1)),
        strict_overlap = page(min_overlap_genes = 50),
        hyper = hyper(),
        hyper_zscore = hyper(output_enrichment = "zscore"),
        hyper_p = hyper(p_value = TRUE),
        hyper_top = hyper(top_percentage = 20, reverse_log_scale = FALSE)
    )
    for (nm in names(configs)) {
        p <- configs[[nm]]
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

test_that("histogram ranks stay within one bin of the exact ranks", {
    set.seed(4)
    G <- 20L; N <- 500L
    M <- Matrix::rsparsematrix(G, N, density = 0.4,
        rand.x = function(n) rexp(n))
    dimnames(M) <- list(paste0("g", seq_len(G)), paste0("c", seq_len(N)))
    n_bins <- 64L
    # the histogram .pe_value_hist() returns, built directly
    mx <- apply(as.matrix(M), 1L, max)
    s <- t(t(as.matrix(M)) * 0) + as.matrix(M) * n_bins / mx
    n <- t(apply(s, 1L, function(v) {
        tabulate(pmin(ceiling(v[v > 0]) - 1, n_bins - 1) + 1L, n_bins)
    }))
    h <- list(n = n, below = cbind(0, t(apply(n, 1L, cumsum)))[, 1:n_bins],
              max = mx, n_zero = N - rowSums(n), n_bins = n_bins)

    approx <- .rank_from_hist(M, h, seq_len(G), "average")
    exact <- t(apply(as.matrix(M), 1L, rank, ties.method = "average"))
    for (g in seq_len(G)) {
        v <- M[g, ]
        on <- v > 0
        b <- pmin(ceiling(v[on] * n_bins / mx[g]) - 1, n_bins - 1) + 1L
        expect_true(all(abs(approx[g, on] - exact[g, on]) <= n[g, b] + 1))
        expect_true(all(approx[g, !on] == exact[g, !on]))
        # order within a feature is kept
        expect_equal(cor(approx[g, on], exact[g, on], method = "spearman"),
                     1, tolerance = 1e-12)
    }
})

test_that("streaming rank enrichment is close to the in-memory method", {
    skip_if_not_installed("Giotto")
    withr::local_options(giottodisk.chunk_size = 300L)
    withr::local_options(giotto.check_valid = FALSE, giotto.verbose = FALSE,
                         giotto.no_python_warn = TRUE)
    set.seed(8)
    G <- 150L; N <- 1000L
    m <- matrix(log2(rpois(G * N, 3) * rexp(G * N) + 1), G, N,
        dimnames = list(paste0("g", seq_len(G)), paste0("c", seq_len(N))))
    dir <- file.path(withr::local_tempdir(), "proj")
    backed <- GiottoClass::createGiottoObject(expression = m,
        spatial_locs = data.table::data.table(cell_ID = colnames(m),
            sdimx = runif(N), sdimy = runif(N)),
        backend = dir)
    pe <- GiottoClass::getExpression(backed, values = "raw",
                                     output = "exprObj")[]
    rm_ <- sapply(1:3, function(j) sample(G))
    dimnames(rm_) <- list(rownames(m), paste0("t", 1:3))

    # "max" ranks are whole numbers, so on 1,000 cells two features of a cell
    # often share one exactly and the within-cell rule settles the tie; any
    # nudge, even +-0.5 on the exact ranks, moves those cells. Hence the
    # looser bar there. Both agree to > 0.9999 on a million cells.
    bar <- c(average = 0.999, max = 0.995)
    for (ties in c("average", "max")) {
        for (out in c("original", "zscore")) {
            p <- Giotto::enrichParam("rank", ties_method = ties,
                                     output_enrichment = out)
            streamed <- GiottoClass::analyzeData(pe, p, sign_matrix = rm_)
            memory <- GiottoClass::analyzeData(m, p, sign_matrix = rm_)
            expect_identical(streamed$cell_ID, memory$cell_ID)
            for (cl in colnames(rm_)) {
                expect_gt(cor(streamed[[cl]], memory[[cl]],
                              method = "spearman"), bar[[ties]])
            }
        }
    }

    expect_error(
        GiottoClass::analyzeData(pe, Giotto::enrichParam("rank",
            p_value = TRUE), sign_matrix = rm_),
        "not supported on a disk-backed store"
    )
})


# --- the chain's phase split -------------------------------------------------
#
# ADR 0002: once `@post_ops` is non-empty, subsequent pushes go there regardless
# of the record's natural phase. `.pe_push_op()` enforces that by erroring on a
# lazy push, so a consumer that pushes `"lazy"` unconditionally fails on any
# store carrying a post op -- which is what `.stream_page()` used to do. No
# other test builds such a store, which is why it survived.

test_that("streaming PAGE works whichever phase the chain is already in", {
    skip_if_not_installed("Giotto")
    f <- .enrich_pe_fixture(G = 300L, N = 200L)
    sm <- .enrich_sign_matrix(rownames(f$m), sizes = c(60L, 40L, 80L))
    pe <- GiottoClass::getExpression(f$backed, values = "raw",
                                     output = "exprObj")[]
    p <- Giotto::enrichParam("PAGE", verbose = FALSE)

    # The same record, placed either side of the split. The chain computes the
    # same thing either way, so the scores must agree.
    #
    # Agree to ~1e-17, not bit-for-bit: on the lazy side Arrow's log2 kernel
    # evaluates it, on the post side R's log1p does. Same function, different
    # implementations, last-ulp disagreement. Asserting identity here would be
    # asserting that two C++ and R math libraries round the same way.
    pe_lazy <- .pe_push_op(pe, list(type = "log", base = 2), phase = "lazy")
    pe_post <- .pe_push_op(pe, list(type = "log", base = 2), phase = "post")
    expect_length(pe_lazy@post_ops, 0L)
    expect_length(pe_post@post_ops, 1L)

    a <- GiottoClass::analyzeData(pe_lazy, p, sign_matrix = sm)
    b <- GiottoClass::analyzeData(pe_post, p, sign_matrix = sm)

    expect_identical(a$cell_ID, b$cell_ID)
    expect_identical(names(a), names(b))
    for (cl in setdiff(names(a), "cell_ID")) {
        expect_equal(a[[cl]], b[[cl]], tolerance = 1e-12, info = cl)
        expect_lt(max(abs(a[[cl]] - b[[cl]])), 1e-14)
    }
})

test_that(".pe_push_op_last appends at the end of the chain", {
    skip_if_not_installed("Giotto")
    f <- .enrich_pe_fixture(G = 40L, N = 60L)
    pe <- GiottoClass::getExpression(f$backed, values = "raw",
                                     output = "exprObj")[]
    op <- list(type = "expm1", base = 2)

    # empty chain -> lazy, so it can still lower to Acero
    lazy <- .pe_push_op_last(pe, op)
    expect_length(lazy@ops, 1L)
    expect_length(lazy@post_ops, 0L)

    # a post op already queued -> the new record has to run after it
    pe2 <- .pe_push_op(pe, list(type = "log", base = 2), phase = "post")
    post <- .pe_push_op_last(pe2, op)
    expect_length(post@ops, 0L)
    expect_identical(vapply(post@post_ops, function(o) o$type, ""),
                     c("log", "expm1"))
})
