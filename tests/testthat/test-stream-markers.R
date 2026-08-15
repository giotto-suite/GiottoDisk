# Tests for streaming pairwise marker detection. The statistic pass
# (`.pe_group_moments`) is checked against a dense reference computed with
# base/Matrix, and the pairwise tail against `stats::t.test` -- independent
# ground truth in both cases, not another run of the same path.

.mk_mat <- function(n_genes = 15L, n_cells = 40L, density = 0.4, seed = 11L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}

# Dense per-group moments. Deliberately naive: loops the groups and uses
# base rowMeans / apply(var) on the densified block, so it shares no code with
# the streaming path.
.ref_moments <- function(mat, groups) {
    d <- as.matrix(mat)
    g <- droplevels(factor(groups))
    lvls <- levels(g)
    n <- setNames(numeric(length(lvls)), lvls)
    means <- vars <- matrix(0, nrow(d), length(lvls),
        dimnames = list(rownames(d), lvls))
    for (k in lvls) {
        cols <- which(!is.na(g) & g == k)
        n[[k]] <- length(cols)
        blk <- d[, cols, drop = FALSE]
        means[, k] <- rowMeans(blk)
        vars[, k] <- if (length(cols) > 1L) {
            apply(blk, 1L, stats::var)
        } else {
            0
        }
    }
    list(n = n, means = means, vars = vars)
}


test_that(".pe_group_moments matches a dense reference on a single store", {
    mat <- .mk_mat()
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b", "c"), length.out = ncol(mat))

    got <- GiottoDisk:::.pe_moments_derive(
        GiottoDisk:::.pe_group_moments(pe, groups))
    ref <- .ref_moments(mat, groups)

    expect_equal(got$n, ref$n)
    expect_equal(got$means, ref$means)
    expect_equal(got$vars, ref$vars)
})


test_that(".pe_group_moments counts all cells, not just stored entries", {
    # A gene that is entirely absent from one group must come back as mean 0
    # and variance 0 over that group's FULL cell count, not as a missing row.
    mat <- .mk_mat(n_genes = 6L, n_cells = 20L, density = 0.5, seed = 3L)
    mat[2L, 1:10] <- 0
    mat <- Matrix::drop0(mat)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b"), each = 10L)

    got <- GiottoDisk:::.pe_moments_derive(
        GiottoDisk:::.pe_group_moments(pe, groups))
    expect_equal(got$n, c(a = 10, b = 10))
    expect_equal(got$means[2L, "a"], 0)
    expect_equal(got$vars[2L, "a"], 0)
    expect_equal(got$means, .ref_moments(mat, groups)$means)
    expect_equal(got$vars, .ref_moments(mat, groups)$vars)
})


test_that(".pe_group_moments excludes NA-group cells from n and the sums", {
    mat <- .mk_mat(n_genes = 8L, n_cells = 24L, seed = 5L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b", NA), length.out = ncol(mat))

    got <- GiottoDisk:::.pe_moments_derive(
        GiottoDisk:::.pe_group_moments(pe, groups))
    ref <- .ref_moments(mat, groups)

    expect_equal(colnames(got$means), c("a", "b"))
    expect_equal(got$n, ref$n)
    expect_equal(got$means, ref$means)
    expect_equal(got$vars, ref$vars)
})


test_that(".pe_group_moments follows a `[` view rather than the whole store", {
    mat <- .mk_mat(n_genes = 10L, n_cells = 30L, seed = 9L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    gsel <- c(3L, 1L, 8L)      # out of order on purpose -- `[` preserves it
    csel <- 5:24
    view <- pe[gsel, csel]
    sub <- mat[gsel, csel, drop = FALSE]
    groups <- rep(c("x", "y"), length.out = length(csel))

    got <- GiottoDisk:::.pe_moments_derive(
        GiottoDisk:::.pe_group_moments(view, groups))
    ref <- .ref_moments(sub, groups)

    expect_equal(rownames(got$means), rownames(sub))
    expect_equal(got$means, ref$means)
    expect_equal(got$vars, ref$vars)
})


test_that(".pe_group_moments agrees across the Acero and chunked paths", {
    # A non-empty @post_ops chain forces the chunked R-side branch. Both
    # branches must produce the same moments for the same values.
    mat <- .mk_mat(n_genes = 12L, n_cells = 36L, seed = 13L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b", "c", "d"), length.out = ncol(mat))

    lazy <- GiottoDisk:::.pe_group_moments(pe, groups)
    expect_length(pe@post_ops, 0L)

    demoted <- GiottoDisk:::.pe_demote_ops(
        GiottoDisk:::.pe_push_op(
            pe, list(type = "multiply", axis = "all", factors = 1),
            phase = "lazy"
        ), 1L
    )
    expect_true(length(demoted@post_ops) > 0L)
    chunked <- GiottoDisk:::.pe_group_moments(demoted, groups)

    expect_equal(chunked$n, lazy$n)
    expect_equal(chunked$means, lazy$means)
    expect_equal(chunked$vars, lazy$vars)
})


test_that(".pe_group_moments validates groups length and level count", {
    mat <- .mk_mat(n_genes = 5L, n_cells = 10L, seed = 2L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    expect_error(
        GiottoDisk:::.pe_group_moments(pe, rep("a", 3L)),
        "one entry per cell"
    )
    expect_error(
        GiottoDisk:::.pe_group_moments(pe, rep("a", 10L)),
        "at least 2 non-empty groups"
    )
})


test_that("pairwise Welch statistics match stats::t.test", {
    # Independent ground truth for the transcribed tail: t.test() on the raw
    # per-group values must reproduce the moment-based statistic and its
    # two-sided p-value exactly.
    mat <- .mk_mat(n_genes = 7L, n_cells = 30L, density = 0.7, seed = 21L)
    groups <- rep(c("a", "b", "c"), each = 10L)
    mom <- .ref_moments(mat, groups)
    d <- as.matrix(mat)

    tt <- GiottoDisk:::.pe_welch(
        host_s2 = mom$vars[, "a"], target_s2 = mom$vars[, "b"],
        host_n = mom$n[["a"]], target_n = mom$n[["b"]]
    )
    lfc <- mom$means[, "a"] - mom$means[, "b"]
    p_out <- GiottoDisk:::.pe_run_t(lfc, tt$err, tt$df)
    log_p <- GiottoDisk:::.pe_choose_lr(p_out$left, p_out$right, "any")

    for (i in seq_len(nrow(d))) {
        ref <- stats::t.test(d[i, groups == "a"], d[i, groups == "b"],
            var.equal = FALSE)
        expect_equal(unname(tt$df[i]), unname(ref$parameter), tolerance = 1e-10)
        expect_equal(unname(lfc[i] / sqrt(tt$err[i])),
            unname(ref$statistic), tolerance = 1e-10)
        expect_equal(exp(log_p[i]), ref$p.value, tolerance = 1e-8)
    }
})


test_that(".pe_logBH matches p.adjust(method = 'BH') in log space", {
    set.seed(4)
    p <- runif(50)^3
    expect_equal(
        exp(GiottoDisk:::.pe_logBH(log(p))),
        stats::p.adjust(p, method = "BH"),
        tolerance = 1e-12
    )
})


test_that("analyzeData(scranMarkersParam) refuses tests it cannot stream", {
    skip_if_not_installed("Giotto")
    mat <- .mk_mat(n_genes = 5L, n_cells = 12L, seed = 6L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b"), 6L)

    expect_error(
        GiottoClass::analyzeData(pe, Giotto::markersParam()),
        "`groups` is required"
    )
    expect_error(
        GiottoClass::analyzeData(
            pe, Giotto::markersParam(test_type = "wilcox"), groups = groups),
        "not available on the streaming backend"
    )
})


test_that("analyzeData(scranMarkersParam) returns one table per group", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("scran")
    mat <- .mk_mat(n_genes = 20L, n_cells = 45L, density = 0.6, seed = 31L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b", "c"), each = 15L)

    res <- GiottoClass::analyzeData(pe, Giotto::markersParam(), groups = groups)
    expect_named(as.list(res), c("a", "b", "c"))
    expect_equal(nrow(res[["a"]]), nrow(mat))
    expect_setequal(rownames(res[["a"]]), rownames(mat))
    expect_true(all(c("Top", "p.value", "FDR") %in% colnames(res[["a"]])))
})


test_that("streaming markers match scran::findMarkers on the same values", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("scran")
    mat <- .mk_mat(n_genes = 25L, n_cells = 60L, density = 0.6, seed = 41L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b", "c"), each = 20L)

    got <- GiottoClass::analyzeData(pe, Giotto::markersParam(), groups = groups)
    ref <- scran::findMarkers(as.matrix(mat), groups = groups)

    # Elementwise, not by correlation: a per-feature scale error is exactly
    # what a correlation gate cannot see.
    for (k in c("a", "b", "c")) {
        g <- got[[k]][rownames(mat), ]
        r <- ref[[k]][rownames(mat), ]
        expect_equal(g$p.value, r$p.value, tolerance = 1e-10)
        expect_equal(g$FDR, r$FDR, tolerance = 1e-10)
        # Integer-valued, so this one must match exactly -- a rank that
        # shifted would mean the ordering diverged, not just the last bit.
        expect_equal(g$Top, r$Top)
        # Effect sizes too, not only the p-values: a per-feature scale error
        # would leave p-values intact while moving every logFC.
        expect_equal(g$summary.logFC, r$summary.logFC, tolerance = 1e-10)
        for (other in setdiff(c("a", "b", "c"), k)) {
            nm <- paste0("logFC.", other)
            expect_equal(g[[nm]], r[[nm]], tolerance = 1e-10, info = nm)
        }
    }
})


test_that("grouped featStats `stats` selection drops unrequested accumulators", {
    skip_if_not_installed("Giotto")
    mat <- .mk_mat(n_genes = 10L, n_cells = 30L, seed = 51L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b", "c"), each = 10L)
    p <- Giotto::analyzeParam("feat_stats")

    full <- GiottoClass::analyzeData(pe, p, groups = groups)
    part <- GiottoClass::analyzeData(pe, p, groups = groups,
        stats = c("sum", "sumsq"))

    # nnz / sum_det were never accumulated, so their columns are absent.
    expect_true(all(c("nr_cells", "perc_cells", "mean_expr_det") %in%
        colnames(full)))
    expect_false(any(c("nr_cells", "perc_cells", "mean_expr_det") %in%
        colnames(part)))
    # Everything the requested accumulators support is present and identical.
    expect_true(all(c("total_expr", "mean_expr", "sumsq", "sd") %in%
        colnames(part)))
    for (nm in c("feats", "group", "n_cells", "total_expr", "mean_expr", "sd")) {
        expect_equal(part[[nm]], full[[nm]], info = nm)
    }
})


test_that("`stats` selection is refused on the ungrouped path", {
    skip_if_not_installed("Giotto")
    mat <- .mk_mat(n_genes = 5L, n_cells = 10L, seed = 52L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(
        GiottoClass::analyzeData(pe, Giotto::analyzeParam("feat_stats"),
            stats = c("sum")),
        "only available on the grouped path"
    )
})


test_that(".pe_pool_moments is exact against a single-group pass", {
    mat <- .mk_mat(n_genes = 12L, n_cells = 30L, seed = 53L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    groups <- rep(c("a", "b", "c"), each = 10L)

    split <- GiottoDisk:::.pe_group_moments(pe, groups)
    pooled <- GiottoDisk:::.pe_pool_moments(
        split, list(all = c("a", "b", "c")))
    # Ground truth: two real groups, then compare the pooled "all" against a
    # dense reference over every cell.
    d <- GiottoDisk:::.pe_moments_derive(pooled)
    dm <- as.matrix(mat)

    expect_equal(unname(pooled$n[["all"]]), ncol(mat))
    expect_equal(unname(d$means[, "all"]), unname(rowMeans(dm)))
    expect_equal(unname(d$vars[, "all"]), unname(apply(dm, 1L, stats::var)))
})


test_that("one-vs-rest on disk matches the in-memory pooled scran run", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("scran")
    mat <- .mk_mat(n_genes = 20L, n_cells = 45L, density = 0.6, seed = 61L)
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    lv <- c("a", "b", "c")
    groups <- rep(lv, each = 15L)

    got <- suppressMessages(GiottoClass::analyzeData(
        pe, Giotto::markersParam(comparison = "one_vs_rest"), groups = groups))
    expect_named(as.list(got), lv)

    # Independent ground truth: scran on the materialized matrix, one pooled
    # two-level comparison per cluster -- what the G-pass in-memory path does.
    for (k in lv) {
        rest <- setdiff(lv, k)
        pooled <- ifelse(groups == k, k, paste0(rest, collapse = "_"))
        ref <- scran::findMarkers(as.matrix(mat), groups = pooled)[[k]]
        g <- got[[k]][rownames(mat), ]
        r <- ref[rownames(mat), ]
        expect_equal(g$p.value, r$p.value, tolerance = 1e-10)
        expect_equal(g$FDR, r$FDR, tolerance = 1e-10)
        # The pooled effect size is what pooling the accumulators has to get
        # right; p-values alone would not catch a bad `rest` group.
        nm <- paste0("logFC.", paste0(rest, collapse = "_"))
        expect_equal(g[[nm]], r[[nm]], tolerance = 1e-10, info = nm)
    }
})
