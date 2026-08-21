# Tests for the grouped form of analyzeData(parquetExprBase, featStatsParam):
# per-(feature, group) accumulators instead of one column of totals.
#
# This is the seam Giotto's gini markers ride on, so the statistics have to
# match an in-memory reference, and a grouping has to mean the same set of
# cells on both backends. The alignment cases are the point of the file: a
# per-cell vector is a payload, and adr/0003 says a payload keyed by view
# position reads the wrong entries once `[` narrows the store.

.grp_mat <- function(n_genes = 12L, n_cells = 30L, density = 0.5, seed = 3L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}

# Reference statistics, cells of each group selected by name.
.ref_grouped <- function(mat, groups) {
    lvls <- levels(droplevels(factor(groups)))
    pick <- function(k) colnames(mat) %in% names(groups)[!is.na(groups) &
        groups == k]
    list(
        lvls = lvls,
        mean = as.numeric(vapply(lvls, function(k) {
            Matrix::rowMeans(mat[, pick(k), drop = FALSE])
        }, numeric(nrow(mat)))),
        nnz = as.numeric(vapply(lvls, function(k) {
            Matrix::rowSums(mat[, pick(k), drop = FALSE] > 0)
        }, numeric(nrow(mat))))
    )
}

.fsg <- function(pe, groups, stats = c("sum", "nnz")) {
    GiottoClass::analyzeData(pe, Giotto::analyzeParam("feat_stats"),
        groups = groups, stats = stats)
}


test_that("grouped featStats matches an in-memory reference", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    grp <- stats::setNames(rep(c("a", "b", "c"), length.out = ncol(mat)),
        colnames(mat))
    ref <- .ref_grouped(mat, grp)
    st <- .fsg(pe, grp)

    # full cross product, groups slowest -- the order the callers reshape on
    expect_equal(nrow(st), nrow(mat) * length(ref$lvls))
    expect_equal(st$feats[seq_len(nrow(mat))], pe@feat_ids)
    expect_equal(unique(st$group[seq_len(nrow(mat))]), ref$lvls[[1L]])

    expect_equal(st$mean_expr, ref$mean)
    expect_equal(as.numeric(st$nr_cells), ref$nnz)
})


test_that("a named grouping is matched on cell ID, not on view position", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(seed = 5L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    grp <- stats::setNames(rep(c("a", "b"), length.out = ncol(mat)),
        colnames(mat))

    # Same assignment, shuffled: identity keying has to make this invariant,
    # where a positional read would repartition the cells entirely.
    set.seed(99L)
    shuffled <- grp[sample(length(grp))]
    expect_false(identical(names(shuffled), pe@cell_ids))

    expect_equal(.fsg(pe, shuffled)$mean_expr, .fsg(pe, grp)$mean_expr)
    expect_equal(.fsg(pe, shuffled)$mean_expr, .ref_grouped(mat, grp)$mean)
})


test_that("a named grouping survives a `[`-narrowed view", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(seed = 8L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    grp <- stats::setNames(rep(c("a", "b"), length.out = ncol(mat)),
        colnames(mat))

    # `[` rewrites the cell-id -> on-disk-key mapping; the grouping is built
    # against the whole object and must keep working unchanged against the
    # narrower view, covering cells the view no longer holds.
    keep <- seq(2L, ncol(mat), by = 3L)
    sub <- pe[, keep]
    expect_lt(length(sub@cell_ids), length(grp))

    got <- .fsg(sub, grp)
    ref <- .ref_grouped(mat[, sub@cell_ids, drop = FALSE],
        grp[sub@cell_ids])

    expect_equal(got$mean_expr, ref$mean)
    expect_equal(as.numeric(got$nr_cells), ref$nnz)

    # and the narrowed answer is genuinely not the full-view one
    expect_false(isTRUE(all.equal(got$mean_expr, .fsg(pe, grp)$mean_expr)))
})


test_that("cells the grouping omits drop out rather than erroring", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(seed = 13L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    grp <- stats::setNames(rep(c("a", "b", "c"), length.out = ncol(mat)),
        colnames(mat))

    # Masking the rest is how a caller narrows to a few groups, so dropping
    # every cell of one group removes that group and leaves the others alone.
    got <- .fsg(pe, grp[grp != "c"])
    expect_setequal(unique(got$group), c("a", "b"))

    full <- .fsg(pe, grp)
    expect_equal(got[got$group == "a"]$mean_expr, full[full$group == "a"]$mean_expr)

    # Names that match nothing is a mistake, not an empty selection.
    expect_error(
        .fsg(pe, stats::setNames(grp, paste0("nope_", names(grp)))),
        "none of its names"
    )
})


test_that("an unnamed grouping stays positional, with a warning", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(seed = 21L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    grp <- stats::setNames(rep(c("a", "b"), length.out = ncol(mat)),
        pe@cell_ids)

    expect_warning(got <- .fsg(pe, unname(grp)), "unnamed")
    expect_equal(got$mean_expr, .fsg(pe, grp)$mean_expr)

    expect_error(
        suppressWarnings(.fsg(pe, unname(grp)[-1L])), "one entry per cell"
    )
})


# ---- cell windowing -------------------------------------------------------
# The grouped pass windows its input because `by_cell` puts a join in front of
# the aggregate whose output is O(nonzeros) (see `.pe_accum_acero_windowed()`).
# Every fixture above fits one window, so these pin the window small enough to
# force a fold across several and assert it is exact, not merely close.

test_that("a grouped pass folds exactly across forced cell windows", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(n_genes = 20L, n_cells = 90L, density = 0.4, seed = 31L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    grp <- stats::setNames(rep(c("a", "b", "c"), length.out = ncol(mat)),
        colnames(mat))
    ref <- .ref_grouped(mat, grp)
    one <- .fsg(pe, grp, stats = c("sum", "sumsq", "nnz"))

    # 90 cells / 13 -> 7 windows. Sums, sums of squares and detection counts
    # are all additive over cells, so windowing must be bit-for-bit invariant.
    for (win in c(13L, 30L, 89L)) {
        withr::local_options(list(giottodisk.chunk_size = win))
        got <- .fsg(pe, grp, stats = c("sum", "sumsq", "nnz"))
        expect_equal(got, one)
        expect_equal(got$mean_expr, ref$mean)
        expect_equal(as.numeric(got$nr_cells), ref$nnz)
    }
})


test_that("windowing a grouped pass is invariant to `[` and to NA groups", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(n_genes = 16L, n_cells = 80L, density = 0.4, seed = 37L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    grp <- stats::setNames(rep(c("a", "b"), length.out = ncol(mat)),
        colnames(mat))
    # NA-group cells are excluded by the join, so they must not reappear when
    # the join is executed once per window instead of once overall.
    grp[seq(1L, length(grp), by = 4L)] <- NA

    sub <- pe[, seq(2L, ncol(mat), by = 2L)]
    one <- .fsg(sub, grp, stats = c("sum", "nnz"))
    withr::local_options(list(giottodisk.chunk_size = 9L))
    expect_equal(.fsg(sub, grp, stats = c("sum", "nnz")), one)
})


test_that("an ungrouped pass is not windowed", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(n_genes = 12L, n_cells = 40L, seed = 41L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # Acero bounds the joinless aggregate on its own; routing it through the
    # window machinery would cost setup per window for no bound.
    called <- 0L
    trace(".pe_accum_acero_windowed", where = asNamespace("GiottoDisk"),
        print = FALSE, tracer = function() called <<- called + 1L)
    on.exit(untrace(".pe_accum_acero_windowed",
        where = asNamespace("GiottoDisk")), add = TRUE)

    # No `stats =` here: the ungrouped verb has a fixed column contract and
    # rejects an accumulator selection.
    withr::local_options(list(giottodisk.chunk_size = 7L))
    st <- GiottoClass::analyzeData(pe, Giotto::analyzeParam("feat_stats"))
    expect_identical(called, 0L)
    expect_equal(st$total_expr, unname(Matrix::rowSums(mat)))
})


test_that("the data.table path folds exactly across cell windows", {
    skip_if_not_installed("Giotto")
    mat <- .grp_mat(n_genes = 18L, n_cells = 84L, density = 0.4, seed = 43L)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # A `log` on @post_ops forces the R-side executor: the chain cannot be
    # lowered, so `.pe_accum_chunked_dt()` runs instead of the Acero plan. That
    # path windows too, and folds its partials the same way -- before this it
    # was the only windowed path with no multi-window coverage.
    pe_post <- .pe_push_op(pe, list(type = "log", base = 2), phase = "post")
    expect_length(pe_post@post_ops, 1L)

    grp <- stats::setNames(rep(c("a", "b", "c"), length.out = ncol(mat)),
        colnames(mat))
    one <- .fsg(pe_post, grp, stats = c("sum", "sumsq", "nnz"))

    for (win in c(11L, 29L, 83L)) {
        withr::local_options(list(giottodisk.chunk_size = win))
        expect_equal(.fsg(pe_post, grp, stats = c("sum", "sumsq", "nnz")), one)
    }

    # Same statistic, computed in memory on the transformed values.
    ref <- .ref_grouped(log2(mat + 1), grp)
    expect_equal(one$mean_expr, ref$mean)
})
