#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-markers ####
#
# analyzeData(parquetExprBase, scranMarkersParam) — streaming pairwise marker
# detection: expression values in, one marker table per group out.
#
# Why it is cheap: scran's `findMarkers(test.type = "t")` never compares cells
# pairwise. It takes ONE pass for per-(feature, group) n / mean / variance, and
# every comparison after that is arithmetic on a features x groups table. That
# pass is a grouped aggregate, so it goes through the shared accumulator
# (design.Rmd, "Expression Statistics") and twenty clusters cost 380
# comparisons that never touch the store.
#
# The pairwise algebra below transcribes scran's, deliberately duplicated
# rather than shared: it is a published statistic — a spec, not GiottoDisk
# machinery — and the alternative was a cross-package helper API that would
# freeze the on-disk/in-mem boundary. Do NOT refactor it into a shared utility.
# If it drifts from scran, re-transcribe it; the parity tests are what catch
# that.
#
# Correctness rests on absent entries meaning zero, which the op registry
# maintains rather than this file checking — see the `add` note in
# R/utils-pestore-ops.R.
#
# Nothing here takes a view on WHICH values these are; `expression_values` is
# resolved upstream by `getExpression()`. Unlike Pearson residual variance, a
# t-test is defined for any real-valued input, so imposing a normalization
# would make this compute a different statistic from the in-memory backend.
#
# Memory: nothing on this path materializes a matrix, so `storeRead()`'s dgc
# slice cap is not a consideration — the `dgcmatrix` in the fallback advice
# below is about the in-memory path a user drops to. Frames are built and
# released per host, which bounds the live set but NOT peak RSS; do not "fix"
# that with a `gc()` in the loop, which halves peak at ~2.7x the runtime.
# Numbers in roadmap.Rmd. Cell count does not enter at all; high cluster counts
# are quadratic in time, where `comparison = "one_vs_rest"` is cheaper.


#' @name analyzeData-scranMarkersParam
#' @rdname analyzeData-scranMarkersParam
#' @title Streaming pairwise marker detection
#' @description
#' [GiottoClass::analyzeData()] method for a [Giotto::scranMarkersParam-class]
#' on a disk-backed expression store.
#'
#' The per-group statistic pass is a grouped aggregate pushed into Acero, so
#' the store is scanned once and the pairwise comparisons that follow never
#' touch it. Results are verified elementwise against
#' \code{\link[scran]{findMarkers}}.
#'
#' `test_type` must be `"t"`: it is the only test that reduces to per-group
#' moments. `"wilcox"` needs the per-cell values to rank, and `"binom"` is not
#' wired up.
#'
#' For those, fall back to the in-memory path. That means materializing the
#' whole store, **so it works only if the whole matrix fits in memory**:
#'
#' ```
#' m <- storeRead(x, output = "dgcmatrix", max_rows = Inf, max_cols = Inf)
#' analyzeData(m, markersParam(method = "scran", test_type = "wilcox"),
#'             groups = groups)
#' ```
#'
#' The cap arguments are required: `storeRead()` refuses a full-size
#' materialization unless one axis is small, which is what forces that fit to
#' be considered rather than discovered. Narrowing with `[` is not an
#' alternative here — a rank test needs every gene and every cell.
#' @param x a `parquetExprBase` store.
#' @param param a [Giotto::scranMarkersParam-class].
#' @param groups vector of cluster assignments, one per cell of the store's
#'   current view. `NA` excludes a cell from every group.
#' @param ... additional arguments (none used).
#' @returns A `SimpleList` of `DataFrame`s, one per group, as
#'   \code{\link[scran]{findMarkers}} returns.
#' @seealso [Giotto::markersParam()], [Giotto::findScranMarkers()]
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "scranMarkersParam"),
    function(x, param, ..., groups = NULL) {
        # Argument validation first: a caller who passed the wrong test_type
        # should hear about that, not be told to install a package they will
        # not end up needing.
        if (is.null(groups)) {
            stop("[analyzeData(parquetExprBase, scranMarkersParam)] `groups` is ",
                 "required: a cluster assignment per cell of the current ",
                 "view.", call. = FALSE)
        }
        test_type <- param$test_type %null% "t"
        if (!identical(test_type, "t")) {
            stop("[analyzeData(parquetExprBase, scranMarkersParam)] test_type = '",
                 test_type, "' is not available on the streaming backend. ",
                 "Only 't' reduces to per-group moments; 'wilcox' needs the ",
                 "per-cell values to rank, and 'binom' needs a per-group ",
                 "detection count that is not wired up yet. Fall back to the ",
                 "in-memory path, which materializes the whole store and so ",
                 "works only if the whole matrix fits in memory: ",
                 "storeRead(x, output = 'dgcmatrix', max_rows = Inf, ",
                 "max_cols = Inf). The cap arguments are required -- ",
                 "storeRead() refuses a full-size materialization unless one ",
                 "axis is small -- and `[`-subsetting is not an alternative, ",
                 "since a rank test needs every gene and every cell.",
                 call. = FALSE)
        }

        # Needed only by the pairwise tail, which reuses
        # `scran::combineMarkers()`. The statistic pass needs neither.
        package_check(pkg_name = "scran", repository = "Bioc")
        package_check(pkg_name = "S4Vectors", repository = "Bioc")

        rlang::inform(
            paste0(
                "[GiottoDisk] Streaming marker detection: pairwise ",
                "t-statistics are computed by a streaming equivalent of ",
                "scran's, verified against it; only combineMarkers() is ",
                "scran's own on this path."
            ),
            .frequency = "regularly",
            .frequency_id = "giottodisk_markers_streaming"
        )

        # ONE statistic pass, whichever comparison was asked for. One-vs-rest
        # then pools in memory -- the accumulators are additive, so the rest
        # group costs arithmetic rather than another scan.
        mom <- .pe_group_moments(x, groups)

        if (identical(param$comparison %null% "pairwise", "one_vs_rest")) {
            return(.pe_markers_one_vs_rest(mom, param))
        }
        .pe_markers_from_moments(.pe_moments_derive(mom), param = param)
    }
)


# ---- the statistic pass ----------------------------------------------------

# Per-(feature, group) accumulators, from the grouped feature-statistics verb
# rather than an aggregate of its own -- going through the verb means any
# backend implementing it supplies markers too.
#
# Asks for `sum` and `sumsq` only, so the `nnz` / `sum_det` work is never done;
# those two plus the group's cell count are the complete minimal sufficient
# statistic for a Welch t-test. Returned raw rather than as mean/sd so
# `.pe_pool_moments()` can combine groups by addition instead of reconstructing
# them through a square root.
#
# `n` comes from the group assignment, never the aggregate: a feature with no
# stored value in a group still has a mean over that group's full cell count.
.pe_group_moments <- function(pe, groups) {
    if (!inherits(pe, "parquetExprBase")) {
        stop("[.pe_group_moments] pe must be a parquetExprBase.",
             call. = FALSE)
    }
    # Length before level count: a mismatched `groups` is the more fundamental
    # error, and checking levels first would report a store of 200 cells as
    # having one group when really the caller passed 3 labels.
    n_cells <- as.integer(pe@n_cells)
    if (length(groups) != n_cells) {
        stop("[.pe_group_moments] `groups` must have one entry per cell of ",
             "the current view (", n_cells, "), got ", length(groups), ".",
             call. = FALSE)
    }
    g <- droplevels(if (is.factor(groups)) groups else factor(groups))
    if (nlevels(g) < 2L) {
        stop("[.pe_group_moments] need at least 2 non-empty groups, got ",
             nlevels(g), ".", call. = FALSE)
    }

    # `featStatsParam` is a class, not a constructor -- `analyzeParam()` is the
    # exported factory for the whole analyze-param family.
    st <- GiottoClass::analyzeData(
        pe, Giotto::analyzeParam("feat_stats"),
        groups = g, stats = c("sum", "sumsq")
    )

    lvls <- levels(g)
    feats <- pe@feat_ids
    shape <- function(v) {
        matrix(v, nrow = length(feats), ncol = length(lvls),
            dimnames = list(feats, lvls))
    }
    # The grouped verb emits groups slowest with features cycling within --
    # column-major order -- so the columns land without a reshape.
    list(
        n = stats::setNames(st$n_cells[match(lvls, st$group)], lvls),
        sums = shape(st$total_expr),
        sumsq = shape(st$sumsq)
    )
}


# Combine groups into larger ones. The accumulators are additive, so this is a
# row-sum over the chosen columns -- no reconstruction, no square root, exact.
# That is the whole reason `.pe_group_moments()` keeps `sums`/`sumsq` instead
# of means and standard deviations.
#
# `sets` is a named list of character vectors naming the groups to merge.
.pe_pool_moments <- function(mom, sets) {
    lvls <- colnames(mom$sums)
    miss <- setdiff(unlist(sets), lvls)
    if (length(miss) > 0L) {
        stop("[.pe_pool_moments] groups not present: ",
             paste(unique(miss), collapse = ", "), call. = FALSE)
    }
    pick <- function(m, gs) {
        if (length(gs) == 1L) m[, gs] else rowSums(m[, gs, drop = FALSE])
    }
    list(
        n = stats::setNames(
            vapply(sets, function(gs) sum(mom$n[gs]), numeric(1L)), names(sets)
        ),
        sums = vapply(sets, function(gs) pick(mom$sums, gs),
            numeric(nrow(mom$sums))),
        sumsq = vapply(sets, function(gs) pick(mom$sumsq, gs),
            numeric(nrow(mom$sumsq)))
    )
}


# Accumulators -> the moments the pairwise tail consumes. Kept separate from
# the pass so pooling can happen first, on the additive form.
.pe_moments_derive <- function(mom) {
    n <- mom$n
    nn <- rep(n, each = nrow(mom$sums))
    means <- mom$sums / rep(n, each = nrow(mom$sums))
    # Clamped for the same reason as the grouped verb: cancellation can push
    # the identity slightly negative and produce an NaN t-statistic.
    vars <- ifelse(nn > 1,
        pmax((mom$sumsq - mom$sums * mom$sums / nn) / (nn - 1), 0), 0)
    dim(vars) <- dim(mom$sumsq)
    dimnames(means) <- dimnames(vars) <- dimnames(mom$sums)
    list(n = n, means = means, vars = vars)
}


# ---- the pairwise tail (transcribed from scran) -----------------------------

# One table per group, each testing that group against the pooled remainder.
#
# Level naming matches what `findScranMarkers(group_1 =, group_2 =)` produces
# in memory, because `findScranMarkers_one_vs_all()`'s post-processing selects
# by it: group 1 is the cluster's own name, group 2 is the other clusters
# pasted with "_".
.pe_markers_one_vs_rest <- function(mom, param) {
    lvls <- colnames(mom$sums)
    out <- lapply(lvls, function(k) {
        rest <- setdiff(lvls, k)
        sets <- stats::setNames(
            list(k, rest), c(k, paste0(rest, collapse = "_"))
        )
        res <- .pe_markers_from_moments(
            .pe_moments_derive(.pe_pool_moments(mom, sets)), param = param
        )
        res[[k]]
    })
    names(out) <- lvls
    S4Vectors::SimpleList(out)
}


.pe_markers_from_moments <- function(mom, param) {
    lvls <- colnames(mom$means)

    # Per host, not all at once. `combineMarkers()` takes the whole
    # `G(G-1)`-long list but immediately splits it on `pairs[, 1]`, and each
    # host's combination reads only its own frames. So calling it once per host
    # is the same computation with `G-1` frames live instead of `G(G-1)`.
    #
    # This is the split scran already performs internally, not a
    # reimplementation of the combination.
    out <- lapply(lvls, function(host) {
        # `setdiff` keeps level order, which is scran's target order within a
        # host (`.reorder_pairwise_output()` sorts by host then target). Column
        # order of the per-comparison stats depends on it.
        targets <- setdiff(lvls, host)
        frames <- lapply(targets, function(target) {
            .pe_pair_stats(mom, host, target, param)
        })
        # `combineMarkers()` checks that rownames agree across the frames it is
        # given and errors otherwise, so the assertion is already made for us.
        res <- .pe_combine(frames, host, targets, param)
        res[[host]]
    })
    names(out) <- lvls
    S4Vectors::SimpleList(out)
}


# One ordered comparison, host vs target, as the statistic frame
# `combineMarkers()` consumes. Three dense numeric vectors of length n_genes --
# no expression values reach this point.
.pe_pair_stats <- function(mom, host, target, param) {
    n <- mom$n
    means <- mom$means
    vars <- mom$vars

    tt <- .pe_welch(
        host_s2 = vars[, host], target_s2 = vars[, target],
        host_n = n[[host]], target_n = n[[target]]
    )
    cur_lfc <- means[, host] - means[, target]
    p_out <- .pe_run_t(cur_lfc, tt$err, tt$df,
        thresh_lfc = param$lfc %null% 0)

    effect <- cur_lfc
    if (isTRUE(param$std_lfc)) {
        pooled_s2 <- ((n[[host]] - 1) * vars[, host] +
            (n[[target]] - 1) * vars[, target]) /
            (n[[host]] + n[[target]] - 2)
        is_zero <- effect == 0
        effect <- effect / sqrt(pooled_s2)
        effect[is_zero] <- 0
    }

    .pe_full_stats(
        effect = effect,
        p = .pe_choose_lr(p_out$left, p_out$right,
            param$direction %null% "any"),
        feats = rownames(means)
    )
}


# One host's comparisons -> its marker table. Returns the one-element
# `SimpleList` `combineMarkers()` gives back for a single `first` level.
.pe_combine <- function(frames, host, targets, param) {
    scran::combineMarkers(
        de.lists = frames,
        pairs = S4Vectors::DataFrame(
            first = rep(host, length(targets)), second = targets
        ),
        pval.field = "log.p.value",
        effect.field = "logFC",
        pval.type = param$pval_type %null% "any",
        min.prop = param$min_prop,
        log.p.in = TRUE,
        log.p.out = isTRUE(param$log_p),
        full.stats = isTRUE(param$full_stats),
        sorted = param$sorted %null% TRUE
    )
}

# Welch standard error and Satterthwaite d.f. Transcribes scran's
# `.get_t_test_stats()`, variance floor included — the floor is why a
# scale-changing transform of the values is not a no-op for near-constant
# features even where the t-statistic is otherwise scale-invariant.
.pe_welch <- function(host_s2, target_s2, host_n, target_n) {
    host_df <- max(0L, host_n - 1L)
    target_df <- max(0L, target_n - 1L)
    host_s2 <- pmax(host_s2, 1e-8)
    target_s2 <- pmax(target_s2, 1e-8)
    if (host_df > 0L && target_df > 0L) {
        host_err <- host_s2 / host_n
        target_err <- target_s2 / target_n
        cur_err <- host_err + target_err
        cur_df <- cur_err^2 /
            (host_err^2 / host_df + target_err^2 / target_df)
    } else {
        cur_err <- cur_df <- NA_real_
    }
    list(err = cur_err, df = cur_df)
}

# Two one-sided log p-values. Transcribes scran's `.run_t_test()`.
.pe_run_t <- function(cur_lfc, cur_err, cur_df, thresh_lfc = 0) {
    thresh_lfc <- abs(thresh_lfc)
    if (thresh_lfc == 0) {
        cur_t <- cur_lfc / sqrt(cur_err)
        left <- stats::pt(cur_t, df = cur_df, lower.tail = TRUE, log.p = TRUE)
        right <- stats::pt(cur_t, df = cur_df, lower.tail = FALSE,
            log.p = TRUE)
    } else {
        lower_t <- (cur_lfc + thresh_lfc) / sqrt(cur_err)
        left <- stats::pt(lower_t, df = cur_df, lower.tail = TRUE,
            log.p = TRUE)
        upper_t <- (cur_lfc - thresh_lfc) / sqrt(cur_err)
        right <- stats::pt(upper_t, df = cur_df, lower.tail = FALSE,
            log.p = TRUE)
    }
    list(left = left, right = right)
}

# Transcribes scran's `.choose_leftright_pvalues()`.
.pe_choose_lr <- function(left, right, direction = "any") {
    switch(direction,
        "up" = right,
        "down" = left,
        pmin(0, pmin(left, right) + log(2))
    )
}

# Benjamini-Hochberg in log space. Transcribes scran's `.logBH()`.
.pe_logBH <- function(log_p) {
    o <- order(log_p)
    repval <- log_p[o] + log(length(o) / seq_along(o))
    repval <- rev(cummin(rev(repval)))
    repval[o] <- repval
    repval
}

# Per-comparison stats frame. Transcribes scran's `.create_full_stats(log.p =
# TRUE)`. `combineMarkers()` reads `logFC` and `log.p.value` from it and
# surfaces `log.FDR` only when `full_stats = TRUE`.
.pe_full_stats <- function(effect, p, feats) {
    p <- as.vector(p)
    S4Vectors::DataFrame(
        logFC = as.vector(effect),
        log.p.value = p,
        log.FDR = .pe_logBH(p),
        check.names = FALSE,
        row.names = feats
    )
}
