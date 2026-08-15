#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-markers ####
#
# analyzeData(parquetExprBase, scranMarkersParam) — streaming pairwise marker
# detection. Expression values in, one marker table per group out; nothing
# partial crosses a package boundary.
#
# Why this is cheap at all: scran's `findMarkers(test.type = "t")` never
# compares cells pairwise. It makes ONE pass over the matrix for
# per-(feature, group) n / mean / variance (`compute_blocked_stats_none` via
# beachmat) and does every pairwise comparison as arithmetic on a
# features x groups table — Welch error, Satterthwaite d.f., logFC, one-sided
# log p-values, then `combineMarkers`. The pass is a group-by aggregate, so it
# lowers into Acero via `.pe_accum_raw(by_cell = )`, and twenty clusters cost
# 380 comparisons that never touch the store.
#
# The pairwise algebra below is a transcription of scran's, deliberately
# duplicated rather than shared. It is a published statistic — a spec, not
# GiottoDisk machinery — and the alternative was a cross-package helper API
# that would have to be exported from Giotto and then kept stable, freezing
# the on-disk/in-mem boundary at whatever shape it had on the day it was
# written. Same call as the duplicated reader orchestration. So: do NOT
# refactor this into a shared utility; if it drifts from scran, the fix is to
# re-transcribe it, and the parity tests are what catch that.
#
# Correctness rests on absent entries meaning zero in the value space being
# aggregated, which the op registry maintains rather than this file checking:
# `multiply` has f(0) = 0, and `log` is log1p with `offset != 1` refused
# precisely to keep it. An `add` op would be the exception, which is why the
# roadmap pairs it with densification into triplet form — once the zero block
# is explicit, `sum` / `sumsq` over stored entries stay exact.
#
# No opinion is taken about which values these are. `expression_values` is
# resolved upstream by `getExpression()`, so this aggregates whatever chain the
# store carries. Unlike Pearson residual variance, a t-test is defined for any
# real-valued input, and a backend that imposed its own normalization would
# compute a different statistic from the in-memory one on the same input.


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
#' wired up. Materialize with `storeRead(x, output = "dgcmatrix")` for those.
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
                 "detection count that is not wired up yet. Workaround: ",
                 "materialize with storeRead(x, output = 'dgcmatrix') and ",
                 "use the in-memory path.", call. = FALSE)
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

# Per-(feature, group) accumulators, from the grouped feature-statistics verb.
# No aggregate of its own: `analyzeData(x, featStatsParam(), groups = )` is the
# statistic pass, and going through the verb rather than reaching past it means
# any backend that implements the verb supplies markers too.
#
# Asks for `sum` and `sumsq` only -- two of the four accumulators, so the
# `nnz` and `sum_det` work is never done. Those two plus the group's cell count
# are the complete minimal sufficient statistic for a Welch t-test.
#
# Returns the raw accumulators rather than mean/sd. Deriving the moments here
# is trivial, and keeping `sum`/`sumsq` lets `.pe_pool_moments()` combine
# groups by plain addition instead of reconstructing them through a square
# root.
#
# Absent entries contribute nothing to either accumulator, so they plus the
# group's FULL cell count give the dense mean and variance exactly. `n` comes
# from the group assignment, never the aggregate, and the grouped verb emits
# the complete cross product so a feature with no stored value in a group
# arrives as a zero row rather than a missing one.
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

# Every ordered pair of groups, Welch-tested on the moments, then combined into
# one table per group by `scran::combineMarkers()` — which is reused rather
# than transcribed, because it never touches expression values.
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
    n <- mom$n
    means <- mom$means
    vars <- mom$vars
    lvls <- colnames(means)
    feats <- rownames(means)
    direction <- param$direction %null% "any"
    lfc_thresh <- param$lfc %null% 0
    std_lfc <- isTRUE(param$std_lfc)

    # scran's pair order: host in level order, target in level order within
    # it. `.reorder_pairwise_output()` produces exactly this, and matching it
    # keeps the per-comparison column order comparable between backends.
    grid <- expand.grid(target = lvls, host = lvls, stringsAsFactors = FALSE)
    grid <- grid[grid$host != grid$target, , drop = FALSE]
    grid <- grid[
        order(match(grid$host, lvls), match(grid$target, lvls)), ,
        drop = FALSE
    ]

    stat_list <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
        host <- grid$host[i]
        target <- grid$target[i]

        tt <- .pe_welch(
            host_s2 = vars[, host], target_s2 = vars[, target],
            host_n = n[[host]], target_n = n[[target]]
        )
        cur_lfc <- means[, host] - means[, target]
        p_out <- .pe_run_t(cur_lfc, tt$err, tt$df, thresh_lfc = lfc_thresh)

        effect <- cur_lfc
        if (isTRUE(std_lfc)) {
            pooled_s2 <- ((n[[host]] - 1) * vars[, host] +
                (n[[target]] - 1) * vars[, target]) /
                (n[[host]] + n[[target]] - 2)
            is_zero <- effect == 0
            effect <- effect / sqrt(pooled_s2)
            effect[is_zero] <- 0
        }

        stat_list[[i]] <- .pe_full_stats(
            effect = effect,
            p = .pe_choose_lr(p_out$left, p_out$right, direction),
            feats = feats
        )
    }

    scran::combineMarkers(
        de.lists = stat_list,
        pairs = S4Vectors::DataFrame(
            first = grid$host, second = grid$target
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
