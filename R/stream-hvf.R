#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-hvf ####
# Streaming HVF-relevant stats for parquetExprStore-backed expression.
# Dispatches via Giotto's analyzeData(x, param) generic. Methods return
# the per-feature stats data.table without performing selection;
# downstream thresholding / selection is a separate step.
#
#   analyzeData(parquetExprStore, covLoessParam)
#       -> per-feature stats including cov_diff (residual COV above
#          a LOESS fit of cov ~ log(mean_expr))
#   analyzeData(parquetExprStore, covGroupsParam)
#       -> per-feature stats including cov_group_zscore (within-bin
#          COV z-score)
#   analyzeData(parquetExprStore, varParam)
#       -> per-feature variance of analytic Poisson Pearson residuals,
#          computed from raw counts (see the varParam section below)
#
# The COV methods share one stats pass, `.stream_gene_stats()`, which
# reduces the triplet stream to n_genes-sized sum / sumsq / nnz vectors;
# LOESS and bin-zscore then run on those. That pass has two shapes: an
# arrow-native aggregate when the recipe is a single norm_libsize_log
# post-op, and a generic per-substore R-side pass otherwise.

# ---- covLoessParam: streaming ---------------------------------------------

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "covLoessParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        stats <- .stream_gene_stats(x, expression_threshold = thr)

        # Match Giotto: drop zero-detection features before fitting
        nr_cells <- cov <- pred_cov <- cov_diff <- mean_expr <- NULL
        stats <- stats[nr_cells > 0]

        loess_fit <- stats::loess(cov ~ log(mean_expr), data = stats)
        stats[, pred_cov := stats::predict(loess_fit, newdata = stats)]
        stats[, cov_diff := cov - pred_cov]
        stats[, pred_cov := NULL]
        data.table::setorder(stats, -cov_diff)
        stats
    }
)


# ---- covGroupsParam: streaming --------------------------------------------

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "covGroupsParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        stats <- .stream_gene_stats(x, expression_threshold = thr)

        # NSE bindings
        nr_cells <- cov <- expr_groups <- cov_group_zscore <- NULL

        # Match Giotto: drop zero-detection features before binning
        stats <- stats[nr_cells > 0]

        # Quantile-bin by mean expression. If too many tied breaks (lots of
        # zero-mean genes), recompute on the strictly positive subset and
        # set the leading break to 0 so all-zero genes still bin into
        # group_1 -- matches Giotto's in-memory .calc_cov_group_hvf.
        n_groups <- as.integer(param$nr_expression_groups)
        prob_seq <- seq(0, 1, by = 1 / n_groups)
        prob_seq[length(prob_seq)] <- 1
        expr_group_breaks <- stats::quantile(stats$mean_expr, probs = prob_seq)
        if (any(duplicated(expr_group_breaks))) {
            m <- stats$mean_expr
            expr_group_breaks <- stats::quantile(m[m > 0], probs = prob_seq)
            expr_group_breaks[[1L]] <- 0
        }

        expr_groups_lbl <- cut(
            stats$mean_expr,
            breaks         = expr_group_breaks,
            labels         = paste0("group_", seq_len(n_groups)),
            include.lowest = TRUE
        )
        stats[, expr_groups := expr_groups_lbl]
        stats[, cov_group_zscore := scale(cov), by = expr_groups]
        stats[, expr_groups := NULL]
        stats
    }
)


# ---- varParam on parquet: analytic Pearson residual variance --------------
#
# `calculateHVF(method = "var_p_resid")` ranks features by the `var` column
# this returns.  In Giotto's in-memory path that is `rowVars()` of whatever
# matrix it was handed, which only equals Pearson residual variance if the
# caller arranged for the expression slot to hold residuals -- with the
# default `expression_values = "normalized"` it is the variance of
# log-normalized values instead.
#
# The streaming backend computes the residual variance directly from RAW
# counts, which is both what the method name means and the only form that
# works here: Pearson residuals are dense (every zero maps to a nonzero
# residual), so a residual matrix would be n_genes x n_cells -- 3.1e9 rows
# for Atera-scale data.  See `.stream_pearson_resid_var()` for how the zero
# block is folded in analytically instead of being materialized.
#
# Consequences worth knowing:
#   * `expression_values` defaults to "raw" here, not "normalized".  Pearson
#     residuals ARE the normalization -- the depth term lives inside
#     `mu = g_i * c_j / T` -- so normalizing first double-corrects.  Measured
#     on synthetic data with 20 injected overdispersed genes: raw counts
#     recover 20/20 with non-variable genes calibrated at 1.000, whereas
#     running the same formula on log2(1+libnorm) values recovers 0/20.
#   * the scale is absolute: 1.0 means "no more variable than Poisson
#     sampling noise", which is what makes `calculateHVF`'s default
#     `var_threshold = 1.5` meaningful rather than arbitrary.

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "varParam"),
    function(x, param, ...) {
        ev <- param$expression_values %null% "raw"
        if (!identical(ev, "raw")) {
            warning("[analyzeData(parquetExprBase, varParam)] ",
                    "expression_values = '", ev, "' is ignored: Pearson ",
                    "residual variance is defined on raw counts, so the ",
                    "streaming backend reads raw values and does not apply ",
                    "the normalization recipe.", call. = FALSE)
        }
        if (!inherits(x, "parquetExprStore")) {
            stop("[analyzeData(parquetExprBase, varParam)] union stores are ",
                 "not supported yet: per-cell totals key on ",
                 "(source_id, row_id) across substores and the gene-axis ",
                 "remap resolves against a single store's @gene_idx.",
                 call. = FALSE)
        }
        .stream_pearson_resid_var(x, size_factors = param$size_factors)
    }
)


# Per-gene variance of analytic (Poisson) Pearson residuals, streamed.
#
#   mu_ij = g_i * c_j / T        z_ij = (x_ij - mu_ij) / sqrt(mu_ij)
#
# with g_i the gene total, c_j the cell total (or a supplied size factor),
# and T the grand total.  Residuals are dense, but the all-zero block folds
# into closed form: treat every cell as zero, then correct only the stored
# nonzeros.
#
#   zeros:  sum_j z_ij^2 = sum_j mu_ij = g_i
#           sum_j z_ij   = -sqrt(g_i / T) * S,      S = sum_j sqrt(c_j)
#   stored nonzero (x > 0):
#           dz  = x / sqrt(mu_ij)
#           dz2 = x^2 / mu_ij - 2x
#
# So the whole statistic needs the gene totals, the cell totals, two scalars,
# and ONE joined pass over the existing triplet stream -- all pushed into
# Acero.  Validated against a dense reference to 2.7e-14; the same algebra
# runs on a dgCMatrix via @i/@p/@x if an in-memory version is ever wanted.
#
# @post_ops is deliberately NOT applied: the values must be counts.

#' @keywords internal
#' @noRd
.stream_pearson_resid_var <- function(pe, size_factors = NULL) {
    # NSE bindings
    row_id <- col_id <- value <- g <- cc <- mu <- dz <- dz2 <- NULL
    sum_dz <- sum_dz2 <- sum_z <- sum_z2 <- var <- feats <- NULL

    n_genes <- as.integer(pe@n_genes)

    # Pass 1: per-gene and per-cell totals (both pushed-down aggregates).
    gt <- storeRead(pe, output = "query") |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(g = sum(value, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()
    if (nrow(gt) == 0L) {
        return(data.table::data.table(feats = pe@feat_ids,
                                      var = numeric(n_genes)))
    }
    ct <- storeRead(pe, output = "query") |>
        dplyr::group_by(row_id) |>
        dplyr::summarise(cc = sum(value, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()

    # Counts sanity: the model assumes integer counts, and nothing about the
    # store can guarantee that (a normalized matrix may have been ingested).
    if (any(abs(gt$g - round(gt$g)) > 1e-8)) {
        warning("[analyzeData(varParam)] expression values do not look like ",
                "integer counts; Pearson residual variance assumes raw ",
                "counts and will be misleading on transformed values.",
                call. = FALSE)
    }

    # Cells with no counts have mu = 0 for every gene, leaving the residual
    # undefined; drop them and shrink n accordingly.
    if (!is.null(size_factors)) {
        sf <- as.numeric(size_factors)
        if (length(sf) != as.integer(pe@n_cells)) {
            stop("[analyzeData(varParam)] size_factors must have one entry ",
                 "per cell (", pe@n_cells, "), got ", length(sf), ".",
                 call. = FALSE)
        }
        cvals <- sf[sf > 0]
    } else {
        cvals <- ct$cc[ct$cc > 0]
    }
    n_eff <- length(cvals)
    if (n_eff < 2L) {
        return(data.table::data.table(feats = pe@feat_ids,
                                      var = numeric(n_genes)))
    }
    Tt <- sum(gt$g)
    S  <- sum(sqrt(cvals))

    # Pass 2: per-nonzero corrections, joined and aggregated in Acero.
    gt_a <- arrow::as_arrow_table(data.frame(
        col_id = as.integer(gt$col_id), g = as.numeric(gt$g)))
    ct_a <- arrow::as_arrow_table(data.frame(
        row_id = as.integer(ct$row_id), cc = as.numeric(ct$cc)))
    corr <- storeRead(pe, output = "query") |>
        dplyr::left_join(gt_a, by = "col_id") |>
        dplyr::left_join(ct_a, by = "row_id") |>
        dplyr::mutate(mu = g * cc / !!Tt) |>
        dplyr::mutate(dz  = value / sqrt(mu),
                      dz2 = value * value / mu - 2 * value) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(sum_dz  = sum(dz,  na.rm = TRUE),
                         sum_dz2 = sum(dz2, na.rm = TRUE)) |>
        dplyr::collect() |> data.table::as.data.table()

    res <- merge(gt, corr, by = "col_id", all.x = TRUE)
    res[is.na(sum_dz),  sum_dz  := 0]
    res[is.na(sum_dz2), sum_dz2 := 0]
    res[, sum_z  := -sqrt(g / Tt) * S + sum_dz]
    res[, sum_z2 := g + sum_dz2]
    res[, var := (sum_z2 - sum_z * sum_z / n_eff) / (n_eff - 1L)]
    # Genes with no counts have undefined residuals; report 0 variance, which
    # is what rowVars() gives for an all-zero row in the in-memory path.
    res[g <= 0, var := 0]

    out <- numeric(n_genes)
    idx <- .pe_remap_col(res$col_id, pe)
    keep <- !is.na(idx)
    out[idx[keep]] <- as.numeric(res$var[keep])

    dt <- data.table::data.table(feats = pe@feat_ids, var = out)
    data.table::setorder(dt, -var)
    dt
}


# Shared tail for both stats paths: derive mean / sd / cov from the summed
# accumulators and assemble the result. Variance via the SS identity,
# σ² = (Σx² − n·μ²)/(n−1).
#
# Also the one place that can cheaply tell whether the values were normalized.
# There is no requirement that they be -- a store may hold values normalized on
# disk, with an empty op chain -- but per-gene totals that are all integers
# almost certainly mean raw counts reached a statistic that is normally taken
# on normalized values. `gene_sum` is already computed, so the check is free.
.gene_stats_dt <- function(pe, gene_sum, gene_sumsq, gene_nnz,
                            n_cells, n_genes, .warn_raw = TRUE) {
    if (isTRUE(.warn_raw)) {
        nzs <- gene_sum[gene_sum > 0]
        if (length(nzs) > 0L && all(abs(nzs - round(nzs)) < 1e-8)) {
            warning("[analyzeData] expression values look like raw integer ",
                    "counts. COV-based feature statistics are normally taken ",
                    "on normalized values; run normalizeGiotto() first, or ",
                    "ignore this if the store already holds normalized data.",
                    call. = FALSE)
        }
    }

    gene_mean <- gene_sum / n_cells
    gene_var  <- if (n_cells > 1L) {
        pmax((gene_sumsq - gene_sum * gene_sum / n_cells) / (n_cells - 1L), 0)
    } else {
        numeric(n_genes)
    }
    gene_sd  <- sqrt(gene_var)
    gene_cov <- ifelse(gene_mean > 0, gene_sd / gene_mean, NaN)

    data.table::data.table(
        feats      = pe@feat_ids,
        nr_cells   = gene_nnz,
        total_expr = gene_sum,
        mean_expr  = gene_mean,
        sd         = gene_sd,
        cov        = gene_cov
    )
}


# ---- Internal: streaming per-gene stats, JIT-normalizing if a recipe -------
#
# No normalization requirement. Any @post_ops present are applied on the way
# through; an empty chain means the on-disk values are used as they are, which
# is correct for a store holding values normalized at write time. Only the
# arrow fast path below needs a specific op, and it tests for that itself.

.stream_gene_stats <- function(pe, expression_threshold = 0) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[.stream_gene_stats] pe must be a parquetExprBase.")

    thr     <- as.numeric(expression_threshold)
    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)

    # Arrow-native eligibility.  The whole aggregate runs in Acero and only
    # the per-gene result (~18k rows) crosses into R, which on Atera (170k
    # cells x 18k genes) is 1.94 s serial against 11.16 s for the best R-side
    # chunking shape we measured.  Requires the norm as arrow compute, so it
    # is limited to a single `norm_libsize_log` post-op; anything else falls
    # through to the generic R-side pass below, which handles any op kind.
    #
    # Union stores fall through too: the (source_id, orig_row_id) join would
    # span substores fine, but `.pe_remap_col` resolves against a single
    # store's @gene_idx.
    can_arrow <- inherits(pe, "parquetExprStore") &&
        (length(pe@post_ops) == 0L ||
         (length(pe@post_ops) == 1L &&
          identical(pe@post_ops[[1L]]$type, "norm_libsize_log")))

    if (can_arrow) {
        return(.stream_gene_stats_arrow(pe, thr,
            n_cells = n_cells, n_genes = n_genes))
    }

    # NSE bindings
    col_id <- value <- s <- s2 <- nz <- NULL

    gene_sum   <- numeric(n_genes)
    gene_sumsq <- numeric(n_genes)
    gene_nnz   <- integer(n_genes)

    # Iterate substores, applying the parent's @post_ops to each collected
    # chunk. Payloads carry source_id / orig_row_id, so per-cell state such
    # as norm scale factors resolves against the right substore.
    #
    # Values are normalized in place with no raw snapshot kept: `total_expr`
    # follows Giotto and reports `rowSums()` of the values it was handed, so
    # it is the NORMALIZED sum (identical to `mean_expr * n_cells`). A
    # `raw_value` copy would add a numeric column over the whole triplet
    # stream (~2.5 GB at 307M rows) for a quantity nothing reads.
    post_ops <- pe@post_ops
    subs <- .exprbase_substores(pe)
    for (sub_entry in subs) {
        sub <- sub_entry$store
        df <- storeRead(sub, output = "query") |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) next
        df <- .pe_apply_post_ops_df(df, post_ops)

        # Per-gene aggregation R-side (data.table by-group)
        agg <- df[, .(
            s  = sum(value, na.rm = TRUE),
            s2 = sum(value * value, na.rm = TRUE),
            nz = sum(value > thr, na.rm = TRUE)
        ), by = col_id]

        if (nrow(agg) > 0L) {
            # `.pe_remap_col` on the substore maps on-disk col_ids to
            # local feat positions; feat_ids align across substores
            # (union invariant), so the same local position indexes the
            # union's @feat_ids axis.
            idx <- .pe_remap_col(agg$col_id, sub)
            keep <- !is.na(idx)
            gene_sum[idx[keep]]       <- gene_sum[idx[keep]] +
                as.numeric(agg$s[keep])
            gene_sumsq[idx[keep]]     <- gene_sumsq[idx[keep]] +
                as.numeric(agg$s2[keep])
            gene_nnz[idx[keep]]       <- gene_nnz[idx[keep]] +
                as.integer(agg$nz[keep])
        }
    }

    .gene_stats_dt(pe, gene_sum, gene_sumsq, gene_nnz, n_cells, n_genes)
}


# Arrow-native per-gene stats.  Single Acero pass: the composed lazy query
# (subset filters + @ops via storeRead) is left-joined to the norm op's
# per-cell scalef table, the normalized value is materialized as a compute
# expression, and the per-gene aggregate is pushed down.  Only the ~n_genes
# result rows cross the arrow -> R boundary.
#
# The norm math here MIRRORS `.pe_apply_post_op_norm_libsize_log_df`; it is
# duplicated because arrow cannot index an R vector positionally inside a
# query, so the per-cell scalef has to arrive as a joinable table.  Keep the
# two in sync -- the sole reason this function is restricted to the
# `norm_libsize_log` op kind.
#
# `total_expr` is the normalized sum `s`, matching the generic pass above.

#' @keywords internal
#' @noRd
.stream_gene_stats_arrow <- function(pe, thr, n_cells, n_genes) {
    # NSE bindings
    row_id <- col_id <- value <- nv <- scalef <- source_id <- NULL
    s <- s2 <- nz <- NULL

    # Empty op chain: the on-disk values ARE the values, so the join and the
    # norm expression drop out and the aggregate runs directly. This is the
    # baked / normalized-at-write case, which would otherwise fall to the
    # generic R-side pass for no reason.
    q <- storeRead(pe, output = "query")
    if (length(pe@post_ops) == 0L) {
        q <- dplyr::mutate(q, nv = value)
    } else {
        op       <- pe@post_ops[[1L]]
        log_flag <- isTRUE(op$log)
        log_base <- op$base %null% 2

        sf <- data.table::as.data.table(op$scalef)
        scalef_tbl <- arrow::as_arrow_table(data.frame(
            source_id   = as.character(sf$source_id),
            orig_row_id = as.integer(sf$orig_row_id),
            scalef      = as.numeric(sf$scalef),
            stringsAsFactors = FALSE
        ))

        q <- q |>
            dplyr::left_join(scalef_tbl,
                by = c("source_id" = "source_id", "row_id" = "orig_row_id")) |>
            dplyr::mutate(nv = value * scalef)
        if (log_flag) {
            q <- dplyr::mutate(q, nv = log1p(nv) / log(!!log_base))
        }
    }

    agg <- q |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(
            s  = sum(nv, na.rm = TRUE),
            s2 = sum(nv * nv, na.rm = TRUE),
            nz = sum(as.integer(nv > !!thr), na.rm = TRUE)
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    gene_sum   <- numeric(n_genes)
    gene_sumsq <- numeric(n_genes)
    gene_nnz   <- integer(n_genes)

    if (nrow(agg) > 0L) {
        idx  <- .pe_remap_col(agg$col_id, pe)
        keep <- !is.na(idx)
        gene_sum[idx[keep]]   <- as.numeric(agg$s[keep])
        gene_sumsq[idx[keep]] <- as.numeric(agg$s2[keep])
        gene_nnz[idx[keep]]   <- as.integer(agg$nz[keep])
    }

    # With a norm recipe the values are normalized by construction and the
    # raw-counts check would be pure noise. With an empty op chain they are
    # whatever is on disk, so the check still applies.
    .gene_stats_dt(pe, gene_sum, gene_sumsq, gene_nnz, n_cells, n_genes,
                   .warn_raw = length(pe@post_ops) == 0L)
}
