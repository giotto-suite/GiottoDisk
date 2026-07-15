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
#
# Implementation: storeRead(pe) returns the arrow query with @ops already
# composed in (norm_libsize_log adds `v_norm`), then a per-gene
# group_by + summarise runs at arrow layer — only the small (n_genes-sized)
# aggregate transfers to R. LOESS / bin-zscore run on those small vectors.
#
# varParam still errors clearly because per-gene variance on a scaled
# (z-scored) matrix requires materialising the dense matrix.

# ---- covLoessParam: streaming ---------------------------------------------

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprStore", param = "covLoessParam"),
    function(x, param, ...) {
        if (!.pe_has_norm_op(x@ops)) {
            stop("[analyzeData(parquetExprStore, covLoessParam)] ",
                 "expression backend has no normalization recipe. Run ",
                 "normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE) ",
                 "first to populate scale factors on the store.",
                 call. = FALSE)
        }

        thr <- param$detection_threshold %null% 0
        stats <- .stream_norm_gene_stats(x, expression_threshold = thr)

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
    signature(x = "parquetExprStore", param = "covGroupsParam"),
    function(x, param, ...) {
        if (!.pe_has_norm_op(x@ops)) {
            stop("[analyzeData(parquetExprStore, covGroupsParam)] ",
                 "expression backend has no normalization recipe. Run ",
                 "normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE) ",
                 "first to populate scale factors on the store.",
                 call. = FALSE)
        }

        thr <- param$detection_threshold %null% 0
        stats <- .stream_norm_gene_stats(x, expression_threshold = thr)

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


# ---- varParam on parquet: clear error -------------------------------------

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprStore", param = "varParam"),
    function(x, param, ...) {
        stop("[analyzeData(parquetExprStore, varParam)] per-feature ",
             "variance on a scaled (z-scored) matrix requires ",
             "materialising the dense matrix and is not supported for ",
             "streaming backends. Use covLoessParam or covGroupsParam.",
             call. = FALSE)
    }
)


# ---- Internal: streaming per-gene stats with JIT normalization ------------

.stream_norm_gene_stats <- function(pe, expression_threshold = 0) {
    if (!inherits(pe, "parquetExprStore"))
        stop("[.stream_norm_gene_stats] pe must be a parquetExprStore.")

    if (!.pe_has_norm_op(pe@ops))
        stop("[.stream_norm_gene_stats] pe has no norm op on @ops.")
    thr <- as.numeric(expression_threshold)

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)

    # NSE bindings
    col_id <- value <- v_norm <- s <- s2 <- nz <- raw_total <- NULL

    # storeRead returns the arrow query with @ops composed in — v_norm is
    # already projected by the norm_libsize_log op. Per-gene aggregation
    # runs at arrow layer; only n_genes-sized result transfers to R.
    agg <- storeRead(pe, output = "query") |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(
            s         = sum(v_norm, na.rm = TRUE),
            s2        = sum(v_norm * v_norm, na.rm = TRUE),
            nz        = sum(v_norm > !!thr, na.rm = TRUE),
            raw_total = sum(value, na.rm = TRUE)
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    gene_sum   <- numeric(n_genes)
    gene_sumsq <- numeric(n_genes)
    gene_nnz   <- integer(n_genes)
    gene_total_raw <- numeric(n_genes)

    if (nrow(agg) > 0L) {
        idx <- .pe_remap_col(agg$col_id, pe)
        keep <- !is.na(idx)
        gene_sum[idx[keep]]       <- as.numeric(agg$s[keep])
        gene_sumsq[idx[keep]]     <- as.numeric(agg$s2[keep])
        gene_nnz[idx[keep]]       <- as.integer(agg$nz[keep])
        gene_total_raw[idx[keep]] <- as.numeric(agg$raw_total[keep])
    }

    gene_mean <- gene_sum / n_cells
    gene_var  <- if (n_cells > 1L) {
        pmax((gene_sumsq - gene_sum * gene_sum / n_cells) / (n_cells - 1L), 0)
    } else {
        numeric(n_genes)
    }
    gene_sd   <- sqrt(gene_var)
    gene_cov  <- ifelse(gene_mean > 0, gene_sd / gene_mean, NaN)

    data.table::data.table(
        feats      = pe@feat_ids,
        nr_cells   = gene_nnz,
        total_expr = gene_total_raw,
        mean_expr  = gene_mean,
        sd         = gene_sd,
        cov        = gene_cov
    )
}
