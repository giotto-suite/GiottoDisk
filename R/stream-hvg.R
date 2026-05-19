#' @include class-parquetExprStore.R
NULL

# stream-hvg ####
# Streaming HVG-relevant stats for parquetExprStore-backed expression.
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
# Implementation: one streaming Arrow pass that applies the JIT
# normalize recipe stored on pe@params$norm and accumulates per-gene
# sum + sum-of-squares, then computes mean / sd / cv on the small
# (n_genes-sized) result vectors and runs LOESS / bin-zscore in memory.
#
# varParam still errors clearly because per-gene variance on a scaled
# (z-scored) matrix requires materialising the dense matrix.

# ---- covLoessParam: streaming ---------------------------------------------

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprStore", param = "covLoessParam"),
    function(x, param, ...) {
        if (is.null(x@params$norm) ||
            is.null(x@params$norm$scale_factors)) {
            stop("[analyzeData(parquetExprStore, covLoessParam)] ",
                 "expression backend has no normalization recipe. Run ",
                 "normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE) ",
                 "first to populate scale factors on the store.",
                 call. = FALSE)
        }

        thr <- param$detection_threshold %null% 0
        stats <- .stream_norm_gene_stats(x, expression_threshold = thr)

        # Reuse Giotto's in-memory LOESS step on the (small) per-gene table
        cov <- pred_cov_feats <- cov_diff <- mean_expr <- NULL
        loess_model <- stats::loess(cov ~ log(mean_expr), data = stats)
        stats$pred_cov_feats <- stats::predict(loess_model, newdata = stats)
        stats[, cov_diff := cov - pred_cov_feats]
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
        if (is.null(x@params$norm) ||
            is.null(x@params$norm$scale_factors)) {
            stop("[analyzeData(parquetExprStore, covGroupsParam)] ",
                 "expression backend has no normalization recipe. Run ",
                 "normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE) ",
                 "first to populate scale factors on the store.",
                 call. = FALSE)
        }

        thr <- param$detection_threshold %null% 0
        stats <- .stream_norm_gene_stats(x, expression_threshold = thr)

        # NSE bindings
        cov <- expr_groups <- cov_group_zscore <- NULL

        # Quantile-bin by mean expression. If too many tied breaks (lots of
        # zero-mean genes), recompute on the strictly positive subset and
        # set the leading break to 0 so all-zero genes still bin into
        # group_1 — matches Giotto's in-memory .calc_cov_group_hvf.
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

    norm <- pe@params$norm %null% list()
    sf   <- norm$scale_factors
    if (is.null(sf))
        stop("[.stream_norm_gene_stats] no scale_factors on pe@params$norm.")
    log_norm <- isTRUE(norm$log)
    log_base <- norm$base %null% 2
    thr      <- as.numeric(expression_threshold)

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)

    gene_sum   <- numeric(n_genes)
    gene_sumsq <- numeric(n_genes)
    gene_nnz   <- integer(n_genes)
    gene_total_raw <- numeric(n_genes)

    # NSE bindings
    row_id <- col_id <- value <- v_norm <- s <- s2 <- nz <- raw_total <- NULL

    # storeRead already filters by pe@cell_idx / pe@gene_idx if set, but we
    # need to walk by SUBSET positions and translate to original parquet
    # row_ids when the store has been subsetted.
    ds_base <- pe@read_fun(pe@path)   # unfiltered — we apply our own
    chunk_size <- as.integer(pe@chunk_size %null% 250000L)

    for (start in seq.int(1L, n_cells, by = chunk_size)) {
        end <- min(start + chunk_size - 1L, n_cells)
        # Translate subset positions [start..end] to original parquet row_ids
        orig_rows <- .pe_orig_row(start:end, pe)

        q <- ds_base |> dplyr::filter(row_id %in% !!orig_rows)
        if (length(pe@gene_idx) > 0L) {
            gi <- pe@gene_idx
            q <- q |> dplyr::filter(col_id %in% !!gi)
        }
        chunk <- q |> dplyr::collect() |> data.table::as.data.table()
        if (nrow(chunk) == 0L) next

        # Remap row_id from original parquet -> subset position so
        # sf[row_id] works (sf is in subset coords)
        chunk[, row_id := .pe_remap_row(row_id, pe)]

        # Apply JIT recipe: scale per cell, then optional log
        chunk[, v_norm := value * sf[row_id]]
        if (log_norm) chunk[, v_norm := log1p(v_norm) / log(log_base)]

        # Per-gene accumulators (use ORIGINAL col_id for grouping then
        # remap once at the end)
        agg <- chunk[, .(
                s         = sum(v_norm),
                s2        = sum(v_norm * v_norm),
                nz        = sum(v_norm > thr),
                raw_total = sum(value)
            ), by = col_id]

        idx <- .pe_remap_col(agg$col_id, pe)
        keep <- !is.na(idx)
        gene_sum[idx[keep]]       <- gene_sum[idx[keep]]       + agg$s[keep]
        gene_sumsq[idx[keep]]     <- gene_sumsq[idx[keep]]     + agg$s2[keep]
        gene_nnz[idx[keep]]       <- gene_nnz[idx[keep]]       + as.integer(agg$nz[keep])
        gene_total_raw[idx[keep]] <- gene_total_raw[idx[keep]] + agg$raw_total[keep]
    }

    gene_mean <- gene_sum / n_cells
    gene_var  <- pmax(gene_sumsq / n_cells - gene_mean * gene_mean, 0)
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
