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
# Implementation: per-substore R-side pass. Arrow query yields raw
# triplets (arrow-side @ops applied, if any). @post_ops (norm) is applied
# R-side after collect. Per-gene stats reduce via data.table by-group.
# LOESS / bin-zscore run on the n_genes-sized aggregate vectors.
#
# varParam still errors clearly because per-gene variance on a scaled
# (z-scored) matrix requires materialising the dense matrix.

# ---- covLoessParam: streaming ---------------------------------------------

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "covLoessParam"),
    function(x, param, ...) {
        if (!.pe_has_norm_op(x)) {
            stop("[analyzeData(parquetExprBase, covLoessParam)] ",
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
    signature(x = "parquetExprBase", param = "covGroupsParam"),
    function(x, param, ...) {
        if (!.pe_has_norm_op(x)) {
            stop("[analyzeData(parquetExprBase, covGroupsParam)] ",
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


# ---- Internal: streaming per-gene stats with JIT normalization ------------
#
# Two paths:
#   (a) Band-parallel via `lapply_flex` when a parallel `future` plan is
#       set AND the store's op state is band-worker-safe (no arrow-side
#       @ops; @post_ops limited to norm_libsize_log which the worker
#       inlines). Each worker handles a contiguous cell range, returns
#       partial per-gene sums, main thread reduces.
#   (b) Serial single-collect otherwise (existing path, unchanged).
#
# Worker contract for (a):
#   * Runs in a fresh R session (mirai_multisession) or in-process
#     (sequential fallback). Cannot see GiottoDisk internal S4 methods —
#     re-opens arrow dataset directly on the store path and computes
#     the norm math inline via a positional scalef vector.
#   * @ops must be empty and @post_ops must contain only known types,
#     verified before dispatch.

# Plain-R worker.  Takes only serializable primitives + vectors — no S4
# objects, no method dispatch — so future.mirai can ship it to fresh
# worker sessions with just arrow/dplyr/data.table loaded.
#
# Shape mirrors scstream::.read_chunk_portable + .hvg_band_w:
#   * Range predicate at the arrow scan (`row_id >= cs, row_id <= ce`)
#     — pushes down cleanly to parquet row-group stats.
#   * If cell_idx has gaps, `findInterval` gives the position window in
#     cell_idx that falls in the raw range; in-mem `%in%` filter drops
#     the dropped cells.
#   * Build a chunk_n × n_genes sparseMatrix from the surviving triplets.
#   * Apply norm in place: `A@x <- A@x * sf[A@i + 1L]`, then log1p — no
#     column allocations, no data.table update-join.
#   * `Matrix::colSums(A)` for normalized sum; `rowsum(A@x^2, col_of_x)`
#     for sumsq; `tabulate(col_of_x[A@x > thr])` for nnz-above-threshold.
.hvg_band_worker <- function(band_cells,
                              sub_path, sub_uid,
                              cell_idx, gene_idx,
                              n_genes,
                              scalef_vec,
                              log_flag, log_base,
                              thr) {
    row_id <- col_id <- source_id <- NULL   # NSE

    empty <- list(
        s = numeric(n_genes), s2 = numeric(n_genes),
        nz = integer(n_genes)
    )

    # Raw row_id range covering the band's cells.  When cell_idx is
    # empty, band positions ARE raw row_ids; otherwise cell_idx maps
    # position → raw row_id.  Because cell_idx is monotonic ascending
    # (populated via `which()` in filterData), the band's first / last
    # positions bracket the full raw range.
    cs_raw <- if (length(cell_idx) > 0L)
        as.integer(cell_idx[band_cells[1L]])
    else
        as.integer(band_cells[1L])
    ce_raw <- if (length(cell_idx) > 0L)
        as.integer(cell_idx[band_cells[length(band_cells)]])
    else
        as.integer(band_cells[length(band_cells)])

    # Range predicate — pushed down to parquet row-group stats.
    ds <- arrow::open_dataset(sub_path)
    q <- ds |>
        dplyr::filter(source_id == !!sub_uid,
                       row_id    >= !!cs_raw,
                       row_id    <= !!ce_raw)
    if (length(gene_idx) > 0L) {
        q <- q |> dplyr::filter(col_id %in% !!gene_idx)
    }
    df <- dplyr::collect(q) |> data.table::as.data.table()
    if (nrow(df) == 0L) return(empty)

    # If cell_idx has gaps in this range, drop rows outside the active
    # window.  cell_idx[band_cells] gives the exact active raw row_ids;
    # findInterval finds the position window in cell_idx that overlaps
    # [cs_raw, ce_raw].
    if (length(cell_idx) > 0L) {
        lo <- findInterval(cs_raw - 1L, cell_idx) + 1L
        hi <- findInterval(ce_raw,      cell_idx)
        if (lo > hi) return(empty)
        active_raw <- as.integer(cell_idx[lo:hi])
        df <- df[row_id %in% active_raw]
        if (nrow(df) == 0L) return(empty)
    }

    # Build a chunk × n_genes sparseMatrix.
    if (length(cell_idx) > 0L) {
        chunk_cells <- sort(unique(df$row_id))
        chunk_n     <- length(chunk_cells)
        cell_map    <- match(df$row_id, chunk_cells)
    } else {
        chunk_cells <- cs_raw:ce_raw
        chunk_n     <- ce_raw - cs_raw + 1L
        cell_map    <- df$row_id - cs_raw + 1L
    }
    gene_map <- if (length(gene_idx) > 0L) {
        match(df$col_id, gene_idx)
    } else {
        as.integer(df$col_id)
    }
    A <- Matrix::sparseMatrix(
        i = cell_map, j = gene_map, x = as.double(df$value),
        dims = c(chunk_n, n_genes), repr = "C"
    )

    # Apply norm in place — matches scstream's `.hvg_band_w`.  scalef_vec
    # is positional in the substore's narrowed cell axis; for cell_idx
    # empty, chunk_cells IS the position, otherwise `match(chunk_cells,
    # cell_idx)` maps raw → narrowed position.
    if (length(scalef_vec) > 0L) {
        sf <- if (length(cell_idx) > 0L) {
            scalef_vec[match(chunk_cells, cell_idx)]
        } else {
            scalef_vec[chunk_cells]
        }
        A@x <- A@x * sf[A@i + 1L]
        if (log_flag) A@x <- log1p(A@x) / log(log_base)
    }

    s  <- as.numeric(Matrix::colSums(A))
    s2 <- numeric(n_genes)
    nz <- integer(n_genes)
    if (length(A@x) > 0L) {
        col_of_x <- rep.int(seq_len(n_genes), diff(A@p))
        csq <- rowsum(A@x * A@x, group = col_of_x, reorder = FALSE)
        grp <- as.integer(rownames(csq))
        s2[grp] <- as.numeric(csq[, 1L])
        if (isTRUE(thr == 0)) {
            nz <- as.integer(tabulate(col_of_x, nbins = n_genes))
        } else {
            nz_cols <- col_of_x[A@x > thr]
            if (length(nz_cols) > 0L) {
                nz <- as.integer(tabulate(nz_cols, nbins = n_genes))
            }
        }
    }

    list(s = s, s2 = s2, nz = nz)
}


.stream_norm_gene_stats <- function(pe, expression_threshold = 0) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[.stream_norm_gene_stats] pe must be a parquetExprBase.")
    if (!.pe_has_norm_op(pe))
        stop("[.stream_norm_gene_stats] pe has no norm op on @post_ops.")

    thr     <- as.numeric(expression_threshold)
    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)

    # Arrow-native eligibility.  The whole aggregate runs in Acero and only
    # the per-gene result (~18k rows) crosses into R, so it beats every
    # R-side chunking shape we measured -- on Atera (170k cells x 18k genes)
    # 1.94 s vs 9.59 s for 8-way band-parallel and 11.16 s via
    # storeRead(dgcmatrix), with no parallelism at all.  Needs the norm
    # expressed as arrow compute, so it is limited to a single
    # `norm_libsize_log` post-op; anything else falls through to the generic
    # R-side executor below, which handles any op kind.
    #
    # Union stores fall through too: the (source_id, orig_row_id) join would
    # span substores fine, but `.pe_remap_col` resolves against a single
    # store's @gene_idx.
    can_arrow <- inherits(pe, "parquetExprStore") &&
        length(pe@post_ops) == 1L &&
        identical(pe@post_ops[[1L]]$type, "norm_libsize_log")

    if (can_arrow) {
        return(.stream_norm_gene_stats_arrow(pe, thr,
            n_cells = n_cells, n_genes = n_genes))
    }

    # NSE bindings
    col_id <- value <- s <- s2 <- nz <- NULL

    gene_sum   <- numeric(n_genes)
    gene_sumsq <- numeric(n_genes)
    gene_nnz   <- integer(n_genes)

    # Iterate substores. Parent's @post_ops applies to each substore's
    # collected chunk (payloads carry source_id / orig_row_id so the
    # update-join resolves per substore).
    post_ops <- pe@post_ops
    subs <- .exprbase_substores(pe)
    for (sub_entry in subs) {
        sub <- sub_entry$store
        df <- storeRead(sub, output = "query") |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) next

        # No raw snapshot: `total_expr` follows Giotto, which reports
        # `rowSums()` of the values it was handed -- i.e. the normalized
        # sum, identical to `mean_expr * n_cells`.  Keeping a `raw_value`
        # copy would cost an extra numeric column over the whole triplet
        # stream (~2.5 GB at 307M rows) for a quantity nothing reads.
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
        total_expr = gene_sum,
        mean_expr  = gene_mean,
        sd         = gene_sd,
        cov        = gene_cov
    )
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
# `total_expr` follows Giotto, which reports `rowSums()` of the values it was
# handed -- the NORMALIZED sum, identical to `mean_expr * n_cells`. So it is
# just `s`; no separate raw aggregate is needed.

#' @keywords internal
#' @noRd
.stream_norm_gene_stats_arrow <- function(pe, thr, n_cells, n_genes) {
    # NSE bindings
    row_id <- col_id <- value <- nv <- scalef <- source_id <- NULL
    s <- s2 <- nz <- NULL

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

    q <- storeRead(pe, output = "query") |>
        dplyr::left_join(scalef_tbl,
            by = c("source_id" = "source_id", "row_id" = "orig_row_id")) |>
        dplyr::mutate(nv = value * scalef)
    if (log_flag) {
        q <- dplyr::mutate(q, nv = log1p(nv) / log(!!log_base))
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


# Band-parallel HVG stats.  Fork-based via `parallel::mclapply` on Unix
# (Linux/macOS): each fork COW-inherits the parent's loaded namespace and
# arrow buffers, so worker startup is ~zero cost and daemon-session
# memory blowup doesn't apply.  On Windows, `mclapply(mc.cores = n > 1)`
# degrades to sequential (documented `parallel` behavior) — HVG stays
# serial there until a fork alternative or a socket-worker path is added.
.stream_norm_gene_stats_bandparallel <- function(pe, thr,
                                                  n_cells, n_genes) {
    hvg_worker <- .hvg_band_worker

    # Only norm_libsize_log is supported here (checked by caller).
    norm_op  <- pe@post_ops[[1L]]
    log_flag <- isTRUE(norm_op$log)
    log_base <- norm_op$base %null% 2

    gene_sum       <- numeric(n_genes)
    gene_sumsq     <- numeric(n_genes)
    gene_nnz       <- integer(n_genes)

    n_workers <- .par_workers()

    for (sub_entry in .exprbase_substores(pe)) {
        sub         <- sub_entry$store
        n_sub       <- as.integer(sub@n_cells)
        scalef_vec  <- .pe_scalef_vecs_for_sub(pe@post_ops, sub@uid)[[1L]]

        # Contiguous cell-axis bands.  ceiling(n_sub / n_workers) per band;
        # last band may be short.
        band_size   <- max(1L, as.integer(ceiling(n_sub / n_workers)))
        band_starts <- seq.int(1L, n_sub, by = band_size)
        bands <- lapply(band_starts, function(s0)
            seq.int(s0, min(s0 + band_size - 1L, n_sub)))

        # Extract primitives so the worker closure doesn't capture sub@.
        sub_path <- sub@path
        sub_uid  <- sub@uid
        cell_idx <- if (length(sub@cell_idx) > 0L) sub@cell_idx else integer(0)
        gene_idx <- if (length(sub@gene_idx) > 0L) sub@gene_idx else integer(0)

        # Self-contained worker closure: pull the worker function's body
        # into a fresh environment parented to globalenv() so serialization
        # doesn't drag the GiottoDisk namespace along.  Pack all captured
        # args explicitly into a per-task list so future.apply doesn't do
        # implicit globals detection on the enclosing frame.  Analogous to
        # scstream's `environment(.read_chunk_portable) <- globalenv()`.
        worker_env <- new.env(parent = globalenv())
        worker_env$hvg_worker <- hvg_worker
        environment(worker_env$hvg_worker) <- worker_env
        task_fn <- function(t) {
            hvg_worker(
                band_cells = t$band_cells,
                sub_path   = t$sub_path,
                sub_uid    = t$sub_uid,
                cell_idx   = t$cell_idx,
                gene_idx   = t$gene_idx,
                n_genes    = t$n_genes,
                scalef_vec = t$scalef_vec,
                log_flag   = t$log_flag,
                log_base   = t$log_base,
                thr        = t$thr
            )
        }
        environment(task_fn) <- worker_env

        tasks <- lapply(bands, function(b) list(
            band_cells = b,
            sub_path   = sub_path,
            sub_uid    = sub_uid,
            cell_idx   = cell_idx,
            gene_idx   = gene_idx,
            n_genes    = n_genes,
            scalef_vec = scalef_vec,
            log_flag   = log_flag,
            log_base   = log_base,
            thr        = thr
        ))

        # Dispatch: skip auto-globals detection (we passed everything
        # in `tasks`), reuse persistent workers set up by the caller
        # via `future::plan()`.
        partials <- if (n_workers > 1L) {
            future.apply::future_lapply(tasks, task_fn,
                future.globals = FALSE,
                future.seed    = NULL)
        } else {
            lapply(tasks, task_fn)
        }

        # mclapply returns try-error objects on worker failure; surface them.
        errs <- vapply(partials, inherits, logical(1L), "try-error")
        if (any(errs)) {
            msg <- attr(partials[[which(errs)[1L]]], "condition")$message
            stop("[.stream_norm_gene_stats_bandparallel] worker failed: ",
                 msg, call. = FALSE)
        }

        for (p in partials) {
            gene_sum       <- gene_sum       + p$s
            gene_sumsq     <- gene_sumsq     + p$s2
            gene_nnz       <- gene_nnz       + p$nz
        }
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
        total_expr = gene_sum,
        mean_expr  = gene_mean,
        sd         = gene_sd,
        cov        = gene_cov
    )
}
