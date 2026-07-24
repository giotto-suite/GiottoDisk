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


# ---- varParam on parquet: clear error -------------------------------------

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "varParam"),
    function(x, param, ...) {
        stop("[analyzeData(parquetExprBase, varParam)] per-feature ",
             "variance on a scaled (z-scored) matrix requires ",
             "materialising the dense matrix and is not supported for ",
             "streaming backends. Use covLoessParam or covGroupsParam.",
             call. = FALSE)
    }
)


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
#   * Compute `total_expr` = `Matrix::colSums(A)` BEFORE norm mutation.
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
        nz = integer(n_genes), raw = numeric(n_genes)
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

    # `total_expr` = per-gene raw sum, computed BEFORE the in-place
    # normalization mutates A@x.
    raw <- as.numeric(Matrix::colSums(A))

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

    list(s = s, s2 = s2, nz = nz, raw = raw)
}


.stream_norm_gene_stats <- function(pe, expression_threshold = 0) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[.stream_norm_gene_stats] pe must be a parquetExprBase.")
    if (!.pe_has_norm_op(pe))
        stop("[.stream_norm_gene_stats] pe has no norm op on @post_ops.")

    thr     <- as.numeric(expression_threshold)
    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)

    # Band-parallel eligibility: workers > 1 (either via option or a
    # future plan), fork-safe platform, no arrow-side @ops, and @post_ops
    # limited to types the worker knows how to inline.
    can_parallel <- .par_workers() > 1L &&
        .Platform$OS.type == "unix" &&
        length(pe@ops) == 0L &&
        length(pe@post_ops) == 1L &&
        identical(pe@post_ops[[1L]]$type, "norm_libsize_log")

    if (can_parallel) {
        return(.stream_norm_gene_stats_bandparallel(pe, thr,
            n_cells = n_cells, n_genes = n_genes))
    }

    # NSE bindings
    col_id <- value <- raw_value <- s <- s2 <- nz <- raw_total <- NULL

    gene_sum       <- numeric(n_genes)
    gene_sumsq     <- numeric(n_genes)
    gene_nnz       <- integer(n_genes)
    gene_total_raw <- numeric(n_genes)

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

        # Snapshot raw before @post_ops mutates value
        df[, raw_value := value]
        df <- .pe_apply_post_ops_df(df, post_ops)

        # Per-gene aggregation R-side (data.table by-group)
        agg <- df[, .(
            s         = sum(value, na.rm = TRUE),
            s2        = sum(value * value, na.rm = TRUE),
            nz        = sum(value > thr, na.rm = TRUE),
            raw_total = sum(raw_value, na.rm = TRUE)
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
            gene_total_raw[idx[keep]] <- gene_total_raw[idx[keep]] +
                as.numeric(agg$raw_total[keep])
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
        total_expr = gene_total_raw,
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
    gene_total_raw <- numeric(n_genes)

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

        run_band <- function(band_cells) {
            hvg_worker(
                band_cells = band_cells,
                sub_path   = sub_path,
                sub_uid    = sub_uid,
                cell_idx   = cell_idx,
                gene_idx   = gene_idx,
                n_genes    = n_genes,
                scalef_vec = scalef_vec,
                log_flag   = log_flag,
                log_base   = log_base,
                thr        = thr
            )
        }

        # Fork on Unix; falls back to sequential lapply elsewhere. Note:
        # empirically the next arrow scan in the parent process takes
        # ~3-4 s longer than a fresh serial run — a mclapply cost we
        # haven't been able to reset via arrow::set_cpu_count(),
        # set_io_thread_count(), memory-pool release_unused, or gc(). The
        # HVG parallel save more than pays for it net.
        partials <- if (.Platform$OS.type == "unix" && n_workers > 1L) {
            parallel::mclapply(bands, run_band,
                mc.cores = n_workers,
                mc.preschedule = TRUE)
        } else {
            lapply(bands, run_band)
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
            gene_total_raw <- gene_total_raw + p$raw
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
        total_expr = gene_total_raw,
        mean_expr  = gene_mean,
        sd         = gene_sd,
        cov        = gene_cov
    )
}
