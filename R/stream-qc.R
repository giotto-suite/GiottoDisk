#' @include class-parquetExprStore.R
NULL

# stream-qc ####
# Streaming QC stats for parquetExprBase-backed expression. Dispatches
# via Giotto's analyzeData(x, param) generic — stats are computed and
# returned; x is not mutated (selection / filtering is a separate step
# under processData).
#
#   analyzeData(parquetExprBase, cellStatsParam) -> per-cell stats data.table
#   analyzeData(parquetExprBase, featStatsParam) -> per-feature stats data.table
#
# Single (`parquetExprStore`) and union (`unionParquetExprStore`) collapse
# to one implementation via `.exprbase_substores()`. Cell stats: per-
# substore per-cell aggregations concatenate in substore order (matches
# the union's `@cell_ids` ordering). Feat stats: per-substore per-feat
# aggregations sum across substores (feat axis aligns by union invariant).

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "cellStatsParam"),
    function(x, param, ...) {
        do.call(.stream_cell_qc_stats,
            c(list(pe = x), as.list(param@param)))
    }
)

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "featStatsParam"),
    function(x, param, ...) {
        do.call(.stream_feat_qc_stats,
            c(list(pe = x), as.list(param@param)))
    }
)


# Internal streaming workers ####
# Single Arrow scan per axis (per substore). Result columns and ordering
# exactly match the in-memory addCellStatistics / addFeatStatistics output,
# so the downstream addCellMetadata / addFeatMetadata calls in those
# functions work without any additional adapter logic.

.stream_cell_qc_stats <- function(pe, detection_threshold = 0, ...) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[.stream_cell_qc_stats] pe must be a parquetExprBase.")

    n_genes <- as.integer(pe@n_genes)
    thr     <- as.numeric(detection_threshold)

    # NSE bindings
    row_id <- col_id <- value <- total_expr <- nr_feats <- NULL

    # Per-substore loop. Cell axes are disjoint across substores (union
    # constructor invariant), so concatenating per-substore result rows
    # in substore order reconstructs the union's cell-axis ordering.
    cell_ids_out <- character(0L)
    cell_total   <- numeric(0L)
    cell_nfeats  <- integer(0L)
    for (sub_entry in .exprbase_substores(pe)) {
        sub <- sub_entry$store
        agg <- storeRead(sub) |>
            dplyr::filter(value > !!thr) |>
            dplyr::group_by(row_id) |>
            dplyr::summarise(
                total_expr = sum(value, na.rm = TRUE),
                nr_feats   = dplyr::n()
            ) |>
            dplyr::collect() |>
            data.table::as.data.table()
        n_sub <- as.integer(sub@n_cells)
        sub_total  <- numeric(n_sub)
        sub_nfeats <- integer(n_sub)
        if (nrow(agg) > 0L) {
            idx <- .pe_remap_row(agg$row_id, sub)
            keep <- !is.na(idx)
            sub_total[idx[keep]]  <- as.numeric(agg$total_expr[keep])
            sub_nfeats[idx[keep]] <- as.integer(agg$nr_feats[keep])
        }
        cell_ids_out <- c(cell_ids_out, sub@cell_ids)
        cell_total   <- c(cell_total, sub_total)
        cell_nfeats  <- c(cell_nfeats, sub_nfeats)
    }

    data.table::data.table(
        cells      = cell_ids_out,
        nr_feats   = cell_nfeats,
        perc_feats = cell_nfeats / n_genes * 100,
        total_expr = cell_total
    )
}


.stream_feat_qc_stats <- function(pe, detection_threshold = 0, ...) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[.stream_feat_qc_stats] pe must be a parquetExprBase.")

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)
    thr     <- as.numeric(detection_threshold)

    # NSE bindings
    col_id <- value <- total_expr <- nr_cells <- NULL

    feat_total  <- numeric(n_genes)
    feat_ncells <- integer(n_genes)

    # Per-substore loop. Feat axis aligns across substores (union invariant
    # — identical `@feat_ids`), so per-substore aggregates sum into the
    # same feat positions.
    for (sub_entry in .exprbase_substores(pe)) {
        sub <- sub_entry$store
        agg <- storeRead(sub) |>
            dplyr::filter(value > !!thr) |>
            dplyr::group_by(col_id) |>
            dplyr::summarise(
                total_expr = sum(value, na.rm = TRUE),
                nr_cells   = dplyr::n()
            ) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(agg) > 0L) {
            idx <- .pe_remap_col(agg$col_id, sub)
            keep <- !is.na(idx)
            feat_total[idx[keep]]  <- feat_total[idx[keep]] +
                as.numeric(agg$total_expr[keep])
            feat_ncells[idx[keep]] <- feat_ncells[idx[keep]] +
                as.integer(agg$nr_cells[keep])
        }
    }

    data.table::data.table(
        feats         = pe@feat_ids,
        nr_cells      = feat_ncells,
        perc_cells    = feat_ncells / n_cells * 100,
        total_expr    = feat_total,
        mean_expr     = feat_total / n_cells,
        mean_expr_det = ifelse(feat_ncells > 0L,
                                feat_total / feat_ncells, NaN)
    )
}
