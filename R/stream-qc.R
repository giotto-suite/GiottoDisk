#' @include class-parquetExprStore.R
NULL

# stream-qc ####
# Streaming QC for parquetExprStore-backed expression. Plugs into Giotto's
# existing processData(x, param) dispatch via two setMethod calls:
#
#   processData(parquetExprStore, cellQcParam) -> per-cell stats data.table
#   processData(parquetExprStore, featQcParam) -> per-feature stats data.table
#
# The cellQcParam / featQcParam classes are defined in Giotto. Once both
# packages are attached, addStatistics(g) / addCellStatistics(g) /
# addFeatStatistics(g) auto-route to streaming code when @exprMat is a
# parquetExprStore — no separate user-facing function name is needed.

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprStore", param = "cellQcParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        .stream_cell_qc_stats(x, detection_threshold = thr)
    }
)

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprStore", param = "featQcParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        .stream_feat_qc_stats(x, detection_threshold = thr)
    }
)


# Internal streaming workers ####
# Single Arrow scan per axis. Result columns and ordering exactly match
# the in-memory addCellStatistics / addFeatStatistics output, so the
# downstream addCellMetadata / addFeatMetadata calls in those functions
# work without any additional adapter logic.

.stream_cell_qc_stats <- function(pe, detection_threshold = 0) {
    if (!inherits(pe, "parquetExprStore"))
        stop("[.stream_cell_qc_stats] pe must be a parquetExprStore.")

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)
    thr     <- as.numeric(detection_threshold)

    # NSE bindings (silence R CMD check)
    row_id <- col_id <- value <- total_expr <- nr_feats <- NULL

    ds <- storeRead(pe)
    cell_agg <- ds |>
        dplyr::filter(value > !!thr) |>
        dplyr::group_by(row_id) |>
        dplyr::summarise(
            total_expr = sum(value, na.rm = TRUE),
            nr_feats   = dplyr::n()
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    cell_total  <- numeric(n_cells)
    cell_nfeats <- integer(n_cells)
    if (nrow(cell_agg) > 0L) {
        idx <- .pe_remap_row(cell_agg$row_id, pe)
        keep <- !is.na(idx)
        cell_total[idx[keep]]  <- as.numeric(cell_agg$total_expr[keep])
        cell_nfeats[idx[keep]] <- as.integer(cell_agg$nr_feats[keep])
    }

    data.table::data.table(
        cells      = pe@cell_ids,
        nr_feats   = cell_nfeats,
        perc_feats = cell_nfeats / n_genes * 100,
        total_expr = cell_total
    )
}


.stream_feat_qc_stats <- function(pe, detection_threshold = 0) {
    if (!inherits(pe, "parquetExprStore"))
        stop("[.stream_feat_qc_stats] pe must be a parquetExprStore.")

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)
    thr     <- as.numeric(detection_threshold)

    row_id <- col_id <- value <- total_expr <- nr_cells <- NULL

    ds <- storeRead(pe)
    feat_agg <- ds |>
        dplyr::filter(value > !!thr) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(
            total_expr = sum(value, na.rm = TRUE),
            nr_cells   = dplyr::n()
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    feat_total  <- numeric(n_genes)
    feat_ncells <- integer(n_genes)
    if (nrow(feat_agg) > 0L) {
        idx <- .pe_remap_col(feat_agg$col_id, pe)
        keep <- !is.na(idx)
        feat_total[idx[keep]]  <- as.numeric(feat_agg$total_expr[keep])
        feat_ncells[idx[keep]] <- as.integer(feat_agg$nr_cells[keep])
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
