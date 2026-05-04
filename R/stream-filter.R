#' @include class-parquetExprStore.R
NULL

# stream-filter ####
# Streaming mask computation for parquetExprStore-backed expression.
# Plugs into Giotto's existing processData(x, param) dispatch via:
#
#   processData(parquetExprStore, filterParam)
#       -> list(feats_keep = <character>, cells_keep = <character>)
#
# Implements Giotto's two-stage filter in two streaming Arrow passes:
#   1. feature mask: count cells per gene with value >= threshold
#   2. cell mask:    count detected genes per cell, restricted to kept genes
#
# This matches Giotto's in-memory two-stage convention exactly (and fixes
# the bug noted in project.md where the original scstream sc_filter used
# pre-computed n_genes_detected without re-counting after the gene mask).

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprStore", param = "filterParam"),
    function(x, param, ...) {
        thr   <- param$expression_threshold
        f_min <- param$feat_det_in_min_cells
        c_min <- param$min_det_feats_per_cell

        .stream_filter_masks(
            x,
            expression_threshold   = thr,
            feat_det_in_min_cells  = f_min,
            min_det_feats_per_cell = c_min
        )
    }
)


.stream_filter_masks <- function(pe,
                                  expression_threshold,
                                  feat_det_in_min_cells,
                                  min_det_feats_per_cell) {
    if (!inherits(pe, "parquetExprStore"))
        stop("[.stream_filter_masks] pe must be a parquetExprStore.")

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)
    thr     <- as.numeric(expression_threshold)
    f_min   <- as.integer(feat_det_in_min_cells)
    c_min   <- as.integer(min_det_feats_per_cell)

    # NSE bindings (silence R CMD check)
    row_id <- col_id <- value <- n_cells_above <- n_feats_above <- NULL

    ds <- storeRead(pe)

    # ---- Pass 1: per-feature count (cells where value >= threshold) -------
    feat_agg <- ds |>
        dplyr::filter(value >= !!thr) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(n_cells_above = dplyr::n()) |>
        dplyr::collect() |>
        data.table::as.data.table()

    feat_count <- integer(n_genes)
    if (nrow(feat_agg) > 0L) {
        idx <- .pe_remap_col(feat_agg$col_id, pe)
        keep <- !is.na(idx)
        feat_count[idx[keep]] <- as.integer(feat_agg$n_cells_above[keep])
    }
    feats_keep_idx_subset <- which(feat_count >= f_min)
    feats_keep_ids <- pe@feat_ids[feats_keep_idx_subset]

    # ---- Pass 2: per-cell count, restricted to kept features --------------
    if (length(feats_keep_idx_subset) == 0L) {
        # nothing kept — short-circuit
        return(list(feats_keep = character(0L),
                    cells_keep = character(0L)))
    }

    # Translate kept-feature subset positions back to ORIGINAL parquet
    # col_ids for the Arrow filter (the parquet payload uses on-disk
    # original indices regardless of whether pe is subsetted).
    feats_keep_orig <- .pe_orig_col(feats_keep_idx_subset, pe)

    cell_agg <- ds |>
        dplyr::filter(value >= !!thr,
                       col_id %in% !!feats_keep_orig) |>
        dplyr::group_by(row_id) |>
        dplyr::summarise(n_feats_above = dplyr::n()) |>
        dplyr::collect() |>
        data.table::as.data.table()

    cell_count <- integer(n_cells)
    if (nrow(cell_agg) > 0L) {
        idx <- .pe_remap_row(cell_agg$row_id, pe)
        keep <- !is.na(idx)
        cell_count[idx[keep]] <- as.integer(cell_agg$n_feats_above[keep])
    }
    cells_keep_idx <- which(cell_count >= c_min)
    cells_keep_ids <- pe@cell_ids[cells_keep_idx]

    list(feats_keep = feats_keep_ids,
         cells_keep = cells_keep_ids)
}
