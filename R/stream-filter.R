#' @include class-parquetExprStore.R
NULL

# stream-filter ####
# Streaming mask computation for parquetExprBase-backed expression.
# Plugs into GiottoClass's filterData(x, param) dispatch via:
#
#   filterData(parquetExprBase, filterParam)
#       -> list(feats_keep = <character>, cells_keep = <character>)
#
# Implements Giotto's two-stage filter in two streaming Arrow passes:
#   1. feature mask: count cells per gene with value >= threshold
#   2. cell mask:    count detected genes per cell, restricted to kept genes
#
# Single (parquetExprStore) and union (unionParquetExprStore) collapse to
# one algorithm via the `.exprbase_substores()` iterator: feat counts sum
# across substores (feat_ids align by union invariant) and per-substore
# cell counts concatenate (substore cell_ids are disjoint).

#' @rdname filterData
#' @export
setMethod("filterData",
    signature(x = "parquetExprBase", param = "filterParam"),
    function(x, param, ...) {
        do.call(.stream_filter_masks,
            c(list(pe = x), as.list(param@param)))
    }
)


.stream_filter_masks <- function(pe,
                                  expression_threshold,
                                  feat_det_in_min_cells,
                                  min_det_feats_per_cell, ...) {
    if (!inherits(pe, "parquetExprBase"))
        stop("[.stream_filter_masks] pe must be a parquetExprBase.")

    thr   <- as.numeric(expression_threshold)
    f_min <- as.integer(feat_det_in_min_cells)
    c_min <- as.integer(min_det_feats_per_cell)

    # NSE bindings (silence R CMD check)
    row_id <- col_id <- value <- n_cells_above <- n_feats_above <- NULL

    subs <- .exprbase_substores(pe)
    n_feats <- as.integer(pe@n_genes)
    feat_count <- integer(n_feats)

    # ---- Pass 1: per-feature count, summed across all substores -----------
    # feat_ids align across substores (union invariant), so each substore's
    # `.pe_remap_col` yields consistent local indices into the union feat
    # axis. For a single parquetExprStore the loop runs once.
    for (sub in subs) {
        store <- sub$store
        agg <- storeRead(store) |>
            dplyr::filter(value >= !!thr) |>
            dplyr::group_by(col_id) |>
            dplyr::summarise(n_cells_above = dplyr::n()) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(agg) > 0L) {
            idx <- .pe_remap_col(agg$col_id, store)
            keep <- !is.na(idx)
            feat_count[idx[keep]] <- feat_count[idx[keep]] +
                as.integer(agg$n_cells_above[keep])
        }
    }
    feats_keep_idx_subset <- which(feat_count >= f_min)
    feats_keep_ids <- pe@feat_ids[feats_keep_idx_subset]

    # ---- Pass 2: per-cell count per substore, restricted to kept feats ----
    if (length(feats_keep_idx_subset) == 0L) {
        return(list(feats_keep = character(0L),
                    cells_keep = character(0L)))
    }
    # Substore cell_ids are disjoint (union constructor invariant for
    # multi-substore; trivially disjoint for single), so concatenating
    # the per-substore kept ids reproduces the union's surviving set in
    # substore order — matches the union's `@cell_ids` concatenation.
    cells_keep_ids <- character(0L)
    for (sub in subs) {
        store <- sub$store
        # Translate kept-feat subset positions to ORIGINAL parquet col_ids
        # for the Arrow filter — payload uses on-disk indices regardless
        # of substore subsetting.
        feats_keep_orig <- .pe_orig_col(feats_keep_idx_subset, store)
        agg <- storeRead(store) |>
            dplyr::filter(value >= !!thr,
                          col_id %in% !!feats_keep_orig) |>
            dplyr::group_by(row_id) |>
            dplyr::summarise(n_feats_above = dplyr::n()) |>
            dplyr::collect() |>
            data.table::as.data.table()
        n_sub_cells <- as.integer(store@n_cells)
        cell_count <- integer(n_sub_cells)
        if (nrow(agg) > 0L) {
            idx <- .pe_remap_row(agg$row_id, store)
            keep <- !is.na(idx)
            cell_count[idx[keep]] <- as.integer(agg$n_feats_above[keep])
        }
        cells_keep_ids <- c(cells_keep_ids,
            store@cell_ids[cell_count >= c_min])
    }

    list(feats_keep = feats_keep_ids,
         cells_keep = cells_keep_ids)
}
