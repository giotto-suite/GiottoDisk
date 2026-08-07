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

    # Filtering reads the STORED values, ignoring anything queued: Giotto's
    # filterGiotto takes counts, and its thresholds are count thresholds. The
    # chain is suppressed explicitly rather than relied on to be absent, since
    # nothing stops a caller from filtering after normalizing.
    #
    # `inclusive = TRUE`: Giotto's `expression_threshold = 1` means "expressed
    # if the count is at least 1", unlike the statistic verbs' strict `>`.
    base <- .pe_chain_none(pe)

    # ---- Pass 1: per-feature detection count ------------------------------
    # One accumulator pass covers a union in a single Acero plan and resolves
    # the feature axis itself, replacing the per-substore loop and its manual
    # `.pe_remap_col` accumulate.
    feat_count <- .stream_expr_accum(base, axis = "feat", thr = thr,
                                     stats = "nnz", inclusive = TRUE)$nnz
    feats_keep_idx_subset <- which(feat_count >= f_min)
    feats_keep_ids <- pe@feat_ids[feats_keep_idx_subset]
    if (length(feats_keep_idx_subset) == 0L) {
        return(list(feats_keep = character(0L),
                    cells_keep = character(0L)))
    }

    # ---- Pass 2: per-cell detection count over the kept features ----------
    # The restriction is a gene subset, so `[` carries it: storeRead turns
    # @gene_idx into the same predicate the explicit `col_id %in% ...` filter
    # used to build, and the accumulator's cell axis is already the union's.
    cell_count <- .stream_expr_accum(base[feats_keep_idx_subset, ],
                                     axis = "cell", thr = thr,
                                     stats = "nnz", inclusive = TRUE)$nnz
    cells_keep_ids <- pe@cell_ids[cell_count >= c_min]

    list(feats_keep = feats_keep_ids,
         cells_keep = cells_keep_ids)
}
