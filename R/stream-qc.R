#' @include class-parquetExprStore.R
NULL

# stream-qc ####

#' @name addStreamStatistics
#' @title Compute QC statistics on a parquet-backed expression matrix (streaming)
#' @description
#' Streaming-compatible analogue of [Giotto::addStatistics()] for giotto
#' objects whose expression backend is a [parquetExprStore-class].
#' Performs **one streaming pass** over the long-format Parquet via
#' Arrow `dplyr::summarise`, producing exactly the columns
#' `addCellStatistics()` and `addFeatStatistics()` produce on an
#' in-memory matrix:
#'
#' Per cell (added to `cellMetaObj`):
#' * `total_expr` — sum of expression across detected genes
#' * `nr_feats`   — number of detected genes (`value > detection_threshold`)
#' * `perc_feats` — `nr_feats / n_genes * 100`
#'
#' Per feature (added to `featMetaObj`):
#' * `total_expr`    — sum of expression across detected cells
#' * `nr_cells`      — number of detected cells
#' * `perc_cells`    — `nr_cells / n_cells * 100`
#' * `mean_expr`     — `total_expr / n_cells`  (mean over ALL cells)
#' * `mean_expr_det` — `total_expr / nr_cells` (mean over DETECTED cells, 0 if none)
#'
#' @param gobject A `giotto` object whose expression backend is a
#'   `parquetExprStore`.
#' @param spat_unit,feat_type Spatial unit / feature type to operate on
#'   (passed through to `getExpression` / `addCellMetadata` /
#'   `addFeatMetadata`).
#' @param detection_threshold A feature is "detected" in a cell if its
#'   value is strictly greater than this. Default `0`.
#' @param verbose Print step messages.
#' @return The `giotto` object with `cellMetaObj` + `featMetaObj` updated.
#' @export
addStreamStatistics <- function(gobject,
    spat_unit           = NULL,
    feat_type           = NULL,
    detection_threshold = 0,
    verbose             = TRUE
) {
    if (!inherits(gobject, "giotto")) {
        stop("[addStreamStatistics] gobject must be a `giotto` object.",
             call. = FALSE)
    }
    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- GiottoClass::set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    expr <- GiottoClass::getExpression(
        gobject   = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output    = "exprObj",
        set_defaults = FALSE
    )
    pe <- slot(expr, "exprMat")
    if (!inherits(pe, "parquetExprStore")) {
        stop("[addStreamStatistics] expression backend is ",
             class(pe)[1], "; this function requires a parquetExprStore.\n",
             "  For dgCMatrix / BPCells / etc. backends, use ",
             "Giotto::addStatistics().", call. = FALSE)
    }

    if (isTRUE(verbose)) message("[addStreamStatistics] streaming QC pass...")
    stats <- .stream_qc_pass(pe, detection_threshold = detection_threshold)

    # Drop any pre-existing stat columns to mirror addCellStatistics behaviour
    cm <- GiottoClass::getCellMetadata(
        gobject, spat_unit = spat_unit, feat_type = feat_type,
        output = "cellMetaObj", copy_obj = TRUE, set_defaults = FALSE
    )
    cm_names <- colnames(cm[])
    drop_c   <- intersect(c("nr_feats", "perc_feats", "total_expr"), cm_names)
    if (length(drop_c) > 0L) {
        cm[][, (drop_c) := NULL]
        gobject <- GiottoClass::setGiotto(gobject, cm, verbose = FALSE)
    }
    fm <- GiottoClass::getFeatureMetadata(
        gobject, spat_unit = spat_unit, feat_type = feat_type,
        output = "featMetaObj", copy_obj = TRUE, set_defaults = FALSE
    )
    fm_names <- colnames(fm[])
    drop_f <- intersect(
        c("nr_cells", "perc_cells", "total_expr", "mean_expr", "mean_expr_det"),
        fm_names
    )
    if (length(drop_f) > 0L) {
        fm[][, (drop_f) := NULL]
        gobject <- GiottoClass::setGiotto(gobject, fm, verbose = FALSE)
    }

    gobject <- GiottoClass::addCellMetadata(
        gobject       = gobject,
        spat_unit     = spat_unit,
        feat_type     = feat_type,
        new_metadata  = stats$cell_stats,
        by_column     = TRUE,
        column_cell_ID = "cells"
    )
    gobject <- GiottoClass::addFeatMetadata(
        gobject       = gobject,
        spat_unit     = spat_unit,
        feat_type     = feat_type,
        new_metadata  = stats$feat_stats,
        by_column     = TRUE,
        column_feat_ID = "feats"
    )

    if (isTRUE(verbose)) message("[addStreamStatistics] done.")
    gobject
}


# .stream_qc_pass ####
# Single Arrow streaming pass that computes per-cell and per-feature
# aggregates. Two grouped summarise() calls are issued (one per axis).
# Arrow handles the chunked scan transparently — peak RAM stays
# proportional to one row-group, not the full matrix.
.stream_qc_pass <- function(pe, detection_threshold = 0) {
    if (!inherits(pe, "parquetExprStore"))
        stop("[.stream_qc_pass] pe must be a parquetExprStore.")

    n_cells <- as.integer(pe@n_cells)
    n_genes <- as.integer(pe@n_genes)
    thr     <- as.numeric(detection_threshold)

    ds <- storeRead(pe)
    # NSE bindings for arrow / data.table — silence R CMD check
    row_id <- col_id <- value <- total_expr <- nr_feats <- nr_cells <- NULL

    # Per-cell aggregate
    cell_agg <- ds |>
        dplyr::filter(value > !!thr) |>
        dplyr::group_by(row_id) |>
        dplyr::summarise(
            total_expr = sum(value, na.rm = TRUE),
            nr_feats   = dplyr::n()
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    # Per-feature aggregate
    feat_agg <- ds |>
        dplyr::filter(value > !!thr) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(
            total_expr = sum(value, na.rm = TRUE),
            nr_cells   = dplyr::n()
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    # Pad to all cells / all genes (entries missing from parquet = empty)
    cell_total  <- numeric(n_cells)
    cell_nfeats <- integer(n_cells)
    if (nrow(cell_agg) > 0L) {
        idx <- as.integer(cell_agg$row_id)
        cell_total[idx]  <- as.numeric(cell_agg$total_expr)
        cell_nfeats[idx] <- as.integer(cell_agg$nr_feats)
    }
    feat_total  <- numeric(n_genes)
    feat_ncells <- integer(n_genes)
    if (nrow(feat_agg) > 0L) {
        idx <- as.integer(feat_agg$col_id)
        feat_total[idx]  <- as.numeric(feat_agg$total_expr)
        feat_ncells[idx] <- as.integer(feat_agg$nr_cells)
    }

    list(
        cell_stats = data.table::data.table(
            cells       = pe@cell_ids,
            nr_feats    = cell_nfeats,
            perc_feats  = cell_nfeats / n_genes * 100,
            total_expr  = cell_total
        ),
        feat_stats = data.table::data.table(
            feats         = pe@feat_ids,
            nr_cells      = feat_ncells,
            perc_cells    = feat_ncells / n_cells * 100,
            total_expr    = feat_total,
            mean_expr     = feat_total / n_cells,
            mean_expr_det = ifelse(feat_ncells > 0L,
                                    feat_total / feat_ncells, 0)
        )
    )
}
