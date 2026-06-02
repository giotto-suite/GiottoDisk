#' @include class-parquetEdgeStore.R
#' @include class-parquetStore.R
NULL

# ============================================================================
# analyzeData(parquetEdgeStore, labelProportionsParam)
#
# Lazy dplyr/arrow pipeline that mirrors the in-mem
# analyzeData(igraph, labelProportionsParam) method in GiottoClass:
#
#   1. Symmetrize edges via dplyr::union_all over the same lazy dataset —
#      two scans, no materialization of the doubled edge set.
#   2. Add self edges from the nodes sidecar (gated by alpha).
#   3. Translate the labels DT (cell_ID, label) to int-ID space via a
#      sidecar filter (small).
#   4. Join + group_by + summarise lazily in arrow.
#   5. Collect (small: groups × labels), back-translate group integer IDs
#      to character cell_IDs via sidecar, compute proportions, dcast.
#
# Carries a unified weight column `w` throughout:
#  - native edge weight if @weight column exists and param$weights = TRUE
#  - synthetic 1.0 otherwise (warning if user requested weights and the
#    edge store has no weight column)
#  - self edges contribute weight = alpha
# Aggregation is always sum(w). Pure adjacency emerges as the
# (alpha = 1, weights = FALSE) special case where all w values are 1.
# ============================================================================

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetEdgeStore",
              param = "labelProportionsParam"),
    function(x, param, ..., labels = NULL) {
        if (is.null(labels)) {
            stop("[analyzeData(parquetEdgeStore, labelProportionsParam)] ",
                 "`labels` data.table is required (cell_ID + label column)",
                 call. = FALSE)
        }
        labels_col <- param$labels
        alpha      <- param$alpha
        weights    <- param$weights

        # NSE vars
        from_id <- to_id <- int_id <- node_id <- weight <- NULL
        cell_ID_int <- group <- label <- w <- n_LPG <- n_NPG <- NULL

        edges <- storeRead(x, output = "arrow")
        use_native_weight <- "weight" %in% names(edges) && isTRUE(weights)
        if (isTRUE(weights) && !use_native_weight) {
            warning("[analyzeData(parquetEdgeStore, labelProportionsParam)] ",
                    "No 'weight' column on edge store; falling back to ",
                    "adjacency.", call. = FALSE)
        }

        # Symmetrize edges (lazy union of two scans).
        if (use_native_weight) {
            forward  <- edges |>
                dplyr::transmute(group = from_id, cell_ID_int = to_id, w = weight)
            backward <- edges |>
                dplyr::transmute(group = to_id,   cell_ID_int = from_id, w = weight)
        } else {
            forward  <- edges |>
                dplyr::transmute(group = from_id, cell_ID_int = to_id, w = 1.0)
            backward <- edges |>
                dplyr::transmute(group = to_id,   cell_ID_int = from_id, w = 1.0)
        }
        sym_edges <- dplyr::union_all(forward, backward)

        # Self edges, gated by alpha.
        if (alpha != 0) {
            self_q <- storeRead(x@nodes) |>
                dplyr::transmute(group = int_id, cell_ID_int = int_id,
                                  w = as.double(alpha))
            sym_edges <- dplyr::union_all(sym_edges, self_q)
        }

        # Translate labels DT to int-ID space via sidecar filter (small).
        # labels: cell_ID + labels_col character columns.
        labels_dt <- data.table::as.data.table(labels)
        used_ids <- unique(labels_dt$cell_ID)
        name_map <- storeRead(x@nodes) |>
            dplyr::filter(node_id %in% !!used_ids) |>
            dplyr::select(int_id, node_id) |>
            dplyr::collect() |>
            data.table::as.data.table()
        labels_int <- merge(labels_dt, name_map,
            by.x = "cell_ID", by.y = "node_id")
        if (nrow(labels_int) == 0L) {
            stop("[analyzeData(parquetEdgeStore, labelProportionsParam)] ",
                 "labels DT had zero overlap with the network's node sidecar.",
                 call. = FALSE)
        }
        labels_arrow <- arrow::arrow_table(
            cell_ID_int = labels_int$int_id,
            label       = labels_int[[labels_col]]
        )

        # Join + aggregate, all lazy in arrow.
        joined <- sym_edges |>
            dplyr::left_join(labels_arrow, by = "cell_ID_int")
        counts_q <- joined |>
            dplyr::group_by(group, label) |>
            dplyr::summarise(n_LPG = sum(w, na.rm = TRUE), .groups = "drop")

        # Collect (small: groups × labels), finish in R.
        counts_dt <- data.table::as.data.table(dplyr::collect(counts_q))
        if (nrow(counts_dt) == 0L) {
            stop("[analyzeData(parquetEdgeStore, labelProportionsParam)] ",
                 "no rows after aggregation.", call. = FALSE)
        }

        # Back-translate group integer IDs to character cell_IDs.
        used_groups <- unique(counts_dt$group)
        group_map <- storeRead(x@nodes) |>
            dplyr::filter(int_id %in% !!used_groups) |>
            dplyr::select(int_id, node_id) |>
            dplyr::collect() |>
            data.table::as.data.table()
        counts_dt <- merge(counts_dt, group_map,
            by.x = "group", by.y = "int_id")
        counts_dt[, group := node_id]
        counts_dt[, node_id := NULL]

        # Compute proportions.
        totals <- counts_dt[, .(n_NPG = sum(n_LPG)), by = group]
        counts_dt <- merge(counts_dt, totals, by = "group")
        counts_dt[, "prop" := n_LPG / n_NPG]

        # Rename label → labels_col for dcast compatibility with the
        # in-mem branch (its wide output uses the label column name).
        data.table::setnames(counts_dt, "label", labels_col)

        data.table::dcast(counts_dt,
            formula = paste("group", labels_col, sep = "~"),
            fill = 0,
            value.var = "prop"
        )
    }
)


# ============================================================================
# analyzeData(parquetGeomBase, labelProportionsParam) — stub
#
# Disk-backed polygon dispatch is now reachable: the GiottoClass giotto-
# class router for group_method = "polygon" delegates to
# analyzeData(<polygon>, labelProportionsParam), and a disk-backed
# gobject hands us a parquetGeomBase instead of an in-mem giottoPolygon.
#
# Real implementation requires a "spatial join" form — pairs of
# (polygon_id, cell_id) that intersect — which the current
# spatRelate(parquetGeomBase, ...) methods don't return (they return a
# narrowed x via filter semantics, not pairs). Once spatRelate gains a
# `form = "join"` mode (the planned but not-yet-implemented case in
# methods-spatRelate.R) this stub becomes a real implementation:
#
#   pairs <- spatRelate(x, y_points, relation = "intersects",
#                       form = "join") |>
#       storeRead(output = "tibble") |>
#       data.table::as.data.table()
#   .lp_aggregate(pairs, labels, param$labels, ...)
#
# Until then: error early with a clear pointer.
# ============================================================================

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetGeomBase",
              param = "labelProportionsParam"),
    function(x, param, ..., labels = NULL, y = NULL) {
        stop("[analyzeData(parquetGeomBase, labelProportionsParam)] ",
             "polygon group_method on disk-backed geometries needs ",
             "spatRelate(form = \"join\"), which is not yet implemented. ",
             "Workaround: materialize the polygons (as.terra() on the ",
             "store, wrap in createGiottoPolygon) and use the in-mem ",
             "analyzeData(giottoPolygon, labelProportionsParam) method.",
             call. = FALSE)
    }
)
