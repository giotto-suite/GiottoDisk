#' @include class-viewCoordinator.R
NULL

# =============================================================================
# methods-resolveSubobject.R — parquetCoordinator dispatch
#
# Registers:
#   1. defaultViewCoordinator method on gsource so that gobjects with a
#      gsource-inheriting `@source` automatically select parquetCoordinator.
#   2. resolveSubobject methods on subobject classes whose internal data
#      slot may hold a parquetStore-inheriting backing.
#
# Cross-storage narrowing follows the heterogeneity-durable shape: a view
# filter's predicate columns may live on a different subobject than the
# target being narrowed, and either side can be backed (parquetBase) or
# in-mem (data.table). The flow is:
#
#   .find_store_with_cols(gobject, cols)
#       walks gobject's tabular subobjects and returns the first one whose
#       data covers `cols`. Tags it kind = "parquetBase" or "data.table".
#
#   .narrow_store_by_predicate(target_store, predicate, gobject)
#       three branches:
#         (a) predicate cols all on target  -> subset() the target store
#             directly (lazy filter op).
#         (b) owner is parquetBase          -> subset() the owner store and
#             queue a lazy join op via target_store[narrowed_owner, on = key,
#             nomatch = NULL]. No I/O at coordinator time.
#         (c) owner is in-mem data.table    -> eager-narrow the data.table,
#             project to the join key, arrow::arrow_table(), queue an
#             id_filter op (existing internal op type — semi-join in arrow,
#             EXISTS subquery in sedona).
# =============================================================================


# defaultViewCoordinator dispatch ####
# Registered against gsource so that any gobject whose @source inherits
# from gsource auto-selects parquetCoordinator. S4 inheritance picks up
# concrete gsource subclasses (gDirSource, etc.) automatically.

#' @rdname parquetCoordinator-class
#' @importFrom GiottoClass defaultViewCoordinator
#' @export
setMethod("defaultViewCoordinator", signature(source = "gsource"),
    function(source, ...) parquetCoordinator()
)


# Helpers ####

# Walk a gobject's tabular subobjects and return the first one whose data
# covers all of `cols`. Considered slots, in the same precedence order
# spatValues uses for column lookup:
#
#   1. cell metadata     (cellMetaObj@metaDT, key = cell_ID)
#   2. feat metadata     (featMetaObj@metaDT, key = feat_ID)
#   3. spatial locations (spatLocsObj@coordinates, key = cell_ID)
#   4. spatial enrichment(spatEnrObj@enrichDT, key = cell_ID)
#
# Matrix-shaped slots (expression, dim_reduction) are skipped — they need
# feat_ID-based routing rather than colname checks.
#
# Returns:
#   list(kind = "parquetBase" | "data.table", source = <store|dt>, key = <id>)
#   NULL if no slot covers all `cols`.
#
#' @keywords internal
#' @noRd
.find_store_with_cols <- function(gobject, cols,
    spat_unit = NULL, feat_type = NULL) {
    stopifnot(is.character(cols), length(cols) > 0L)

    .colnames_any <- function(src) {
        if (inherits(src, "parquetBase")) colnames(src)
        else if (data.table::is.data.table(src)) names(src)
        else character(0L)
    }
    .wrap <- function(src, key) {
        if (inherits(src, "parquetBase")) {
            list(kind = "parquetBase", source = src, key = key)
        } else if (data.table::is.data.table(src)) {
            list(kind = "data.table", source = src, key = key)
        } else {
            NULL
        }
    }
    .check <- function(src, key) {
        if (is.null(src)) return(NULL)
        if (!all(cols %in% .colnames_any(src))) return(NULL)
        .wrap(src, key)
    }

    cm <- tryCatch(GiottoClass::getCellMetadata(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type,
        output = "cellMetaObj", copy_obj = FALSE, set_defaults = TRUE
    ), error = function(e) NULL)
    if (!is.null(cm)) {
        hit <- .check(cm@metaDT, "cell_ID")
        if (!is.null(hit)) return(hit)
    }

    fm <- tryCatch(GiottoClass::getFeatureMetadata(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type,
        output = "featMetaObj", copy_obj = FALSE, set_defaults = TRUE
    ), error = function(e) NULL)
    if (!is.null(fm)) {
        hit <- .check(fm@metaDT, "feat_ID")
        if (!is.null(hit)) return(hit)
    }

    sl <- tryCatch(GiottoClass::getSpatialLocations(
        gobject = gobject, spat_unit = spat_unit,
        output = "spatLocsObj", copy_obj = FALSE, set_defaults = TRUE
    ), error = function(e) NULL)
    if (!is.null(sl)) {
        hit <- .check(sl@coordinates, "cell_ID")
        if (!is.null(hit)) return(hit)
    }

    se <- tryCatch(GiottoClass::getSpatialEnrichment(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type,
        output = "spatEnrObj", copy_obj = FALSE, set_defaults = TRUE
    ), error = function(e) NULL)
    if (!is.null(se)) {
        hit <- .check(se@enrichDT, "cell_ID")
        if (!is.null(hit)) return(hit)
    }

    # Expression — gene names as "columns".
    expr <- tryCatch(GiottoClass::getExpression(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type,
        output = "exprObj", set_defaults = TRUE
    ), error = function(e) NULL)

    if (!is.null(expr) && inherits(expr@exprMat, "parquetExprStore") &&
        all(cols %in% expr@exprMat@feat_ids)) {
        # parquetExprStore-specific path: read long-format triplets via
        # storeRead, pivot to wide. storeRead applies any pending @ops
        # (e.g. norm_libsize_log) inside the arrow query, so the long_dt
        # carries a projected `v_norm` column when normalization is
        # recorded — `.expr_store_gene_slice` picks v_norm if present,
        # else raw `value`.
        wide_dt <- .expr_store_gene_slice(expr@exprMat, cols)
        return(.wrap(wide_dt, "cell_ID"))
    }

    # Generic matrix-like fallback for in-mem and lazy backends:
    # plain matrix, dgCMatrix, DelayedArray, BPCells IterableMatrix —
    # anything that supports `m[feats, ]` and `as.matrix()`. Realization
    # composes whatever lazy ops the backend has queued
    # (DelayedArray DelayedOps; BPCells transforms), so the returned
    # values reflect prior normalize / log / scale steps.
    if (!is.null(expr)) {
        m <- expr@exprMat
        rn <- tryCatch(rownames(m), error = function(e) NULL)
        if (!is.null(rn) && length(rn) > 0L && all(cols %in% rn)) {
            m_slice <- m[cols, , drop = FALSE]
            # Realize to dense (genes × cells); transpose to cells × genes.
            dense <- as.matrix(m_slice)
            wide <- t(dense)
            wide_dt <- data.table::as.data.table(wide,
                keep.rownames = "cell_ID")
            return(.wrap(wide_dt, "cell_ID"))
        }
    }

    NULL
}


# Materialize a "gene slice" of a parquetExprStore as a wide
# data.table. Result has `cell_ID` plus one column per requested gene
# (values = expression for that gene, 0 for cells with no recorded
# value). Used by `.find_store_with_cols` to expose gene names as
# columns to `.compute_narrowed_ids` for predicate evaluation.
#
# Cost: reads only the rows of the long-format parquet whose `col_id`
# is in the requested gene set. For atlas-scale this is bounded by
# `n_cells × |genes_requested|` (typically 1M × few genes — tens of
# MB), much smaller than the full expression matrix.
#
#' @keywords internal
#' @noRd
.expr_store_gene_slice <- function(pe, feat_ids) {
    narrowed <- pe[feat_ids, , drop = FALSE]
    n_rows <- length(narrowed@feat_ids)

    # storeRead(output = "dgcmatrix") applies any pending @ops
    # (norm_libsize_log etc.) inside the arrow query and materializes
    # a gene x cell sparseMatrix with dimnames = (feat_ids, cell_ids).
    # Zero-expression cells are already represented via the matrix's
    # full cell_id dimnames — no backfill needed.
    #
    # Asymmetric guard: the feat axis is the narrow axis (the gene
    # slice the caller asked for); cell axis is the full population.
    # Override max_rows = n_rows so any requested feature count
    # passes the guard regardless of n_cells.
    m <- storeRead(narrowed, output = "dgcmatrix", max_rows = n_rows)

    # Pivot to cells x feats data.table. Densification is fine here:
    # the plot / predicate consumer needs every cell's value, and the
    # feat axis is bounded (typically a few dozen).
    data.table::data.table(
        cell_ID = colnames(m),
        as.matrix(Matrix::t(m))
    )
}


# Materialize the narrowed owner's id-column for a single predicate.
# Returns `list(ids_tab, key)`. Building block for the view-level
# intersection — not cached itself; the caller caches the final
# intersected result instead.
#
#' @keywords internal
#' @noRd
.compute_narrowed_ids <- function(predicate, gobject,
    spat_unit = NULL, feat_type = NULL) {
    cols <- all.vars(predicate)
    owner <- .find_store_with_cols(gobject, cols,
        spat_unit = spat_unit, feat_type = feat_type)
    if (is.null(owner)) {
        stop("[.compute_narrowed_ids] no subobject covers predicate ",
            "cols: ", paste(cols, collapse = ", "), call. = FALSE)
    }
    ids_df <- if (identical(owner$kind, "parquetBase")) {
        narrowed <- subset(owner$source, predicate, quote = FALSE)
        storeRead(narrowed[, owner$key, drop = FALSE], output = "tibble")
    } else {
        mask <- eval(predicate, envir = owner$source, enclos = baseenv())
        as.data.frame(unique(owner$source[mask, owner$key, with = FALSE]))
    }
    list(
        ids_tab = arrow::arrow_table(unique(ids_df)),
        key = owner$key
    )
}


# Resolve `cells in region` for one viewCrop step against the gobject's
# spatial locations, handling both single-giotto and giottoMulti shapes.
#
# * single giotto: `getSpatialLocations(g, output = "spatLocsObj")`
#   returns ONE spatLocsObj. Project through `space` once via
#   `.apply_space_to_subobj`, then check centroid-in-region.
#
# * giottoMulti: same getter returns a NAMED LIST of per-sample
#   spatLocsObjs. Per sample: scope the space to that sample (rename
#   key to `:default:` so `.apply_space_to_subobj` picks it up),
#   project, check centroid-in-region. Union surviving cell_IDs across
#   children — joint cellMeta has globally unique cell_IDs, so union
#   composes correctly.
#
#' @keywords internal
#' @noRd
.cells_in_region_for_view <- function(sl, gobject, space, region,
    relation = "intersects") {
    if (inherits(sl, "spatLocsObj")) {
        sl <- .apply_space_to_subobj(sl, gobject, space)
        return(.cells_in_region_dt(sl@coordinates, region, relation))
    }
    # Multi: named list of spatLocsObj, one per sample.
    if (is.list(sl)) {
        ids <- character(0L)
        for (samp_name in names(sl)) {
            child_sl <- sl[[samp_name]]
            if (!inherits(child_sl, "spatLocsObj")) next
            child_space <- .scope_space_to_sample_local(space, samp_name)
            child_sl <- .apply_space_to_subobj(child_sl, gobject,
                child_space)
            ids <- c(ids, .cells_in_region_dt(child_sl@coordinates,
                region, relation))
        }
        return(unique(ids))
    }
    stop("[.cells_in_region_for_view] unexpected sl type: ",
        toString(class(sl)), call. = FALSE)
}


# Cells in a region by centroid (numeric/SpatExtent fast path,
# SpatVector via AABB pre-filter + terra::is.related). Mirrors
# GiottoClass's unexported `.cells_in_region`. Used by the cache
# helper to fold viewCrop steps into the surviving cell_ID set.
#
#' @keywords internal
#' @noRd
.cells_in_region_dt <- function(sl_dt, region, relation = "intersects") {
    if (is.null(region)) return(sl_dt$cell_ID)
    if (!inherits(region, "SpatVector")) {
        ext <- if (inherits(region, "SpatExtent")) region[]
            else as.numeric(region)
        in_ext <- sl_dt$sdimx >= ext[[1L]] & sl_dt$sdimx <= ext[[2L]] &
                  sl_dt$sdimy >= ext[[3L]] & sl_dt$sdimy <= ext[[4L]]
        return(sl_dt$cell_ID[in_ext])
    }
    bbox <- terra::ext(region)[]
    in_bbox <- sl_dt$sdimx >= bbox[[1L]] & sl_dt$sdimx <= bbox[[2L]] &
               sl_dt$sdimy >= bbox[[3L]] & sl_dt$sdimy <= bbox[[4L]]
    candidates <- sl_dt[in_bbox, ]
    if (nrow(candidates) == 0L) return(character())
    pts <- terra::vect(
        as.matrix(candidates[, c("sdimx", "sdimy"), with = FALSE]),
        type = "points")
    surv <- terra::is.related(pts, region, relation)
    candidates$cell_ID[surv]
}


# Compute the intersected surviving cell_ID arrow Table for a view's
# filter + crop steps. Caches under a single slot `"surviving_cell_ids"`
# in `.cache` — mirrors `dataTableCoordinator`'s
# `.cached_surviving_cell_ids` convention (one final answer per
# resolution scope, regardless of how many steps contributed).
#
# Filter steps: each predicate's owner is located via
# `.find_store_with_cols`, narrowed, and the surviving id column
# intersected with the running set.
#
# viewCrop steps: cell_IDs whose centroid (from gobject's spatial_locs,
# projected through `space` if provided) satisfies (region, relation)
# get intersected. This is the in-mem coordinator's strategy and is
# what makes one `id_filter` op suffice for both filter and crop in a
# cache-mode resolve. For giottoPoints (non-cell-keyed) the cache path
# isn't applicable — those still go through spatRelate on the points'
# own geom.
#
# Returns NULL when the view has no contributing steps. Errors if any
# filter resolves against a non-`cell_ID` key (multi-key intersection
# out of scope; matches the in-mem coordinator's assumption).
#
# Trade-off: this path always materializes the owner side. For
# atlas-scale workflows where the owner is a parquetBase store and
# memory is tight, callers should pass `.cache = NULL` so the
# per-step paths in `.narrow_*_by_predicate` keep their lazy `[`-join
# branch.
#
#' @keywords internal
#' @noRd
.surviving_cell_ids_arrow <- function(view, gobject, .cache,
    spat_unit = NULL, feat_type = NULL) {
    ck <- "surviving_cell_ids"
    if (!is.null(.cache) && exists(ck, envir = .cache, inherits = FALSE)) {
        return(get(ck, envir = .cache))
    }
    filter_steps <- Filter(function(s) inherits(s, "viewFilter"),
        view@steps)
    crop_steps   <- Filter(function(s) inherits(s, "viewCrop"),
        view@steps)
    if (length(filter_steps) == 0L && length(crop_steps) == 0L) {
        if (!is.null(.cache)) assign(ck, NULL, envir = .cache)
        return(NULL)
    }

    surviving <- NULL  # arrow Table when non-NULL

    for (step in filter_steps) {
        narrowed <- .compute_narrowed_ids(step@predicate, gobject,
            spat_unit = spat_unit, feat_type = feat_type)
        if (!identical(narrowed$key, "cell_ID")) {
            stop("[.surviving_cell_ids_arrow] only cell_ID-keyed filter ",
                "steps are supported in cache mode (got '",
                narrowed$key, "')", call. = FALSE)
        }
        surviving <- if (is.null(surviving)) {
            narrowed$ids_tab
        } else {
            arrow::arrow_table(dplyr::collect(
                dplyr::semi_join(surviving, narrowed$ids_tab,
                    by = "cell_ID")
            ))
        }
    }

    if (length(crop_steps) > 0L) {
        # Predicate frame is read from view@space (the frame the crop
        # region was drawn in). Independent of any output space the
        # caller may have requested -- output space is applied by
        # `.apply_space_to_subobj` separately, on the actual subobject.
        pred_space <- if (!is.na(view@space)) {
            GiottoClass:::.resolve_space(gobject, view@space)
        } else NULL
        sl <- tryCatch(GiottoClass::getSpatialLocations(gobject,
            output = "spatLocsObj", spat_unit = spat_unit),
            error = function(e) NULL)
        if (is.null(sl)) {
            warning("[.surviving_cell_ids_arrow] viewCrop steps skipped: ",
                "no spatial locations available", call. = FALSE)
        } else {
            for (step in crop_steps) {
                step_ids <- .cells_in_region_for_view(sl, gobject,
                    pred_space, step@region, step@relation)
                step_tab <- arrow::arrow_table(data.frame(
                    cell_ID = step_ids, stringsAsFactors = FALSE))
                surviving <- if (is.null(surviving)) {
                    step_tab
                } else {
                    arrow::arrow_table(dplyr::collect(
                        dplyr::semi_join(surviving, step_tab,
                            by = "cell_ID")
                    ))
                }
            }
        }
    }

    if (!is.null(.cache)) assign(ck, surviving, envir = .cache)
    surviving
}


# Push a view-filter predicate onto `target_store` as a lazy op, narrowing
# across subobject boundaries when needed. Three branches:
#
#   (a) same-store `subset()` — predicate cols all on target;
#   (b) cross-store lazy `[`-join for parquetBase owners (no I/O at
#       coordinator time) — the scale-out path for 100M+ cells;
#   (c) cross-store eager-narrow-then-id_filter for DT owners.
#
# The cache path doesn't pass through here — it's handled at the
# `.push_view_to_pstore` level, which intersects all filter steps once
# per view and queues a single `id_filter`.
#
#' @keywords internal
#' @noRd
.narrow_store_by_predicate <- function(target_store, predicate, gobject,
    target_key = NULL, spat_unit = NULL, feat_type = NULL) {
    cols <- all.vars(predicate)

    if (length(cols) > 0L && all(cols %in% colnames(target_store))) {
        return(subset(target_store, predicate, quote = FALSE))
    }

    owner <- .find_store_with_cols(gobject, cols,
        spat_unit = spat_unit, feat_type = feat_type)
    if (is.null(owner)) {
        stop("[.narrow_store_by_predicate] no subobject covers ",
            "predicate cols: ", paste(cols, collapse = ", "),
            call. = FALSE)
    }

    join_by <- if (is.null(target_key) || identical(target_key, owner$key)) {
        owner$key
    } else {
        stats::setNames(owner$key, target_key)
    }

    if (identical(owner$kind, "parquetBase")) {
        narrowed_owner <- subset(owner$source, predicate, quote = FALSE)
        return(target_store[narrowed_owner,
            on = join_by, nomatch = NULL])
    }

    mask <- eval(predicate, envir = owner$source, enclos = baseenv())
    narrowed_dt <- owner$source[mask, owner$key, with = FALSE]
    narrowed_dt <- unique(narrowed_dt)
    ids_tab <- arrow::arrow_table(narrowed_dt)
    target_store@ops <- c(target_store@ops, list(list(
        type = "id_filter",
        ids_tab = ids_tab,
        by = join_by
    )))
    target_store
}


# Walk a view's steps and push each onto a parquetStore-inheriting store
# via the existing lazy-op API. Filter steps route through
# .narrow_store_by_predicate so they handle cross-storage cases. Returns
# the store with lazy ops queued. Pure lazy — no I/O for parquetBase
# owners; eager subset only for in-mem data.table owners.
#
# viewCrop / viewSampleSelect / spaceTransform pushdown land in
# follow-ups.
#
#' @keywords internal
#' @noRd
# Derive the composite 3x3 post-multiply affine matrix of a giottoSpace
# for a given gobject by applying its transforms to three known basis
# points. Returns NULL when the space has no steps applicable to this
# gobject (no matching sample key, or empty step list).
#
# Post-multiply convention (matches GiottoClass affine2d): [x, y, 1] %*% M
# yields [x', y', 1], so:
#   M[1, 1] = coef on x in x';  M[1, 2] = coef on x in y'
#   M[2, 1] = coef on y in x';  M[2, 2] = coef on y in y'
#   M[3, 1] = tx;               M[3, 2] = ty
# Derivation:
#   (0,0) -> (tx, ty)
#   (1,0) -> (M[1,1] + tx, M[1,2] + ty)
#   (0,1) -> (M[2,1] + tx, M[2,2] + ty)
#' @keywords internal
#' @noRd
.space_composite_affine <- function(space, gobject) {
    if (is.null(space)) return(NULL)
    key <- GiottoClass:::.space_sample_key_for(gobject, space)
    if (is.null(key)) return(NULL)
    steps <- space@samples[[key]]
    if (length(steps) == 0L) return(NULL)
    test_pts <- terra::vect(
        matrix(c(0, 0, 1, 0, 0, 1), ncol = 2L, byrow = TRUE),
        type = "points")
    out_pts <- test_pts
    for (step in steps) {
        out_pts <- do.call(step@op, c(list(x = out_pts), step@args))
    }
    p <- terra::crds(out_pts)
    M <- matrix(0, nrow = 3L, ncol = 3L)
    M[1L, 1L] <- p[2L, 1L] - p[1L, 1L]
    M[1L, 2L] <- p[2L, 2L] - p[1L, 2L]
    M[2L, 1L] <- p[3L, 1L] - p[1L, 1L]
    M[2L, 2L] <- p[3L, 2L] - p[1L, 2L]
    M[3L, 1L] <- p[1L, 1L]
    M[3L, 2L] <- p[1L, 2L]
    M[3L, 3L] <- 1
    M
}

# Project a SpatVector region from `from_space`'s frame to `to_space`'s
# frame. Composition: y_native = inv(M_from) @ y_from;  y_to = M_to @ y_native.
# Either space may be NULL -- a NULL space means "the gobject's native
# frame," i.e. no transform on that side.
#' @keywords internal
#' @noRd
.project_region_between_spaces <- function(y, gobject,
    from_space = NULL, to_space = NULL) {
    if (is.null(from_space) && is.null(to_space)) return(y)
    if (!inherits(y, "SpatVector")) {
        if (is.numeric(y))         y <- terra::as.polygons(terra::ext(y))
        if (inherits(y, "SpatExtent")) y <- terra::as.polygons(y)
    }
    m_from <- .space_composite_affine(from_space, gobject)
    m_to   <- .space_composite_affine(to_space,   gobject)
    # If both affines collapse to the identical matrix the composition is
    # identity -- skip the round-trip to avoid unnecessary float drift.
    if (!is.null(m_from) && !is.null(m_to) &&
        isTRUE(all.equal(m_from, m_to))) {
        return(y)
    }
    if (!is.null(m_from)) {
        y <- GiottoClass::affine(y, m_from, inv = TRUE)
    }
    if (!is.null(m_to)) {
        y <- GiottoClass::affine(y, m_to)
    }
    y
}

.push_view_to_pstore <- function(store, view, gobject,
    target_key = NULL, space = NULL, spat_unit = NULL, feat_type = NULL,
    .cache = NULL) {
    if (is.null(view) || length(view@steps) == 0L) return(store)

    # The `space` arg is the OUTPUT frame: when non-NULL the caller has
    # already composed it into the store's @post_ops via
    # .apply_space_to_subobj. The PREDICATE frame is read from view@space
    # and used here to project any crop regions into the same frame as
    # the geom column that will be tested -- i.e. into the OUTPUT frame,
    # since that's what the geom column ends up in after @post_ops apply.
    predicate_space <- if (!is.na(view@space)) {
        GiottoClass:::.resolve_space(gobject, view@space)
    } else NULL

    # Cache path: a single arrow Table holding cell_IDs that survive
    # ALL filter + crop steps in the view (via the in-mem coordinator's
    # spatial_locs-centroid pattern for crops). Queues ONE id_filter
    # on the target and short-circuits — the per-step loop below is
    # only used when no cache is provided. Cache mode requires a
    # cell-keyed target (target_key cell_ID or poly_ID conventionally
    # equal); callers whose target isn't cell-keyed (giottoPoints, the
    # current non-aggregated case) explicitly pass `.cache = NULL` to
    # route viewCrop through spatRelate on the points' geom instead.
    if (!is.null(.cache)) {
        surv <- .surviving_cell_ids_arrow(view, gobject, .cache,
            spat_unit = spat_unit, feat_type = feat_type)
        if (!is.null(surv)) {
            join_by <- if (is.null(target_key) ||
                identical(target_key, "cell_ID")) {
                "cell_ID"
            } else {
                stats::setNames("cell_ID", target_key)
            }
            store@ops <- c(store@ops, list(list(
                type = "id_filter",
                ids_tab = surv,
                by = join_by
            )))
        }
        return(store)
    }

    for (step in view@steps) {
        if (inherits(step, "viewFilter")) {
            store <- .narrow_store_by_predicate(
                target_store = store,
                predicate = step@predicate,
                gobject = gobject,
                target_key = target_key,
                spat_unit = spat_unit,
                feat_type = feat_type
            )
        } else if (inherits(step, "viewCrop")) {
            if (!inherits(store, "parquetGeomBase")) {
                warning("[push_view_to_pstore] viewCrop step skipped: ",
                    "target store is not parquetGeomBase (no geom ",
                    "column to evaluate the predicate on)", call. = FALSE)
                next
            }
            region <- step@region
            # spatRelate accepts SpatVector / sf / giottoPolygon /
            # spatLocsObj / WKT / parquetGeomBase. Promote numeric /
            # SpatExtent to a polygon SpatVector first.
            y <- region
            if (is.numeric(region)) {
                y <- terra::as.polygons(terra::ext(region))
            } else if (inherits(region, "SpatExtent")) {
                y <- terra::as.polygons(region)
            }
            # Project region from predicate frame → output frame so the
            # spat_relate eval (which runs against geom in the output
            # frame, via @post_ops) sees both sides in the same frame.
            y <- .project_region_between_spaces(y, gobject,
                from_space = predicate_space, to_space = space)
            store <- spatRelate(store, y, relation = step@relation)
        }
        # viewSampleSelect / space transform pushdown follow
    }
    store
}


# Arrow-bridge narrowing for in-mem data.table targets. Wraps `target_dt`
# in an arrow Table (zero-copy for shared column buffers), applies the
# predicate via arrow's dplyr backend, and collects back to data.table.
#
# Three branches mirror .narrow_store_by_predicate:
#   (a) same-store: filter target_arrow directly, collect, setDT.
#   (b) parquet owner: lazy-narrow the owner via subset(); materialize
#       ONLY the surviving join-key column via storeRead(output = "tibble");
#       semi_join into target_arrow; collect.
#   (c) data.table owner: eager-narrow the owner data.table; project to
#       join key; semi_join into target_arrow; collect.
#
# Predicate cols spanning multiple owners are not supported in v1 — first
# match wins via .find_store_with_cols.
#
#' @keywords internal
#' @noRd
.narrow_dt_via_arrow <- function(target_dt, predicate, gobject,
    key = "cell_ID", spat_unit = NULL, feat_type = NULL) {
    cols <- all.vars(predicate)
    target_arrow <- arrow::arrow_table(target_dt)

    if (length(cols) > 0L && all(cols %in% names(target_dt))) {
        out <- dplyr::collect(
            dplyr::filter(target_arrow, !!predicate)
        )
        return(data.table::setDT(out))
    }

    owner <- .find_store_with_cols(gobject, cols,
        spat_unit = spat_unit, feat_type = feat_type)
    if (is.null(owner)) {
        stop("[.narrow_dt_via_arrow] no subobject covers predicate cols: ",
            paste(cols, collapse = ", "), call. = FALSE)
    }
    if (!identical(owner$key, key)) {
        stop("[.narrow_dt_via_arrow] target key '", key,
            "' does not match owner key '", owner$key,
            "' — cross-key narrowing is out of v1 scope", call. = FALSE)
    }

    if (identical(owner$kind, "parquetBase")) {
        narrowed_owner <- subset(owner$source, predicate, quote = FALSE)
        ids_tbl <- storeRead(narrowed_owner[, owner$key, drop = FALSE],
            output = "tibble")
        ids_arrow <- arrow::arrow_table(ids_tbl)
        out <- dplyr::collect(
            dplyr::semi_join(target_arrow, ids_arrow, by = owner$key)
        )
        return(data.table::setDT(out))
    }

    mask <- eval(predicate, envir = owner$source, enclos = baseenv())
    narrowed_dt <- unique(owner$source[mask, owner$key, with = FALSE])
    ids_arrow <- arrow::arrow_table(narrowed_dt)
    out <- dplyr::collect(
        dplyr::semi_join(target_arrow, ids_arrow, by = owner$key)
    )
    data.table::setDT(out)
}


# Apply a giottoSpace's per-sample transform steps to a subobject.
#
# Dispatch routing:
#
# * Backed `giottoPolygon` / `giottoPoints` (whose @spatVector inherits
#   parquetBase): dispatch transforms DIRECTLY on the inner store. The
#   wrapper methods on giottoPolygon/giottoPoints (`.shift_gpoly`,
#   `.do_gpoly + terra::spin`, `.affine_sv`, etc.) are SpatVector-
#   specific and call terra:: functions that have no parquetGeomBase
#   methods. The parquetGeomBase has its own transform methods (`spin`,
#   `affine`, `spatShift`, `rescale`, `shear`, `flip`, `t`) that
#   compose into @post_ops via matrix product — multi-step recipes
#   accumulate into a single affine2d at storeRead time.
#
# * In-mem / DT-backed subobjects (spatLocsObj, in-mem giottoPolygon):
#   dispatch on the subobj itself; existing GiottoClass methods mutate
#   coordinates / SpatVector eagerly.
#
# Affine matrix coercion: a `spaceTransform` whose `@op == "affine"`
# records `args = list(y = <matrix>)`. parquetGeomBase has an
# `(parquetGeomBase, affine2d)` method but NOT `(parquetGeomBase,
# matrix)` (would dispatch fail otherwise). Wrap the matrix in an
# affine2d before dispatch.
#
# Sample-key resolution: use the `:default:` sentinel for single-
# giotto contexts; the single key if just one is present; otherwise
# a no-op.
#
# TODO: When GiottoClass exposes `.apply_space_to_subobj` (or adds an
# `applySpace` generic dispatching on coordinator), drop this duplicate
# and replace the parquetGeomBase routing with a coordinator-specific
# method override.
#
#' @keywords internal
#' @noRd
# Narrow a giottoSpace's @samples list to a single sample key, renaming
# it to the `:default:` sentinel. Mirror of GiottoClass's unexported
# `.scope_space_to_sample`. Used by the multi-sample crop path in
# `.surviving_cell_ids_arrow` so each child's spatLocsObj is projected
# through that sample's transforms before centroid-in-region check.
#
#' @keywords internal
#' @noRd
.scope_space_to_sample_local <- function(space, sample_name) {
    if (is.null(space)) return(NULL)
    keys <- names(space@samples)
    pick <- if (sample_name %in% keys) sample_name
        else if (":default:" %in% keys) ":default:"
        else NULL
    if (is.null(pick)) return(NULL)
    out <- space
    out@samples <- stats::setNames(list(space@samples[[pick]]), ":default:")
    out
}


.apply_space_to_subobj <- function(subobj, gobject, space) {
    if (is.null(space)) return(subobj)
    keys <- names(space@samples)
    if (length(keys) == 0L) return(subobj)
    key <- if (":default:" %in% keys) {
        ":default:"
    } else if (length(keys) == 1L) {
        keys[[1L]]
    } else {
        return(subobj)
    }
    steps <- space@samples[[key]]
    backed_geom <- .hasSlot(subobj, "spatVector") &&
        inherits(subobj@spatVector, "parquetBase")
    for (step in steps) {
        args <- step@args
        if (backed_geom) {
            if (identical(step@op, "affine") &&
                inherits(args$y, "matrix")) {
                # ANY,missing affine method wraps a matrix into an affine2d
                args$y <- affine(args$y)
            }
            subobj@spatVector <- do.call(step@op,
                c(list(x = subobj@spatVector), args))
        } else {
            subobj <- do.call(step@op, c(list(x = subobj), args))
        }
    }
    subobj
}


# Walk a view's filter + crop steps and apply each to a data.table target
# via the arrow bridge (filters) or the cells-in-region helper (crops).
# `viewSampleSelect` / `spaceTransform` follow the same pattern in future
# revisions.
#
#' @keywords internal
#' @noRd
.push_view_to_dt <- function(dt, view, gobject, key = "cell_ID",
    spat_unit = NULL, feat_type = NULL, .cache = NULL) {
    if (is.null(view) || length(view@steps) == 0L) return(dt)
    if (is.null(dt) || nrow(dt) == 0L) return(dt)

    # Cache path: one semi_join against the view-wide intersection.
    if (!is.null(.cache)) {
        surv <- .surviving_cell_ids_arrow(view, gobject, .cache,
            spat_unit = spat_unit, feat_type = feat_type)
        if (is.null(surv)) return(dt)
        if (!identical(key, "cell_ID")) {
            stop("[push_view_to_dt] cache mode supports cell_ID-keyed ",
                "targets only (got key = '", key, "')", call. = FALSE)
        }
        out <- dplyr::collect(
            dplyr::semi_join(arrow::arrow_table(dt), surv, by = "cell_ID")
        )
        return(data.table::setDT(out))
    }

    # No-cache path: per-step narrowing. Filter steps go through the
    # arrow-bridge predicate helper. Crop steps narrow by cell_IDs whose
    # centroid satisfies the predicate against the region -- the same
    # semantics as the in-memory dataTableCoordinator path (centroid-in-
    # region), and the cache path (.surviving_cell_ids_arrow). The
    # predicate frame is `view@space`; the output frame (any `space=`
    # arg on the getter) is applied separately by the resolveSubobject
    # method via `.apply_space_to_subobj`.
    crop_steps <- Filter(function(s) inherits(s, "viewCrop"), view@steps)
    have_crops <- length(crop_steps) > 0L
    if (have_crops && !identical(key, "cell_ID")) {
        warning("[push_view_to_dt] viewCrop steps skipped on non-cell-",
            "keyed target (key = '", key, "')", call. = FALSE)
        have_crops <- FALSE
    }
    pred_space <- if (have_crops && !is.na(view@space)) {
        GiottoClass:::.resolve_space(gobject, view@space)
    } else NULL
    sl <- if (have_crops) {
        tryCatch(GiottoClass::getSpatialLocations(gobject,
            output = "spatLocsObj", spat_unit = spat_unit),
            error = function(e) NULL)
    } else NULL
    if (have_crops && is.null(sl)) {
        warning("[push_view_to_dt] viewCrop steps skipped: no spatial ",
            "locations available", call. = FALSE)
        have_crops <- FALSE
    }

    for (step in view@steps) {
        if (inherits(step, "viewFilter")) {
            dt <- .narrow_dt_via_arrow(dt, step@predicate, gobject,
                key = key, spat_unit = spat_unit, feat_type = feat_type)
        } else if (inherits(step, "viewCrop") && have_crops) {
            step_ids <- .cells_in_region_for_view(sl, gobject, pred_space,
                step@region, step@relation)
            dt <- dt[cell_ID %in% step_ids]
        }
    }
    dt
}


# resolveSubobject methods (parquetCoordinator) ####
#
# DT-target subobjects (cellMetaObj / featMetaObj / spatLocsObj /
# spatEnrObj) currently hold their data in strict-typed data.table slots
# that cannot hold a parquetStore directly. Until the slot types widen,
# parquetCoordinator handles them via the arrow bridge
# (.push_view_to_dt / .narrow_dt_via_arrow) rather than falling through
# to dataTableCoordinator — preserves the cross-storage path (DT target,
# parquet owner) end-to-end. Same-store predicates collapse to a single
# arrow filter + collect, which is comparable to data.table on 1M-row
# tables and consistent with the rest of the parquetCoordinator pipeline.

#' @rdname parquetCoordinator-class
#' @importFrom GiottoClass resolveSubobject
#' @export
setMethod("resolveSubobject",
    signature(subobj = "cellMetaObj",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        if (is.null(view) || length(view@steps) == 0L) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        if (inherits(subobj@metaDT, "dataStore")) {
            subobj@metaDT <- .push_view_to_pstore(subobj@metaDT, view,
                gobject = gobject, space = space,
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        } else {
            subobj@metaDT <- .push_view_to_dt(subobj@metaDT, view, gobject,
                key = "cell_ID",
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        }
        subobj
    }
)

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveSubobject",
    signature(subobj = "featMetaObj",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        if (is.null(view) || length(view@steps) == 0L) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        if (inherits(subobj@metaDT, "dataStore")) {
            subobj@metaDT <- .push_view_to_pstore(subobj@metaDT, view,
                gobject = gobject, space = space,
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        } else {
            subobj@metaDT <- .push_view_to_dt(subobj@metaDT, view, gobject,
                key = "feat_ID",
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        }
        subobj
    }
)

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveSubobject",
    signature(subobj = "spatLocsObj",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        # OUTPUT-space transform: only applied when the caller explicitly
        # passed `space=`. `view@space` is NOT consulted here -- it only
        # affects the predicate frame inside the crop-step handlers.
        if (!is.null(space)) {
            subobj <- .apply_space_to_subobj(subobj, gobject, space)
        }
        if (is.null(view) || length(view@steps) == 0L) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        if (inherits(subobj@coordinates, "dataStore")) {
            subobj@coordinates <- .push_view_to_pstore(
                subobj@coordinates, view, gobject = gobject,
                space = space, spat_unit = spat_unit, .cache = .cache)
        } else {
            subobj@coordinates <- .push_view_to_dt(
                subobj@coordinates, view, gobject,
                key = "cell_ID", spat_unit = spat_unit,
                .cache = .cache)
        }
        subobj
    }
)

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveSubobject",
    signature(subobj = "spatEnrObj",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        if (is.null(view) || length(view@steps) == 0L) return(subobj)
        if (is.null(subobj@enrichDT)) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        if (inherits(subobj@enrichDT, "dataStore")) {
            subobj@enrichDT <- .push_view_to_pstore(subobj@enrichDT, view,
                gobject = gobject, space = space,
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        } else {
            subobj@enrichDT <- .push_view_to_dt(subobj@enrichDT, view,
                gobject, key = "cell_ID",
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        }
        subobj
    }
)


# giottoPolygon / giottoPoints have @spatVector = "ANY" so they CAN hold a
# parquetGeomStore today. When backed, push view filters via the existing
# lazy-op API. When in-mem (terra SpatVector / sf), fall through to
# dataTableCoordinator -- which only applies viewCrop (geometry) on those
# subobject classes.
#
# Join-key convention: a polygon's `poly_ID` aligns with the cells
# spat_unit's `cell_ID`, so cross-store narrowing against a cellMeta-keyed
# owner uses `target_key = "poly_ID"` to map (target.poly_ID =
# owner.cell_ID) in the queued join op.

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveSubobject",
    signature(subobj = "giottoPolygon",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        if (!inherits(subobj@spatVector, "parquetBase")) {
            # In-mem SpatVector path: defer to dataTableCoordinator
            # (handles space + crop geometrically).
            return(callNextMethod())
        }
        # OUTPUT-space transform: only applied when caller passed `space=`.
        # `view@space` is consumed inside `.push_view_to_pstore` as the
        # predicate frame -- the crop region is interpreted there before
        # being projected into the output frame for the spat_relate eval.
        if (!is.null(space)) {
            subobj <- .apply_space_to_subobj(subobj, gobject, space)
        }
        if (is.null(view) || length(view@steps) == 0L) return(subobj)
        # TODO: decouple semantic routing (centroid vs polygon-geom) from
        # cache memoization in `.push_view_to_pstore`. The right routing
        # decision is per-step, based on the predicate relation (centroid
        # is approximately correct only for intersects/disjoint; within /
        # contains / covers / overlaps / touches / crosses all require
        # the geom). Target storage kind (DT vs parquet) is NOT the right
        # discriminator -- a DT-target subobject can still narrow by a
        # geom-mode predicate, because the geom evaluation runs on the
        # gobject's polygon source and the resulting cell_ID set narrows
        # the DT downstream. Cache is purely an optimization layer under
        # the centroid path; it should not gate semantics. Once routed
        # per-step, the two patches below (`force cache for polygons`
        # here, and `force NULL cache for points` in the points
        # resolveSubobject) both go away.
        #
        # Force the centroid-based cell_ID narrow path even when the
        # caller didn't pass a cache. The polygon resolveSubobject is
        # cell-aggregatable (poly_ID == cell_ID by convention) so the
        # surviving cell_ID set computed from spatLocs centroids applies
        # directly. Allocating a one-shot cache here routes
        # `.push_view_to_pstore` through `.surviving_cell_ids_arrow` --
        # the same path the batched `materialize(g, view)` takes -- so
        # the polygon narrow always agrees with the spatLocs / cellMeta
        # / expression narrow for the same view recipe. (The per-step
        # spat_relate-on-geom path is the right default for giottoPoints
        # below, where poly_ID == cell_ID does NOT hold.)
        .cache <- list(...)$.cache %null% new.env(parent = emptyenv())
        subobj@spatVector <- .push_view_to_pstore(subobj@spatVector, view,
            gobject = gobject, target_key = "poly_ID", space = space,
            spat_unit = GiottoClass::spatUnit(subobj),
            .cache = .cache)
        subobj
    }
)

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveSubobject",
    signature(subobj = "giottoPoints",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        if (!inherits(subobj@spatVector, "parquetBase")) {
            return(callNextMethod())
        }
        # OUTPUT-space transform: only applied when caller passed `space=`.
        if (!is.null(space)) {
            subobj <- .apply_space_to_subobj(subobj, gobject, space)
        }
        if (is.null(view) || length(view@steps) == 0L) return(subobj)
        .cache <- list(...)$.cache
        # viewCrop reaches parity with in-mem via parquetGeomBase's crop().
        #
        # TODO: viewFilter on giottoPoints is unreachable in real
        # workflows today — `viewFilter` is cell-centric by convention
        # (the in-mem method skips it). Same-store predicates like
        # `feature_name == "GENE1"` work mechanically here but no view
        # consumer drives feature subsetting through viewFilter — would
        # belong in a future `viewFeatSelect` step. Cross-store
        # narrowing via cellMeta requires an aggregated points store
        # carrying cell_ID, which the current pipeline doesn't produce.
        # Code is in place for when either of these gaps closes.
        # Bypass the cache path: giottoPoints' geom is point-level, not
        # cell-aggregated, so the cell_ID surviving set computed against
        # spatial_locs doesn't directly apply. viewCrop on points goes
        # through spatRelate on the points' own geom; viewFilter is
        # unreachable today (see TODO above).
        subobj@spatVector <- .push_view_to_pstore(subobj@spatVector, view,
            gobject = gobject, target_key = "cell_ID",
            feat_type = GiottoClass::featType(subobj),
            .cache = NULL)
        subobj
    }
)


# Matrix-shaped subobjects (exprObj, dimObj) ####
#
# These narrow by cell_ID via the matrix indexing path rather than the
# @ops queue. parquetExprStore has its own `[, j]` operator that updates
# @cell_idx lazily (no Parquet rewrite, no I/O). dimObj@coordinates is
# typically an in-mem matrix with cell_IDs as rownames — `coords[ids, ]`
# is the in-mem narrow.
#
# Both methods use the same cached surviving cell_ID set computed by
# `.surviving_cell_ids_arrow` (filter + crop folded together). The
# materialization of the cell_ID vector is mandatory here — matrix
# indexing can't accept an arrow Table, so we collect once. The cache
# means this collect happens once per `materialize()` call regardless
# of how many matrix-shaped subobjects share it.

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveSubobject",
    signature(subobj = "exprObj",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        if (is.null(view) || length(view@steps) == 0L) return(subobj)
        if (!inherits(subobj@exprMat, "parquetExprStore")) {
            # In-mem matrix / dgCMatrix / etc. — delegate to
            # dataTableCoordinator's exprObj method.
            return(callNextMethod())
        }
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        surv <- .surviving_cell_ids_arrow(view, gobject, .cache,
            spat_unit = spat_unit, feat_type = feat_type)
        if (is.null(surv)) return(subobj)

        surv_ids <- dplyr::collect(surv)$cell_ID
        keep <- intersect(surv_ids, subobj@exprMat@cell_ids)
        if (length(keep) == 0L) {
            # Empty narrow — give back a fresh zero-cell view of the store
            subobj@exprMat <- subobj@exprMat[, integer(0L)]
        } else {
            subobj@exprMat <- subobj@exprMat[, keep]
        }
        subobj
    }
)

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveSubobject",
    signature(subobj = "dimObj",
              coordinator = "parquetCoordinator"),
    function(subobj, gobject, view, space, coordinator, ...) {
        # dimObj@coordinates is `ANY` but typically a matrix with
        # cell_IDs as rownames; backed (parquetBase) coordinates aren't
        # currently produced by the pipeline. Delegate to
        # dataTableCoordinator for in-mem matrices.
        if (inherits(subobj@coordinates, "dataStore")) {
            # Future: queue id_filter via @ops if backed. For now,
            # fall through — no consumer produces this shape today.
            return(callNextMethod())
        }
        return(callNextMethod())
    }
)
