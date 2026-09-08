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


# Centroid table for the crop steps, in the PREDICATE frame.
#
# This is a fetch helper only. The predicate itself is GiottoClass's
# `.cells_in_crop_step()`, which owns the crop semantics for both
# `geom` arms — including the `disjoint` AABB exclusion and the
# rectangle fast path. GiottoDisk deliberately does NOT keep its own
# copy of that: an earlier duplicate here silently diverged (it had
# neither the disjoint fix nor a geom arm), so a backed object answered
# a different question than an in-memory one for the same recipe.
#
# What is local is the FETCH, because `spatLocsObj@coordinates` may hold
# a store, which GiottoClass's `.get_projected_spatlocs()` cannot read.
#
# * single giotto: `getSpatialLocations(g, output = "spatLocsObj")`
#   returns ONE spatLocsObj. Project through `space` once, then hand
#   over the coordinates.
#
# * giottoMulti: the same getter returns a NAMED LIST of per-sample
#   spatLocsObjs. Per sample: scope the space to that sample, project,
#   and prefix cell_IDs with `<sample>::` so the result speaks the joint
#   `@cell_metadata` vocabulary. Rows are stacked, not unioned by ID —
#   joint cell_IDs are globally unique, so the stack is the joint table.
#
#' @keywords internal
#' @noRd
.projected_spatlocs_dt <- function(gobject, space, spat_unit = NULL) {
    cell_ID <- NULL # NSE
    sl <- tryCatch(GiottoClass::getSpatialLocations(gobject,
        output = "spatLocsObj", spat_unit = spat_unit),
        error = function(e) NULL)
    if (is.null(sl)) return(NULL)

    .coords <- function(x) {
        co <- x@coordinates
        if (inherits(co, "dataStore")) {
            co <- data.table::setDT(storeRead(co, output = "tibble"))
        }
        data.table::as.data.table(co)
    }

    if (inherits(sl, "spatLocsObj")) {
        return(.coords(.apply_space_to_subobj(sl, gobject, space)))
    }
    if (!is.list(sl)) {
        stop("[.projected_spatlocs_dt] unexpected spatial locations type: ",
            toString(class(sl)), call. = FALSE)
    }
    parts <- lapply(names(sl), function(nm) {
        child <- sl[[nm]]
        if (!inherits(child, "spatLocsObj")) return(NULL)
        child <- .apply_space_to_subobj(child, gobject,
            GiottoClass:::.scope_space_to_sample(space, nm))
        dt <- .coords(child)
        dt[, cell_ID := paste(nm, cell_ID, sep = "::")]
        dt
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) return(NULL)
    data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
}


# Surviving-cell_ID arrow Tables, memoized per resolution scope.
#
# THREE cached quantities, not one, and every one of them is
# TARGET-INDEPENDENT -- which is the property that makes a single
# `.cache` safe to share across every subobject in one `materialize()`:
#
#   "filter_ids"         cell_IDs surviving all filter steps
#   "crop_ids:centroid"  cell_IDs surviving the crop steps that declared
#                        `geom = "centroid"`
#   "crop_ids:poly"      cell_IDs surviving the crop steps that declared
#                        `geom = "poly"`
#
# The split exists because different targets consume different subsets. A
# backed `giottoPolygon` store evaluates its own poly-arm crops lazily on
# its `geom` column -- that is the whole pushdown win -- while still
# needing the filter and centroid arms as an eager id set. Folding
# everything into one number would force the polygon store to materialize
# a predicate it could have pushed down.
#
# Filter steps: each predicate's owner is located via
# `.find_store_with_cols`, narrowed, and the surviving id column
# intersected with the running set. This is the arrow / cross-storage
# path, which is why it is not delegated to GiottoClass.
#
# Crop steps: delegated to `GiottoClass:::.cells_in_crop_step()`, which
# owns crop semantics for both `geom` arms. Only the centroid table is
# fetched locally, because a backed `@coordinates` needs `storeRead`.
#
# Each returns NULL when the view has no contributing steps of that kind,
# meaning "unconstrained" rather than "empty". Errors if a filter
# resolves against a non-`cell_ID` key -- multi-key intersection is out
# of scope, matching the in-mem coordinator's assumption.
#
# Trade-off unchanged from the original: the filter path materializes the
# owner side. For atlas-scale workflows where the owner is a parquetBase
# store and memory is tight, callers pass `.cache = NULL` so the per-step
# paths in `.narrow_*_by_predicate` keep their lazy `[`-join branch. That
# is the ONLY thing `.cache` decides -- memoization, plus an eager/lazy
# choice for filters. It is never a semantic switch; crop semantics are
# read off the step's declared `geom`.

#' @keywords internal
#' @noRd
.memo <- function(.cache, key, compute) {
    if (!is.null(.cache) && exists(key, envir = .cache, inherits = FALSE)) {
        return(get(key, envir = .cache))
    }
    val <- compute()
    if (!is.null(.cache)) assign(key, val, envir = .cache)
    val
}

# Intersect two id sets, either of which may be NULL (= unconstrained).
#' @keywords internal
#' @noRd
.intersect_ids_arrow <- function(a, b) {
    if (is.null(a)) return(b)
    if (is.null(b)) return(a)
    arrow::arrow_table(dplyr::collect(
        dplyr::semi_join(a, b, by = "cell_ID")))
}

#' @keywords internal
#' @noRd
.ids_arrow <- function(ids) {
    arrow::arrow_table(data.frame(cell_ID = ids, stringsAsFactors = FALSE))
}

#' @keywords internal
#' @noRd
.view_filter_ids_arrow <- function(view, gobject, .cache,
    spat_unit = NULL, feat_type = NULL) {
    .memo(.cache, "filter_ids", function() {
        steps <- GiottoClass:::.view_steps_of(view, "filter")
        if (length(steps) == 0L) return(NULL)
        surviving <- NULL
        for (step in steps) {
            # Q7 records the predicate deparsed, so it comes back as a
            # string and has to be re-parsed before `all.vars()` /
            # `eval()` downstream.
            narrowed <- .compute_narrowed_ids(str2lang(step$predicate),
                gobject, spat_unit = spat_unit, feat_type = feat_type)
            if (!identical(narrowed$key, "cell_ID")) {
                stop("[.view_filter_ids_arrow] only cell_ID-keyed filter ",
                    "steps are supported in cache mode (got '",
                    narrowed$key, "')", call. = FALSE)
            }
            surviving <- .intersect_ids_arrow(surviving, narrowed$ids_tab)
        }
        surviving
    })
}

#' @keywords internal
#' @noRd
.view_crop_ids_arrow <- function(view, gobject, .cache, coordinator,
    geom = c("centroid", "poly"), spat_unit = NULL) {
    geom <- match.arg(geom)
    .memo(.cache, paste0("crop_ids:", geom), function() {
        steps <- Filter(function(s) identical(s$geom, geom),
            GiottoClass:::.view_steps_of(view, "crop"))
        if (length(steps) == 0L) return(NULL)

        # Predicate frame is read from `view$space` -- the frame the crop
        # region was drawn in. Independent of any output space the caller
        # requested; that one is applied by `.apply_space_to_subobj` on
        # the actual subobject.
        pred_space <- if (!is.na(view$space)) {
            GiottoClass:::.resolve_space(gobject, view$space)
        } else NULL

        sl_dt <- if (identical(geom, "centroid")) {
            .projected_spatlocs_dt(gobject, pred_space,
                spat_unit = spat_unit)
        } else NULL
        if (identical(geom, "centroid") && is.null(sl_dt)) {
            warning("[.view_crop_ids_arrow] crop steps skipped: no ",
                "spatial locations available", call. = FALSE)
            return(NULL)
        }

        surviving <- NULL
        for (step in steps) {
            ids <- GiottoClass:::.cells_in_crop_step(gobject, step, sl_dt,
                pred_space, coordinator, spat_unit = spat_unit)
            surviving <- .intersect_ids_arrow(surviving, .ids_arrow(ids))
        }
        surviving
    })
}

# The full eager answer: filters plus BOTH crop arms. What a target with
# no geometry of its own (cellMeta, spatLocs, spatEnr, expression) needs,
# and it collapses to a single queued `id_filter`.
#' @keywords internal
#' @noRd
.surviving_cell_ids_arrow <- function(view, gobject, .cache, coordinator,
    spat_unit = NULL, feat_type = NULL) {
    .memo(.cache, "surviving_cell_ids", function() {
        ids <- .view_filter_ids_arrow(view, gobject, .cache,
            spat_unit = spat_unit, feat_type = feat_type)
        for (g in c("centroid", "poly")) {
            ids <- .intersect_ids_arrow(ids,
                .view_crop_ids_arrow(view, gobject, .cache, coordinator,
                    geom = g, spat_unit = spat_unit))
        }
        ids
    })
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
    steps <- space$samples[[key]]
    if (length(steps) == 0L) return(NULL)
    # Probe carrier is a spatLocsObj, not a bare SpatVector: GiottoClass
    # implements all seven transform generics on spatLocsObj but not on
    # SpatVector (`spatShift` has no SpatVector method), so a SpatVector
    # probe dies on the most common step there is. spatLocsObj is also
    # the carrier the centroid path already uses, so the derived matrix
    # is measured through the same methods the real data goes through
    # rather than a parallel implementation.
    probe <- GiottoClass::createSpatLocsObj(
        data.table::data.table(
            cell_ID = c("o", "x", "y"),
            sdimx = c(0, 1, 0), sdimy = c(0, 0, 1)),
        name = "probe", verbose = FALSE)
    for (step in steps) {
        probe <- do.call(step$op, c(list(x = probe), step$args))
    }
    p <- as.matrix(probe@coordinates[, c("sdimx", "sdimy")])
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

# Push a view's steps onto a parquetStore-inheriting store as lazy ops.
#
# ROUTING (A7). Two independent decisions, and neither is made by
# `.cache`:
#
#   * a CROP step is routed by its declared `geom`. `geom = "poly"` on a
#     store that carries the cell polygon is evaluated lazily by the
#     store itself, via `spatRelate()` on its `geom` column -- no
#     materialization, engine chosen later at `storeRead()`. Everything
#     else reduces to an eager cell_ID set (`.view_crop_ids_arrow`) and
#     is queued as an `id_filter`.
#
#   * a FILTER step is routed by `.cache`, which is purely an
#     eager/lazy performance choice: with a cache all filters fold into
#     one memoized `id_filter`; without one each narrows the store on
#     its own, which keeps `.narrow_store_by_predicate`'s lazy
#     cross-store `[`-join available for atlas-scale owners.
#
# `cell_keyed` says whether one store row IS one cell -- true for
# cellMeta / spatLocs / a cell-polygon store, false for a transcript
# points store. It is the reason points behave differently, and it is
# NOT a proxy for storage kind. Before A7 this distinction was smuggled
# through `.cache` by two per-class patches at the call sites (force a
# cache for polygons, force `NULL` for points), which made the cache
# decide semantics: a backed `giottoPolygon` could not honour
# `geom = "poly"` at all, because the forced cache always answered the
# centroid question instead. Both patches are gone.
#
#' @keywords internal
#' @noRd
.push_view_to_pstore <- function(store, view, gobject, coordinator,
    target_key = NULL, space = NULL, spat_unit = NULL, feat_type = NULL,
    cell_keyed = TRUE, .cache = NULL) {
    if (is.null(view) || length(view$steps) == 0L) return(store)

    join_by <- if (is.null(target_key) ||
        identical(target_key, "cell_ID")) {
        "cell_ID"
    } else {
        stats::setNames("cell_ID", target_key)
    }
    .queue_ids <- function(store, ids_tab) {
        if (is.null(ids_tab)) return(store)
        store@ops <- c(store@ops, list(list(
            type = "id_filter", ids_tab = ids_tab, by = join_by)))
        store
    }

    # The `space` arg is the OUTPUT frame: when non-NULL the caller has
    # already composed it into the store's @post_ops via
    # .apply_space_to_subobj. The PREDICATE frame is read from
    # `view$space` and used here to project any crop regions into the
    # same frame as the geom column that will be tested -- i.e. into the
    # OUTPUT frame, since that is what the geom column ends up in once
    # @post_ops apply.
    predicate_space <- if (!is.na(view$space)) {
        GiottoClass:::.resolve_space(gobject, view$space)
    } else NULL

    # ---- crop steps ---------------------------------------------------
    # A store can evaluate a crop on its own geometry when it HAS one and
    # that geometry is what the step is asking about: the cell polygon
    # for a cell-keyed store under `geom = "poly"`, or the points'
    # own geometry for a non-cell-keyed store (where a crop means
    # "clip these points", the same thing the in-memory path does).
    crop_steps <- GiottoClass:::.view_steps_of(view, "crop")
    can_push <- inherits(store, "parquetGeomBase")
    lazy_crops <- Filter(
        function(s) can_push && (!cell_keyed || identical(s$geom, "poly")),
        crop_steps)

    for (step in lazy_crops) {
        # `.materialize_crop_region` turns the recorded WKT back into a
        # SpatVector; spatRelate then takes it directly.
        y <- GiottoClass:::.materialize_crop_region(step$region)
        # Project region from predicate frame -> output frame so the
        # spat_relate eval (which runs against geom in the output frame,
        # via @post_ops) sees both sides in the same frame.
        y <- .project_region_between_spaces(y, gobject,
            from_space = predicate_space, to_space = space)
        store <- spatRelate(store, y, relation = step$relation)
    }

    # Crops the store could not push down become an eager cell_ID set.
    # Only meaningful for a cell-keyed target: on a points store a
    # cell_ID set does not address rows.
    if (length(crop_steps) > 0L) {
        if (cell_keyed) {
            eager_arms <- setdiff(c("centroid", "poly"),
                vapply(lazy_crops, `[[`, character(1L), "geom"))
            for (g in eager_arms) {
                store <- .queue_ids(store,
                    .view_crop_ids_arrow(view, gobject, .cache, coordinator,
                        geom = g, spat_unit = spat_unit))
            }
        } else if (!can_push) {
            warning("[push_view_to_pstore] crop step skipped: ",
                "target store is not parquetGeomBase (no geom ",
                "column to evaluate the predicate on)", call. = FALSE)
        }
    }

    # ---- filter steps -------------------------------------------------
    if (!is.null(.cache)) {
        return(.queue_ids(store,
            .view_filter_ids_arrow(view, gobject, .cache,
                spat_unit = spat_unit, feat_type = feat_type)))
    }
    for (step in GiottoClass:::.view_steps_of(view, "filter")) {
        store <- .narrow_store_by_predicate(
            target_store = store,
            predicate = str2lang(step$predicate),
            gobject = gobject,
            target_key = target_key,
            spat_unit = spat_unit,
            feat_type = feat_type
        )
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
# Affine matrix coercion: a transform step whose `op == "affine"` records
# `args = list(y = <matrix>)`. parquetGeomBase has an
# `(parquetGeomBase, affine2d)` method but NOT `(parquetGeomBase,
# matrix)` (dispatch would fail). Wrap the matrix in an affine2d before
# dispatch.
#
# Sample-key resolution: use the `:default:` sentinel for single-
# giotto contexts; the single key if just one is present; otherwise
# a no-op.
#
# This deliberately shadows nothing: GiottoClass has an internal of the
# same name, but it dispatches transforms on the subobject wrapper, which
# is exactly what a backed geometry must NOT do. Its signature takes a
# `coordinator`, so the eventual fix is for it to become a generic that
# parquetCoordinator can override; until then this stays local.
#
#' @keywords internal
#' @noRd
.apply_space_to_subobj <- function(subobj, gobject, space) {
    if (is.null(space)) return(subobj)
    keys <- names(space$samples)
    if (length(keys) == 0L) return(subobj)
    key <- if (":default:" %in% keys) {
        ":default:"
    } else if (length(keys) == 1L) {
        keys[[1L]]
    } else {
        return(subobj)
    }
    steps <- space$samples[[key]]
    backed_geom <- .hasSlot(subobj, "spatVector") &&
        inherits(subobj@spatVector, "parquetBase")
    for (step in steps) {
        args <- step$args
        if (backed_geom) {
            if (identical(step$op, "affine") &&
                inherits(args$y, "matrix")) {
                # ANY,missing affine method wraps a matrix into an affine2d
                args$y <- affine(args$y)
            }
            subobj@spatVector <- do.call(step$op,
                c(list(x = subobj@spatVector), args))
        } else {
            subobj <- do.call(step$op, c(list(x = subobj), args))
        }
    }
    subobj
}


# Walk a view's filter + crop steps and apply each to a data.table
# target via the arrow bridge (filters) or the crop-step id helper.
#
# A data.table target has no geometry of its own, so crops ALWAYS reduce
# to a cell_ID set here -- both `geom` arms, routed by
# `.view_crop_ids_arrow`. There is no lazy alternative to weigh, which
# is why this function has no `cell_keyed` argument: a non-cell-keyed DT
# target (featMeta, keyed by feat_ID) cannot be addressed by a cell_ID
# set at all, and says so.
#
#' @keywords internal
#' @noRd
.push_view_to_dt <- function(dt, view, gobject, coordinator,
    key = "cell_ID", spat_unit = NULL, feat_type = NULL, .cache = NULL) {
    cell_ID <- NULL # NSE
    if (is.null(view) || length(view$steps) == 0L) return(dt)
    if (is.null(dt) || nrow(dt) == 0L) return(dt)

    # Cache path: one semi_join against the view-wide intersection.
    # Only applies to cell_ID-keyed targets — current view steps narrow
    # on the cell axis. Non-cell keys (e.g. feat_metadata's feat_ID) fall
    # through to the no-cache path, which handles them correctly: crop
    # steps are skipped with a warning, filter steps go through the arrow
    # bridge keyed appropriately.
    if (!is.null(.cache) && identical(key, "cell_ID")) {
        surv <- .surviving_cell_ids_arrow(view, gobject, .cache,
            coordinator, spat_unit = spat_unit, feat_type = feat_type)
        if (is.null(surv)) return(dt)
        out <- dplyr::collect(
            dplyr::semi_join(arrow::arrow_table(dt), surv, by = "cell_ID")
        )
        return(data.table::setDT(out))
    }

    # No-cache path: per-step narrowing. Filter steps go through the
    # arrow-bridge predicate helper. Crop steps go through the same
    # `.view_crop_ids_arrow` the cache path uses, so the two agree by
    # construction rather than by two parallel implementations.
    crop_steps <- GiottoClass:::.view_steps_of(view, "crop")
    if (length(crop_steps) > 0L && !identical(key, "cell_ID")) {
        warning("[push_view_to_dt] crop steps skipped on non-cell-",
            "keyed target (key = '", key, "')", call. = FALSE)
        crop_steps <- list()
    }
    if (length(crop_steps) > 0L) {
        for (g in c("centroid", "poly")) {
            ids <- .view_crop_ids_arrow(view, gobject, .cache, coordinator,
                geom = g, spat_unit = spat_unit)
            if (is.null(ids)) next
            keep <- dplyr::collect(ids)$cell_ID
            dt <- dt[cell_ID %in% keep]
        }
    }
    for (step in GiottoClass:::.view_steps_of(view, "filter")) {
        dt <- .narrow_dt_via_arrow(dt, str2lang(step$predicate), gobject,
            key = key, spat_unit = spat_unit, feat_type = feat_type)
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
        if (is.null(view) || length(view$steps) == 0L) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        if (inherits(subobj@metaDT, "dataStore")) {
            subobj@metaDT <- .push_view_to_pstore(subobj@metaDT, view,
                gobject = gobject, coordinator = coordinator, space = space,
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        } else {
            subobj@metaDT <- .push_view_to_dt(subobj@metaDT, view, gobject,
                coordinator = coordinator, key = "cell_ID",
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
        if (is.null(view) || length(view$steps) == 0L) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        if (inherits(subobj@metaDT, "dataStore")) {
            subobj@metaDT <- .push_view_to_pstore(subobj@metaDT, view,
                gobject = gobject, coordinator = coordinator, space = space,
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        } else {
            subobj@metaDT <- .push_view_to_dt(subobj@metaDT, view, gobject,
                coordinator = coordinator, key = "feat_ID",
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
        # passed `space=`. `view$space` is NOT consulted here -- it only
        # affects the predicate frame inside the crop-step handlers.
        if (!is.null(space)) {
            subobj <- .apply_space_to_subobj(subobj, gobject, space)
        }
        if (is.null(view) || length(view$steps) == 0L) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        if (inherits(subobj@coordinates, "dataStore")) {
            subobj@coordinates <- .push_view_to_pstore(
                subobj@coordinates, view, gobject = gobject,
                coordinator = coordinator, space = space,
                spat_unit = spat_unit, .cache = .cache)
        } else {
            subobj@coordinates <- .push_view_to_dt(
                subobj@coordinates, view, gobject,
                coordinator = coordinator, key = "cell_ID",
                spat_unit = spat_unit, .cache = .cache)
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
        if (is.null(view) || length(view$steps) == 0L) return(subobj)
        if (is.null(subobj@enrichDT)) return(subobj)
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        if (inherits(subobj@enrichDT, "dataStore")) {
            subobj@enrichDT <- .push_view_to_pstore(subobj@enrichDT, view,
                gobject = gobject, coordinator = coordinator, space = space,
                spat_unit = spat_unit, feat_type = feat_type,
                .cache = .cache)
        } else {
            subobj@enrichDT <- .push_view_to_dt(subobj@enrichDT, view,
                gobject, coordinator = coordinator, key = "cell_ID",
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
        # `view$space` is consumed inside `.push_view_to_pstore` as the
        # predicate frame -- the crop region is interpreted there before
        # being projected into the output frame for the spat_relate eval.
        if (!is.null(space)) {
            subobj <- .apply_space_to_subobj(subobj, gobject, space)
        }
        if (is.null(view) || length(view$steps) == 0L) return(subobj)
        # A cell polygon store IS cell-keyed (poly_ID == cell_ID by
        # convention), so it consumes the same eager cell_ID set as
        # cellMeta / spatLocs / expression for the filter and centroid
        # arms -- and pushes `geom = "poly"` crops down onto its own geom
        # column, which is the one place in the pipeline where the exact
        # polygon predicate is free.
        #
        # `.cache` is passed straight through. It used to be forced to a
        # fresh environment here, which routed every crop through the
        # centroid answer and made `geom = "poly"` unreachable on a
        # backed polygon; the routing now comes off the step.
        subobj@spatVector <- .push_view_to_pstore(subobj@spatVector, view,
            gobject = gobject, coordinator = coordinator,
            target_key = "poly_ID", space = space,
            spat_unit = GiottoClass::spatUnit(subobj),
            .cache = list(...)$.cache)
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
        if (is.null(view) || length(view$steps) == 0L) return(subobj)
        # A crop reaches parity with in-mem via parquetGeomBase's
        # spat_relate op on the points' own geometry.
        #
        # `cell_keyed = FALSE` is the whole story for points: one store
        # row is one transcript, not one cell, so a surviving cell_ID set
        # does not address rows here and a crop means "clip these
        # points". This used to be expressed as `.cache = NULL`, which
        # got the right behaviour for the wrong reason and made the cache
        # look like a semantic switch.
        #
        # TODO: a filter step on giottoPoints is unreachable in real
        # workflows today — filters are cell-centric by convention (the
        # in-mem method skips them). Same-store predicates like
        # `feature_name == "GENE1"` work mechanically here but no view
        # consumer drives feature subsetting through a filter step; that
        # would belong to a future feature-select step. Cross-store
        # narrowing via cellMeta requires an aggregated points store
        # carrying cell_ID, which the current pipeline doesn't produce.
        # Code is in place for when either gap closes.
        subobj@spatVector <- .push_view_to_pstore(subobj@spatVector, view,
            gobject = gobject, coordinator = coordinator,
            target_key = "cell_ID",
            feat_type = GiottoClass::featType(subobj),
            cell_keyed = FALSE, .cache = NULL)
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
        if (is.null(view) || length(view$steps) == 0L) return(subobj)
        if (!inherits(subobj@exprMat, "parquetExprStore")) {
            # In-mem matrix / dgCMatrix / etc. — delegate to
            # dataTableCoordinator's exprObj method.
            return(callNextMethod())
        }
        .cache <- list(...)$.cache
        spat_unit <- GiottoClass::spatUnit(subobj)
        feat_type <- GiottoClass::featType(subobj)
        surv <- .surviving_cell_ids_arrow(view, gobject, .cache,
            coordinator, spat_unit = spat_unit, feat_type = feat_type)
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
