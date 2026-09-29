#' @include class-viewCoordinator.R
NULL

# =============================================================================
# methods-resolveSubobject.R — parquetCoordinator dispatch
#
# Registers:
#   1. defaultViewCoordinator on gsource, so a gobject with a
#      gsource-inheriting `@source` selects parquetCoordinator.
#   2. resolveKeep for parquetCoordinator: a view evaluated into the op's
#      surviving cell set as a LAZY arrow query -- the accumulated narrowing,
#      run only when a target is read.
#   3. resolve leaf methods on the subobject classes whose data slot may hold
#      a parquetStore-inheriting backing. Each queues the `arrow` form of the
#      set as an `id_filter`; in-memory data falls through to GiottoClass's
#      dataTableCoordinator leaves, which read the `vector` form.
#
# A filter's predicate columns may live on a different subobject than the
# target being narrowed, and either side can be backed or in memory.
# `.find_store_with_cols()` finds the owner; a backed owner narrows lazily
# with `subset()` and contributes its key column as a query, an in-memory
# owner narrows eagerly and contributes an arrow Table. Crop steps reduce to
# a cell_ID set through `.crop_step_ids()`. The steps are chained with
# semi-joins, so the whole view is one plan that nothing executes until a
# target store is collected -- on the arrow engine inside that target's own
# query, on duckdb / sedona as a registered id table.
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
#   5. polygon attributes(giottoPolygon@spatVector, key = poly_ID)
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

    # Polygon attributes. poly_ID aligns with the cells spat_unit's cell_ID
    # by convention, so a polygon store can answer for the cell axis.
    gp <- tryCatch(GiottoClass::getPolygonInfo(gobject, name = spat_unit,
        return_giottoPolygon = TRUE, verbose = FALSE),
        error = function(e) NULL)
    if (inherits(gp, "giottoPolygon")) {
        sv <- gp@spatVector
        src <- if (inherits(sv, "SpatVector")) {
            data.table::as.data.table(terra::values(sv))
        } else sv
        hit <- .check(src, "poly_ID")
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

# One filter step -> the owner's surviving key column, still lazy when the
# owner is backed. Returns `list(ids, key)`, where `ids` is an arrow query
# (backed owner, nothing read yet) or an arrow Table (in-memory owner, which
# is already in memory, so narrowing it eagerly costs no I/O).
#
#' @keywords internal
#' @noRd
.narrowed_ids <- function(predicate, gobject, spat_unit = NULL,
                          feat_type = NULL) {
    cols <- all.vars(predicate)
    owner <- .find_store_with_cols(gobject, cols,
        spat_unit = spat_unit, feat_type = feat_type)
    if (is.null(owner)) {
        stop("[resolveKeep] no subobject covers predicate cols: ",
            paste(cols, collapse = ", "), call. = FALSE)
    }
    key <- owner$key
    ids <- if (identical(owner$kind, "parquetBase")) {
        narrowed <- subset(owner$source, predicate, quote = FALSE)
        # `[, key]` narrows the read; a geometry store still carries its
        # index columns, so the key is selected explicitly
        q <- storeRead(narrowed[, key], output = "query")
        dplyr::distinct(dplyr::select(q, dplyr::all_of(key)))
    } else {
        mask <- eval(predicate, envir = owner$source, enclos = baseenv())
        arrow::arrow_table(as.data.frame(
            unique(owner$source[mask, key, with = FALSE])))
    }
    # a polygon owner speaks the cell vocabulary under another column name
    if (identical(key, "poly_ID")) {
        ids <- dplyr::rename(ids, cell_ID = "poly_ID")
        key <- "cell_ID"
    }
    list(ids = ids, key = key)
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
.projected_points <- function(gobject, space, spat_unit = NULL) {
    cell_ID <- NULL # NSE
    sl <- tryCatch(GiottoClass::getSpatialLocations(gobject,
        output = "spatLocsObj", spat_unit = spat_unit),
        error = function(e) NULL)
    if (is.null(sl)) return(NULL)

    # A backed `@coordinates` has to be read before the geometry can be
    # built; everything after that is the same points `SpatVector` the
    # in-memory path uses, so the predicate itself is shared code.
    .materialize <- function(x) {
        co <- x@coordinates
        if (inherits(co, "dataStore")) {
            x@coordinates <- data.table::setDT(
                storeRead(co, output = "tibble"))
        }
        x
    }

    if (inherits(sl, "spatLocsObj")) {
        return(GiottoClass::as.points(.materialize(
            .apply_space_to_subobj(sl, space))))
    }
    if (!is.list(sl)) {
        stop("[.projected_points] unexpected spatial locations type: ",
            toString(class(sl)), call. = FALSE)
    }
    parts <- lapply(names(sl), function(nm) {
        child <- sl[[nm]]
        if (!inherits(child, "spatLocsObj")) return(NULL)
        # `[` owns the sample-resolution rule
        child <- .materialize(.apply_space_to_subobj(child,
            space[nm]))
        dt <- data.table::copy(child[])
        dt[, cell_ID := paste(nm, cell_ID, sep = "::")]
        child[] <- dt
        child
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) return(NULL)
    # fold as spatLocsObjs -- a data.table rbind -- and convert once, so
    # there is one terra allocation rather than one per child
    GiottoClass::as.points(Reduce(rbind2, parts))
}

# Fetch the polygon source in the predicate frame, for a `geom = "poly"`
# crop. Returns a `giottoPolygon`, so `spatRelate()` dispatches on it and
# a backed `@spatVector` brings its own engines.
#
# giottoMulti: polygons live per child, so each child's are fetched,
# space-scoped, and their poly_IDs prefixed to match the joint cell
# vocabulary the resulting ID set is intersected against.
#' @keywords internal
#' @noRd
.projected_polys <- function(gobject, space, spat_unit = NULL) {
    one <- function(g, samp) {
        gp <- tryCatch(GiottoClass::getPolygonInfo(g, name = spat_unit,
            return_giottoPolygon = TRUE, verbose = FALSE),
            error = function(e) NULL)
        if (!inherits(gp, "giottoPolygon")) return(NULL)
        .apply_space_to_subobj(gp,
            if (is.null(space) || is.na(samp)) space else space[samp])
    }

    if (!inherits(gobject, "giottoMulti")) return(one(gobject, NA_character_))

    parts <- lapply(names(gobject@objects), function(nm) {
        gp <- one(gobject@objects[[nm]], nm)
        if (is.null(gp)) return(NULL)
        sv <- gp@spatVector
        sv$poly_ID <- paste(nm, terra::values(sv)$poly_ID, sep = "::")
        gp@spatVector <- sv
        gp@unique_ID_cache <- terra::values(sv)$poly_ID
        gp
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) return(NULL)
    if (length(parts) == 1L) return(parts[[1L]])
    do.call(rbind, parts)
}

# One crop step -> the cell_IDs that survive it.
#
# Both arms are one public expression: narrow a carrier with
# `spatRelate()`, then read its IDs. `geom` picks WHICH carrier -- the
# cells' centroids or their polygons -- and nothing else. The polygon
# carrier may be backed, in which case `spatRelate()` dispatches to the
# store's own sedona / duckdb / terra engines.
#' @keywords internal
#' @noRd
.crop_step_ids <- function(gobject, step, carriers) {
    region <- terra::vect(step$region)
    switch(step$geom,
        centroid = {
            pts <- carriers$points(step$space)
            if (is.null(pts)) return(NULL)
            GiottoClass::spatRelate(pts, region,
                relation = step$relation)$cell_ID
        },
        poly = {
            polys <- carriers$polys(step$space)
            if (is.null(polys)) {
                stop(sprintf(paste0(
                    "[crop] geom = \"poly\" was requested (relation '%s'), ",
                    "but this object has no polygon source to evaluate it ",
                    "on.\nEither add polygons (`setPolygonInfo()`) or use ",
                    "geom = \"centroid\"."), step$relation), call. = FALSE)
            }
            GiottoClass::spatIDs(GiottoClass::spatRelate(polys, region,
                relation = step$relation))
        },
        stop("[crop] unknown geom '", step$geom, "'", call. = FALSE)
    )
}

#' @keywords internal
#' @noRd
.ids_arrow <- function(ids) {
    arrow::arrow_table(data.frame(cell_ID = ids, stringsAsFactors = FALSE))
}

# Carriers for the crop arms, plus the frame each one is projected into.
#
# The predicate frame is named BY THE STEP -- the frame that step's region
# coordinates were read in. It is independent of any OUTPUT space the
# caller asked for; that one is applied by `.apply_space_to_subobj` on the
# subobject itself. Conflating the two is what made a `space =`-bound view
# return rotated coordinates from a plain getter.
#
# Built lazily and memoised PER FRAME: two crop steps in one view may name
# different spaces, and steps sharing a frame -- the common case -- share
# one build. A `NULL` build is a real answer ("no spatial locations"), so
# it is cached too and the warning fires once per frame, not once per step.
#' @keywords internal
#' @noRd
.crop_carriers <- function(gobject, spat_unit = NULL) {
    memo <- new.env(parent = emptyenv())
    memoised <- function(kind, space_name, build) {
        key <- paste0(kind, ":", if (is.na(space_name)) "" else space_name)
        if (!exists(key, envir = memo, inherits = FALSE)) {
            assign(key, build(.step_space(gobject, space_name)),
                envir = memo)
        }
        base::get(key, envir = memo)
    }
    list(
        space = function(space_name) .step_space(gobject, space_name),
        points = function(space_name) {
            out <- memoised("pts", space_name, function(sp) {
                .projected_points(gobject, sp, spat_unit = spat_unit)
            })
            if (is.null(out)) {
                warning("[view] crop step skipped: no spatial locations ",
                    "available", call. = FALSE)
            }
            out
        },
        polys = function(space_name) {
            memoised("poly", space_name, function(sp) {
                .projected_polys(gobject, sp, spat_unit = spat_unit)
            })
        }
    )
}

#' @keywords internal
#' @noRd
.step_space <- function(gobject, space_name) {
    if (is.null(space_name) || is.na(space_name)) return(NULL)
    GiottoClass::giottoSpace(gobject, space_name)
}

# The whole cell-axis answer for a view, as ONE lazy arrow plan: walk the
# recorded steps IN ORDER and semi-join what each leaves standing onto the
# running result. Nothing is read here -- a backed filter owner contributes a
# query -- so the plan runs once per target that is collected, inside that
# target's own read.
#
# Ordered, not gathered by kind. Order is information -- a read-time
# collapse needs it, and gathering by type discards it.
#
# NULL means "unconstrained", not "empty". Errors if a filter resolves
# against a non-`cell_ID` key -- multi-key intersection is out of scope,
# matching the in-memory coordinator. v1 models the CELL axis only; feature
# and subcellular-point axes are deferred -- see adr/0015.
#' @keywords internal
#' @noRd
.surviving_ids_query <- function(view, gobject, spat_unit = NULL,
                                 feat_type = NULL) {
    carriers <- .crop_carriers(gobject, spat_unit = spat_unit)
    surviving <- NULL
    for (step in view@steps) {
        ids <- switch(step$type,
            filter = {
                # Q7 records the predicate deparsed, so it comes back as a
                # string and is re-parsed before `all.vars()` / `eval()`.
                narrowed <- .narrowed_ids(str2lang(step$predicate), gobject,
                    spat_unit = spat_unit, feat_type = feat_type)
                if (!identical(narrowed$key, "cell_ID")) {
                    stop("[resolveKeep] only cell_ID-keyed filter steps ",
                        "are supported (got '", narrowed$key, "')",
                        call. = FALSE)
                }
                narrowed$ids
            },
            crop = {
                got <- .crop_step_ids(gobject, step, carriers)
                if (is.null(got)) NULL else .ids_arrow(got)
            },
            # `samples` narrows children, not cells; resolved before any
            # per-child work reaches here.
            NULL
        )
        if (is.null(ids)) next
        surviving <- if (is.null(surviving)) ids else
            dplyr::semi_join(surviving, ids, by = "cell_ID")
    }
    surviving
}


# resolveKeep (parquetCoordinator) ####
#
# The set carries both forms: `arrow`, the lazy plan every backed leaf queues
# as its `id_filter`, and `vector`, which the dataTableCoordinator leaves read
# when an in-memory subobject inside a backed gobject falls through to them.
# `vector` can only be had by running the plan, so it is installed as a
# promise and collected only if such a leaf reads it.

#' @rdname parquetCoordinator-class
#' @importFrom GiottoClass resolveKeep
#' @export
setMethod("resolveKeep", signature(coordinator = "parquetCoordinator"),
    function(coordinator, gobject, view, spat_unit = NULL, feat_type = NULL,
             ...) {
        if (is.null(view) || length(view@steps) == 0L) return(NULL)
        q <- .surviving_ids_query(view, gobject, spat_unit = spat_unit,
            feat_type = feat_type)
        if (is.null(q)) return(NULL)
        keep <- structure(list2env(list(arrow = q), parent = emptyenv()),
            class = "viewKeep")
        delayedAssign("vector", dplyr::collect(q)$cell_ID,
            assign.env = keep)
        keep
    }
)


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
# is exactly what a backed geometry must NOT do.
#
#' @keywords internal
#' @noRd
.apply_space_to_subobj <- function(subobj, space) {
    if (is.null(space)) return(subobj)
    # `[[` owns sample resolution. NA means "no sample
    # identity", which is what a single `giotto` -- or an already-scoped
    # handle from `space[<child>]` -- presents.
    steps <- space[[NA_character_]]
    if (length(steps) == 0L) return(subobj)
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


# resolve leaf methods (parquetCoordinator) ####
#
# Each method narrows a BACKED slot and hands anything in memory to the
# inherited dataTableCoordinator leaf with `callNextMethod()`, which reads
# the `vector` form. featMetaObj and dimObj register nothing here: feature
# metadata is not cell-keyed, and no pipeline produces backed dim-reduction
# coordinates, so both inherit the in-memory leaf outright.

# Queue the op's surviving set on a backed store. `by` maps the store's key
# column onto the set's `cell_ID`.
#' @keywords internal
#' @noRd
.queue_keep <- function(store, keep, by = "cell_ID") {
    if (is.null(keep)) return(store)
    store@ops <- c(store@ops, list(list(
        type = "id_filter", ids_tab = keep$arrow, by = by)))
    store
}

#' @rdname parquetCoordinator-class
#' @importFrom GiottoClass resolve
#' @export
setMethod("resolveRecipe",
    signature(x = "cellMetaObj", coordinator = "parquetCoordinator"),
    function(x, coordinator, keep = NULL, space = NULL, view = NULL, ...) {
        if (!inherits(x@metaDT, "dataStore")) return(callNextMethod())
        x@metaDT <- .queue_keep(x@metaDT, keep)
        x
    }
)

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveRecipe",
    signature(x = "spatEnrObj", coordinator = "parquetCoordinator"),
    function(x, coordinator, keep = NULL, space = NULL, view = NULL, ...) {
        if (!inherits(x@enrichDT, "dataStore")) return(callNextMethod())
        x@enrichDT <- .queue_keep(x@enrichDT, keep)
        x
    }
)

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveRecipe",
    signature(x = "spatLocsObj", coordinator = "parquetCoordinator"),
    function(x, coordinator, keep = NULL, space = NULL, view = NULL, ...) {
        if (!inherits(x@coordinates, "dataStore")) return(callNextMethod())
        x <- .apply_space_to_subobj(x, space)
        x@coordinates <- .queue_keep(x@coordinates, keep)
        x
    }
)

# A cell polygon store IS cell-keyed: poly_ID aligns with the cells
# spat_unit's cell_ID by convention, so it narrows by the same set, mapping
# poly_ID onto cell_ID.

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveRecipe",
    signature(x = "giottoPolygon", coordinator = "parquetCoordinator"),
    function(x, coordinator, keep = NULL, space = NULL, view = NULL, ...) {
        if (!inherits(x@spatVector, "parquetBase")) return(callNextMethod())
        x <- .apply_space_to_subobj(x, space)
        x@spatVector <- .queue_keep(x@spatVector, keep,
            by = c(poly_ID = "cell_ID"))
        x
    }
)

# Points are not cell-keyed: one row is one transcript, so a cell_ID set
# does not address rows, and a crop means "clip these points" -- a
# `spat_relate` op on the store's own geometry, in the output frame. The
# region is projected there from the frame its step was drawn in, looked up
# in `spaces`. Filter steps are skipped, as the in-memory leaf skips them:
# filters are cell-centric, and a feature-select step would be its own axis.

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveRecipe",
    signature(x = "giottoPoints", coordinator = "parquetCoordinator"),
    function(x, coordinator, keep = NULL, space = NULL, view = NULL,
             spaces = NULL, ...) {
        if (!inherits(x@spatVector, "parquetBase")) return(callNextMethod())
        x <- .apply_space_to_subobj(x, space)
        for (step in if (is.null(view)) list() else view@steps) {
            if (!identical(step$type, "crop")) next
            y <- GiottoClass::project_region(step$region,
                from_space = .step_frame(step, spaces), to_space = space)
            x@spatVector <- spatRelate(x@spatVector, y,
                relation = step$relation)
        }
        x
    }
)

# The frame a crop step's region was drawn in, or NULL for the native frame.
#' @keywords internal
#' @noRd
.step_frame <- function(step, spaces) {
    nm <- step$space
    if (is.null(nm) || is.na(nm)) return(NULL)
    sp <- spaces[[nm]]
    if (is.null(sp)) {
        stop(sprintf(paste0("[resolve] a crop step was drawn in space ",
            "'%s', which is not registered on this object"), nm),
            call. = FALSE)
    }
    sp
}

# A backed expression store narrows through `[, j]`, which updates
# `@cell_idx` lazily -- no rewrite, no I/O -- but takes a character vector,
# so it reads the `vector` form and collects the plan once.

#' @rdname parquetCoordinator-class
#' @export
setMethod("resolveRecipe",
    signature(x = "exprObj", coordinator = "parquetCoordinator"),
    function(x, coordinator, keep = NULL, space = NULL, view = NULL, ...) {
        if (!inherits(x@exprMat, "parquetExprStore")) return(callNextMethod())
        if (is.null(keep)) return(x)
        ids <- intersect(keep$vector, x@exprMat@cell_ids)
        x@exprMat <- if (length(ids) == 0L) {
            x@exprMat[, integer(0L)]
        } else {
            x@exprMat[, ids]
        }
        x
    }
)
