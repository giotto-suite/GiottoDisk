#' @include class-viewCoordinator.R
NULL

# =============================================================================
# methods-resolveSubobject.R — parquetCoordinator dispatch
#
# `resolveKeep()` evaluates a view into the op's surviving cell set as a
# lazy arrow plan; the `resolveRecipe()` leaf methods queue it on each backed
# store as an `id_filter` and hand in-memory data to the inherited leaves.
# =============================================================================


# defaultViewCoordinator dispatch ####
# Any gsource subclass selects parquetCoordinator.

#' @rdname parquetCoordinator-class
#' @importFrom GiottoClass defaultViewCoordinator
#' @export
setMethod("defaultViewCoordinator", signature(source = "gsource"),
    function(source, ...) parquetCoordinator()
)


# Helpers ####

# The first subobject whose data covers all of `cols`, in spatValues'
# precedence: cell metadata, feature metadata, spatial locations, spatial
# enrichment, polygon attributes (key poly_ID), then expression by gene name.
# Returns `list(kind, source, key)`, or NULL.
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


# A wide data.table of the requested genes (`cell_ID` plus one column per
# gene), so a filter can reference genes as columns. Reads only those genes'
# rows: bounded by n_cells x genes requested.
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

# Centroids for the crop steps, in the predicate frame. Only the fetch is
# local -- a backed `@coordinates` must be read first -- and the predicate is
# GiottoClass's `spatRelate()`, so both backends answer the same question.
# On a multi, cell_IDs are prefixed `<sample>::` to match the joint
# vocabulary.
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

# The polygon source for a `geom = "poly"` crop, in the predicate frame, as
# a giottoPolygon so a backed `@spatVector` brings its own engines. On a
# multi, poly_IDs are prefixed to match the joint vocabulary.
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

# One crop step -> the cell_IDs that survive it. `geom` picks the carrier
# (centroids or polygons); `spatRelate()` evaluates it on whichever engine
# the carrier has.
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

# Carriers for the crop arms, in the frame each step names (independent of
# the output space). Built lazily and memoised per frame; a NULL build ("no
# spatial locations") is cached too, so its warning fires once per frame.
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

# The view's surviving cell set as one lazy arrow plan: the steps, in order,
# semi-joined onto a running result. Nothing is read here; the plan runs
# inside each target's read. NULL means unconstrained, not empty. Cell axis
# only (adr/0015).
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
# `arrow` is the plan backed leaves queue; `vector` is what an in-memory
# subobject falling through to the dataTableCoordinator leaves reads. It can
# only be had by running the plan, so it is a promise.

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


# Apply a space's transform steps to a subobject. A backed polygon / points
# store is transformed directly (its methods compose into `@post_ops`),
# because the wrapper methods call terra; anything else dispatches on the
# subobject. An `affine` step's matrix is wrapped in an affine2d, the only
# form the store method takes.
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
# Each method narrows a backed slot and hands in-memory data to the
# inherited leaf. featMetaObj and dimObj inherit the in-memory leaf outright.

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

# Points are not cell-keyed, so a crop clips the store's own geometry in the
# output frame, projected from the frame its step names. Filter steps are
# skipped, as in memory.

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
