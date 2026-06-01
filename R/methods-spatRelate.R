#' @include class-parquetStore.R
NULL

# spatRelate on parquetGeomBase: lazy filter via spatial predicate ####
#
# Queues a "spat_relate" op carrying the query geometry (inline WKT, or a
# reference to another parquetGeomStore) and the predicate name. Evaluated
# at storeRead time by the SQL compile in `.pstore_to_sedona`. The arrow
# backend has no native spatial predicates; storeRead errors loudly on that
# path and directs callers to `output = "sedona"`. A tile-streaming arrow
# implementation is possible but not implemented.
#
# Phase 4a scope: filter form only (semi-join semantic -- narrow x by whether
# any feature of y satisfies the predicate). Store/store + form="join" is
# deferred to Phase 5.


# Predicate set (matches sf/sedona naming; terra's "coveredby" is normalized
# via `.terra_relation_name`).
.SPATRELATE_PREDICATES <- c(
    "intersects", "touches", "crosses", "overlaps",
    "within", "contains", "covers", "covered_by", "disjoint"
)

.validate_spatrelate_relation <- function(relation) {
    relation <- match.arg(relation, choices = .SPATRELATE_PREDICATES)
    relation
}

# terra uses "coveredby" (no underscore); normalize from user-facing name.
.terra_relation_name <- function(relation) {
    if (identical(relation, "covered_by")) "coveredby" else relation
}

# Choose the SRID literal for `ST_GeomFromText` to match the store's geom
# CRS. Sedonadb requires both sides of a spatial predicate to share a CRS.
# Biology data is typically written with no CRS (`crs = ""`); sedonadb's
# GeoParquet reader assigns the spec default (`ogc:crs84`, internally
# `lnglat()`) to absent or `null` CRS. `ST_GeomFromText('...', 4326)`
# produces a geometry with CRS `epsg:4326`, which sedonadb also normalizes
# to `lnglat()` -- so the two sides match. For stores written with an
# explicit `EPSG:NNNN` CRS, pull the SRID from the store params.
.spatrelate_store_srid <- function(store) {
    crs <- if (.hasSlot(store, "params")) store@params$crs else NULL
    if (is.null(crs) || !nzchar(crs)) return(4326L)
    if (identical(crs, "OGC:CRS84")) return(4326L)
    m <- regmatches(crs, regexec("EPSG:([0-9]+)", crs))[[1L]]
    if (length(m) >= 2L) return(as.integer(m[[2L]]))
    # Unrecognized CRS string -- emit without explicit SRID and surface
    # any downstream mismatch as a sedonadb error.
    NA_integer_
}

# sedona / DuckDB SQL function for a predicate.
.sql_relation_fn <- function(relation) {
    switch(relation,
        "intersects" = "ST_Intersects",
        "touches"    = "ST_Touches",
        "crosses"    = "ST_Crosses",
        "overlaps"   = "ST_Overlaps",
        "within"     = "ST_Within",
        "contains"   = "ST_Contains",
        "covers"     = "ST_Covers",
        "covered_by" = "ST_CoveredBy",
        "disjoint"   = "ST_Disjoint",
        stop(sprintf("[spat_relate] unknown relation: %s", relation), call. = FALSE)
    )
}


# Canonical entry: WKT string ####
#
# All other (parquetGeomBase, *) methods normalize to this one. Queues the op
# with the WKT inline and re-validates the predicate.

#' @name spatRelate
#' @rdname spatRelate
#' @inheritParams GiottoClass::spatRelate
#' @param y query geometry. For `parquetGeomBase` x: a single WKT string,
#'   a `SpatVector` (up to ~1000 features, see below), an `sf`/`sfc`, a
#'   `giottoPolygon`/`giottoPoints`/`spatLocsObj`, or another
#'   `parquetGeomBase` for the store/store path.
#' @details
#' The `parquetGeomBase` methods queue a lazy `"spat_relate"` op on `@ops`.
#' Evaluation happens at [storeRead()] time via the sedona path, which emits
#' `ST_<predicate>(geom, ...)` SQL against the parquet dataset. The arrow
#' backend has no native spatial predicates and errors loudly when a
#' `spat_relate` op is present — use `output = "sedona"` for spatial filters.
#'
#' Inline geometry inputs (`SpatVector`, `sf`, single WKT) carry the query
#' as a WKT string on the op. To keep ops compact and avoid pathological
#' serialization, the number of features in a `SpatVector`/`sf` query is
#' capped at `getOption("giottodisk.spatrelate_inline_max", 1000L)`. Larger
#' query sets should be written to a `parquetGeomStore` first; the store/
#' store path then handles them at scale.
#'
#' Single-feature inputs are passed through as-is. Multi-feature inputs are
#' unioned (`terra::aggregate` / `sf::st_union`) into one MULTIPOLYGON /
#' GEOMETRYCOLLECTION before serializing to WKT. Attributes are not carried
#' through (WKT is geometry-only); use the store/store path with
#' `form = "join"` (Phase 5) for attribute carry.
#' @export
setMethod(
    "spatRelate",
    signature(x = "parquetGeomBase", y = "character"),
    function(x, y, relation = "intersects", ...) {
        relation <- .validate_spatrelate_relation(relation)
        checkmate::assert_string(y, min.chars = 1L)
        x@ops <- c(x@ops, list(list(
            type     = "spat_relate",
            y_wkt    = y,
            y_store  = NULL,
            relation = relation,
            form     = "filter"
        )))
        x
    }
)

# SpatVector → WKT ####

#' @rdname spatRelate
#' @export
setMethod(
    "spatRelate",
    signature(x = "parquetGeomBase", y = "SpatVector"),
    function(x, y, relation = "intersects", ...) {
        n <- nrow(y)
        threshold <- getOption("giottodisk.spatrelate_inline_max", 1000L)
        if (n == 0L) {
            stop("[spatRelate] `y` has 0 features", call. = FALSE)
        }
        if (n > threshold) {
            stop(sprintf(
                "[spatRelate] inline SpatVector has %d features (max %d).\n  Convert to a parquetGeomStore via `storeWrite()` for the store/store path.",
                n, threshold
            ), call. = FALSE)
        }
        y_use <- if (n == 1L) y else terra::aggregate(y)
        wkt <- terra::geom(y_use, wkt = TRUE)
        spatRelate(x, wkt, relation = relation, ...)
    }
)

# sf / sfc → WKT ####

#' @rdname spatRelate
#' @export
setMethod(
    "spatRelate",
    signature(x = "parquetGeomBase", y = "sf"),
    function(x, y, relation = "intersects", ...) {
        GiottoUtils::package_check("sf")
        n <- nrow(y)
        threshold <- getOption("giottodisk.spatrelate_inline_max", 1000L)
        if (n == 0L) {
            stop("[spatRelate] `y` has 0 features", call. = FALSE)
        }
        if (n > threshold) {
            stop(sprintf(
                "[spatRelate] inline sf has %d features (max %d).\n  Convert to a parquetGeomStore via `storeWrite()` for the store/store path.",
                n, threshold
            ), call. = FALSE)
        }
        geom_col <- sf::st_geometry(y)
        wkt <- if (n == 1L) {
            sf::st_as_text(geom_col)
        } else {
            sf::st_as_text(sf::st_union(geom_col))
        }
        spatRelate(x, wkt, relation = relation, ...)
    }
)

# Unwrap giottoPolygon / giottoPoints (inner is SpatVector or parquetGeomStore) ####

#' @rdname spatRelate
#' @export
setMethod(
    "spatRelate",
    signature(x = "parquetGeomBase", y = "giottoPolygon"),
    function(x, y, relation = "intersects", ...) {
        spatRelate(x, y[], relation = relation, ...)
    }
)

#' @rdname spatRelate
#' @export
setMethod(
    "spatRelate",
    signature(x = "parquetGeomBase", y = "giottoPoints"),
    function(x, y, relation = "intersects", ...) {
        spatRelate(x, y[], relation = relation, ...)
    }
)

# spatLocsObj → points SpatVector → re-dispatch ####

#' @rdname spatRelate
#' @export
setMethod(
    "spatRelate",
    signature(x = "parquetGeomBase", y = "spatLocsObj"),
    function(x, y, relation = "intersects", ...) {
        # Promote points-as-data.table to a SpatVector of points so the
        # SpatVector method's size-check + WKT serialization can run.
        spatRelate(x, GiottoClass::as.points(y), relation = relation, ...)
    }
)

# Store/store: queue a reference to y's store ####

#' @rdname spatRelate
#' @param form `character(1)`. Either `"filter"` (default; semi-join semantic
#'   -- narrow x to features satisfying the predicate against any y) or
#'   `"join"` (bring in y's columns; not yet implemented).
#' @export
setMethod(
    "spatRelate",
    signature(x = "parquetGeomBase", y = "parquetGeomBase"),
    function(x, y, relation = "intersects",
             form = c("filter", "join"), ...) {
        relation <- .validate_spatrelate_relation(relation)
        form <- match.arg(form)
        if (form == "join") {
            stop("[spatRelate] form='join' (store/store join with attribute carry) is not yet implemented.",
                call. = FALSE)
        }
        x@ops <- c(x@ops, list(list(
            type     = "spat_relate",
            y_wkt    = NULL,
            y_store  = y,
            relation = relation,
            form     = form
        )))
        x
    }
)

# Rejection: in-memory x against stored y ####

#' @rdname spatRelate
#' @export
setMethod(
    "spatRelate",
    signature(x = "ANY", y = "parquetGeomBase"),
    function(x, y, ...) {
        stop("[spatRelate] in-memory x against stored y is unsupported.\n",
            "  Either materialize y first via `storeRead(y, output = \"terra\")`,\n",
            "  or wrap x in a parquetGeomStore for lazy evaluation.",
            call. = FALSE)
    }
)


# Evaluation: trim + sedona + collect ids ####
#
# Called by `.pbase_storeread_processing` when it hits a `spat_relate` op.
# Builds a trimmed copy of the store with @ops[1:i-1] (any prior spat_relate
# ops swapped for `id_filter` ops carrying cached ids), runs sedonadb on
# that to evaluate this op's predicate, and returns the surviving id rows
# as an arrow Table. The caller then `semi_join`s the arrow query with
# those ids -- arrow stays lazy past the spatial step.

.spat_relate_narrow <- function(store, i, cache) {
    GiottoUtils::package_check("sedonadb",
        repository = "github:apache/sedona-db/r/sedonadb")

    op <- store@ops[[i]]
    if (is.null(op$y_wkt) && !is.null(op$y_store)) {
        stop("[spat_relate] store-store form is not yet wired in the ",
            "narrow path; use `output = \"sedona\"` directly for now.",
            call. = FALSE)
    }

    # id cols: row_index is universal on parquetGeomBase; tile stores
    # additionally need tile_index since row_index resets per tile.
    id_cols <- if (inherits(store, "parquetGeomTileStore")) {
        c("row_index", "tile_index")
    } else {
        "row_index"
    }

    # Trim @ops to everything before this op. Replace any prior spat_relate
    # ops with id_filter ops carrying their cached ids, so the sedona compile
    # doesn't re-evaluate spatial predicates we already have answers for.
    prior_ops <- if (i > 1L) store@ops[seq_len(i - 1L)] else list()
    prior_ops <- lapply(seq_along(prior_ops), function(j) {
        op_j <- prior_ops[[j]]
        if (identical(op_j$type, "spat_relate")) {
            cached <- cache[[as.character(j)]]
            list(type = "id_filter", ids_tab = cached, by = names(cached))
        } else {
            op_j
        }
    })

    trim_store <- store
    trim_store@ops <- prior_ops
    # User-facing field narrowing (`[, j]` -> `store@fields`) must NOT
    # apply to the internal sedona evaluation -- we need `geom` available
    # to compute the spatial predicate regardless of what the caller asked
    # for in the final output projection. Set @fields to NULL so
    # `.pstore_lazy_fields` returns NULL and the underlying SELECT is `*`.
    if (.hasSlot(trim_store, "fields")) trim_store@fields <- NULL

    # Compile prior state to sedona, then layer this op's predicate via SQL.
    # Using sd_sql against the trim store's registered view is symmetric with
    # the existing `.pstore_to_sedona` SQL-building style and avoids the NSE
    # gymnastics of sd_filter with a constructed call.
    sdf <- storeRead(trim_store, output = "sedona")
    # `sdf` already reflects the trim_store's WHERE/SELECT/affine. Register
    # it under a fresh view name so the predicate we layer here applies on
    # top of the trim_store's filtered result, not the underlying raw
    # parquet union.
    view_name <- tolower(paste0("gd_srnarrow_", .make_uid()))
    sedonadb::sd_to_view(sdf, view_name, overwrite = TRUE)

    sql_pred <- .sql_relation_fn(op$relation)
    wkt_escaped <- gsub("'", "''", op$y_wkt, fixed = TRUE)
    srid <- .spatrelate_store_srid(store)
    geom_sql <- if (is.na(srid)) {
        sprintf("ST_GeomFromText('%s')", wkt_escaped)
    } else {
        sprintf("ST_GeomFromText('%s', %d)", wkt_escaped, srid)
    }
    select_sql <- paste(sprintf('"%s"', id_cols), collapse = ", ")
    sql <- sprintf(
        'SELECT %s FROM "%s" WHERE %s(geom, %s)',
        select_sql, view_name, sql_pred, geom_sql
    )
    ids_df <- sedonadb::sd_collect(sedonadb::sd_sql(sql))
    arrow::as_arrow_table(ids_df)
}
