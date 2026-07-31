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

# Accepted values for the `engine` arg / option. NULL passes through
# (caller-supplied default; resolver handles option + auto fallback).
.SPATRELATE_ENGINES <- c("sedona", "duckdb", "terra", "auto")
.validate_spatrelate_engine <- function(engine) {
    if (is.null(engine)) return(NULL)
    match.arg(engine, choices = .SPATRELATE_ENGINES)
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
#' @param engine `character(1)` or `NULL`. Spatial-query engine to use for
#'   the predicate evaluation. One of `"sedona"`, `"duckdb"`, `"terra"`, or
#'   `"auto"`. `NULL` (default) defers to the package option (see Details).
#' @details
#' The `parquetGeomBase` methods queue a lazy `"spat_relate"` op on `@ops`.
#' At [storeRead()] time the op is evaluated by one of three engines:
#'
#' * `"sedona"` — emits `ST_<predicate>(geom, ...)` SQL against parquet via
#'   SedonaDB; lazy, parquet-stats pushdown applies. Requires `{sedonadb}`.
#' * `"duckdb"` — same shape via DuckDB's spatial extension; lazy. Requires
#'   `{duckdb}` + `{dbplyr}`.
#' * `"terra"` — deps-free fallback. Tile stores stream per-tile via
#'   `tilework::tileApply`; non-tile stores materialize the trim and run
#'   `terra::relate`. Slower than the SQL engines at scale but always
#'   available since `terra` is a hard import.
#'
#' Engine selection precedence (highest first):
#'
#' 1. `engine` arg on the `spatRelate()` call. Wins over the option.
#' 2. `getOption("giottodisk.spatial_query_engine")` — set globally or
#'    via `GiottoUtils::gwith_options(...)`.
#' 3. Default `"auto"`: picks the best available — sedona > duckdb > terra.
#'    If auto falls through to terra (no SQL engine installed), a one-shot
#'    `rlang::inform` message suggests installing `{sedonadb}` or
#'    `{duckdb}` for performance.
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
    function(x, y, relation = "intersects", engine = NULL, ...) {
        relation <- .validate_spatrelate_relation(relation)
        checkmate::assert_string(y, min.chars = 1L)
        engine <- .validate_spatrelate_engine(engine)
        x@ops <- c(x@ops, list(list(
            type     = "spat_relate",
            y_wkt    = y,
            y_store  = NULL,
            relation = relation,
            form     = "filter",
            engine   = engine
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
             form = c("filter", "join"), engine = NULL, ...) {
        relation <- .validate_spatrelate_relation(relation)
        form <- match.arg(form)
        engine <- .validate_spatrelate_engine(engine)
        if (form == "join") {
            stop("[spatRelate] form='join' (store/store join with attribute carry) is not yet implemented.",
                call. = FALSE)
        }
        x@ops <- c(x@ops, list(list(
            type     = "spat_relate",
            y_wkt    = NULL,
            y_store  = y,
            relation = relation,
            form     = form,
            engine   = engine
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


# Evaluation: trim + run engine + collect ids ####
#
# Called by `.pbase_storeread_processing` when it hits a `spat_relate` op.
# Builds a trimmed copy of the store with @ops[1:i-1] (any prior spat_relate
# ops swapped for `id_filter` ops carrying cached ids), runs the selected
# engine on that to evaluate this op's predicate, and returns the surviving
# id rows as an arrow Table. The caller then `semi_join`s the arrow query
# with those ids -- arrow stays lazy past the spatial step.
#
# Engine is resolved per call:
#   1. If `getOption("giottodisk.spatial_query_engine")` is set to a
#      specific engine, use it.
#   2. If set to "auto" (or unset), pick the best available:
#      sedona > duckdb > terra.
# Tests can pin a specific engine via `withr::with_options(...)`.

.spat_relate_narrow <- function(store, i, cache) {
    op <- store@ops[[i]]
    if (is.null(op$y_wkt) && !is.null(op$y_store)) {
        stop("[spat_relate] store-store form is not yet wired in the ",
            "narrow path; use `output = \"sedona\"` directly for now.",
            call. = FALSE)
    }

    # id cols: source_id is included as a safety backstop for any future
    # union / multi-source geom store (it's free -- the SQL engines inject
    # it as a per-tile literal, arrow auto-promotes from the hive
    # partition). row_index is the universal per-store identifier; tile
    # stores additionally need tile_index since row_index resets per tile.
    id_cols <- if (inherits(store, "parquetGeomTileStore")) {
        c("source_id", "row_index", "tile_index")
    } else {
        c("source_id", "row_index")
    }

    # Trim @ops to everything before this op. Replace any prior spat_relate
    # ops with id_filter ops carrying their cached ids, so the engine does
    # not re-evaluate spatial predicates we already have answers for.
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
    # apply to the internal evaluation -- we need `geom` available to
    # compute the spatial predicate regardless of what the caller asked
    # for in the final output projection. Set @fields to NULL so
    # `.pstore_lazy_fields` returns NULL and the underlying SELECT is `*`.
    if (.hasSlot(trim_store, "fields")) trim_store@fields <- NULL

    engine <- .resolve_spat_relate_engine(op$engine)
    switch(engine,
        "sedona" = .spat_relate_narrow_sedona(trim_store, op, id_cols, store),
        "duckdb" = .spat_relate_narrow_duckdb(trim_store, op, id_cols, store),
        "terra"  = .spat_relate_narrow_terra(trim_store, op, id_cols, store),
        stop(sprintf("[spat_relate] unknown engine: '%s'", engine), call. = FALSE)
    )
}


# Internal indirection so tests can mock the "is this engine installed?"
# check without touching base::requireNamespace.
.spat_engine_available <- function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
}


# Resolve the spatial-query engine to use for spat_relate narrowing.
# Honors (in precedence order):
#   1. `op_engine` — per-call engine passed to `spatRelate(..., engine = ...)`
#   2. `giottodisk.spatial_query_engine` option — user / session default
#   3. "auto" — best available (sedona > duckdb > terra)
# When "auto" falls through to terra (because no SQL spatial backend is
# installed) we emit a one-shot `rlang::inform` nudging the user toward a
# faster engine -- but only when the user hasn't already made a
# deliberate engine choice via the arg or the option.
.resolve_spat_relate_engine <- function(op_engine = NULL) {
    # Per-call engine arg wins over the option.
    eff <- op_engine %||% getOption("giottodisk.spatial_query_engine", "auto")
    if (!identical(eff, "auto")) return(eff)
    if (.spat_engine_available("sedonadb")) return("sedona")
    if (.spat_engine_available("duckdb"))   return("duckdb")
    rlang::inform(
        paste0(
            "spat_relate is using the terra engine (the default fallback ",
            "when no SQL spatial backend is installed). For ad-hoc ",
            "spatial queries, sedonadb or duckdb are typically faster. ",
            "Install one and set ",
            "`options(giottodisk.spatial_query_engine = \"sedona\")` or ",
            "`\"duckdb\"` to silence this message."
        ),
        .frequency = "once",
        .frequency_id = "giottodisk.spat_relate_terra_nudge"
    )
    "terra"
}


# Sedona engine for spat_relate narrowing.
# Compiles trim_store via `.pstore_to_sedona`, registers the resulting
# sdf as a view, layers this op's predicate via raw SQL, collects only
# the id cols.
.spat_relate_narrow_sedona <- function(trim_store, op, id_cols, store) {
    GiottoUtils::package_check("sedonadb",
        repository = "github:apache/sedona-db/r/sedonadb")
    sdf <- storeRead(trim_store, output = "sedona")
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


# DuckDB engine for spat_relate narrowing.
# Mirrors the sedona path: compile trim_store via `.pstore_to_duckdb`
# (which already returns a tbl_dbi over the trim state), then layer this
# op's predicate via raw SQL against the underlying view name. DuckDB's
# spatial extension's `ST_GeomFromText` takes (VARCHAR [, BOOLEAN]) and
# has no CRS concept, so SRID is omitted.
.spat_relate_narrow_duckdb <- function(trim_store, op, id_cols, store) {
    GiottoUtils::package_check("duckdb")
    GiottoUtils::package_check("dbplyr")
    tbl <- storeRead(trim_store, output = "duckdb")
    conn <- tbl$src$con
    inner_view <- dbplyr::remote_name(tbl)
    sql_pred <- .sql_relation_fn(op$relation)
    wkt_escaped <- gsub("'", "''", op$y_wkt, fixed = TRUE)
    geom_sql <- sprintf("ST_GeomFromText('%s')", wkt_escaped)
    select_sql <- paste(sprintf('"%s"', id_cols), collapse = ", ")
    sql <- sprintf(
        'SELECT %s FROM "%s" WHERE %s(geom, %s)',
        select_sql, inner_view, sql_pred, geom_sql
    )
    ids_df <- DBI::dbGetQuery(conn, sql)
    arrow::as_arrow_table(ids_df)
}


# Terra engine for spat_relate narrowing (deps-free fallback).
# Materializes the trim store's id_cols + geom via the arrow-output path,
# decodes WKB to a SpatVector, runs `terra::relate` against the query
# geom, returns surviving id rows.
#
# This is the deps-free path -- terra is a hard import of GiottoDisk, so
# this works on any install. The trade-off is that the whole trim is
# pulled into memory (sedona / duckdb stream parquet directly). For
# tile-store workloads at atlas scale, prefer a SQL backend; for ad-hoc
# in-memory queries and dev workflows the materialization cost is fine.
.spat_relate_narrow_terra <- function(trim_store, op, id_cols, store) {
    relation <- .terra_relation_name(op$relation)
    y_wkt <- op$y_wkt

    # `trim_store@fields` was set to NULL by the caller so the SQL
    # engines see all on-disk cols; we likewise need the partition /
    # row_index cols visible here. `omit_internals = FALSE` keeps them.

    # Per-tile (or per-store for non-tile) predicate evaluation. WKB
    # decode via `terra::vect(as.list(geom))`: `as.list` strips the
    # arrow_binary class, terra accepts WKB list-of-raw natively, no
    # `wk` dep. Faster than going through output = "terra" since we
    # skip the `terra::values(sv) <- data` attribute attachment that
    # `.pgstore_to_spatial` does (we filter on geom only and pluck
    # id_cols from the source tibble).
    eval_chunk <- function(df) {
        if (is.null(df) || nrow(df) == 0L) return(NULL)
        if (!"geom" %in% names(df)) {
            stop("[spat_relate][terra] expected 'geom' col missing from chunk",
                call. = FALSE)
        }
        x_sv <- terra::vect(as.list(df$geom))
        y_sv <- terra::vect(y_wkt)
        rel <- terra::relate(x_sv, y_sv, relation = relation)
        keep <- if (is.matrix(rel)) {
            as.logical(rowSums(rel, na.rm = TRUE) > 0L)
        } else {
            as.logical(rel)
        }
        if (!any(keep)) return(NULL)
        df[keep, id_cols, drop = FALSE]
    }

    if (inherits(store, "parquetGeomTileStore")) {
        # Tile-stream via tilework::tileApply -- bounded memory per tile.
        # tileApply iterates tiles, calls storeRead per tile with
        # get_params_x, and passes the per-tile data to FUN. trim_store's
        # @ops apply per-tile because storeRead applies them at each call.
        # Per-tile reads don't carry `source_id` (that's a hive partition
        # col added by the store-level read path); we inject it after the
        # eval since it's a known constant for the substore.
        tile_id_cols <- setdiff(id_cols, "source_id")
        eval_tile <- function(df) {
            if (is.null(df) || nrow(df) == 0L) return(NULL)
            if (!"geom" %in% names(df)) {
                stop("[spat_relate][terra] expected 'geom' col missing from tile",
                    call. = FALSE)
            }
            x_sv <- terra::vect(as.list(df$geom))
            y_sv <- terra::vect(y_wkt)
            rel <- terra::relate(x_sv, y_sv, relation = relation)
            keep <- if (is.matrix(rel)) {
                as.logical(rowSums(rel, na.rm = TRUE) > 0L)
            } else {
                as.logical(rel)
            }
            if (!any(keep)) return(NULL)
            df[keep, tile_id_cols, drop = FALSE]
        }
        results <- tilework::tileApply(
            trim_store,
            tiles = trim_store@tiles,
            FUN = function(tile_data, ...) eval_tile(tile_data),
            get_params_x = list(output = "tibble", omit_internals = FALSE)
        )
        results <- results[!vapply(results, is.null, logical(1L))]
        if (length(results) == 0L) {
            empty <- as.data.frame(
                matrix(integer(0L), ncol = length(id_cols),
                    dimnames = list(NULL, id_cols)))
            return(arrow::as_arrow_table(empty))
        }
        # `data.table::rbindlist` + a single arrow conversion is faster
        # than building per-tile arrow Tables and `concat_tables`-ing
        # them: per-tile arrow has fixed per-conversion overhead that
        # dominates for many-tile workloads (benched 5x-275x faster
        # across n_tiles=10..1000, k_per_tile=100..100k).
        combined <- data.table::rbindlist(results)
        combined[, source_id := trim_store@uid]
        arrow::as_arrow_table(combined[, ..id_cols])
    } else {
        # Non-tile parquetGeomStore: materialize whole. As the user
        # noted, callers who need streaming at this scale should be on
        # a tile store anyway.
        df <- storeRead(trim_store, output = "tibble", omit_internals = FALSE)
        out <- eval_chunk(df)
        if (is.null(out)) {
            return(arrow::as_arrow_table(df[integer(0L), id_cols, drop = FALSE]))
        }
        arrow::as_arrow_table(out)
    }
}
