
#' @name storeRead
#' @title Read a `dataStore`
#' @description
#' Read from a `dataStore` inheriting object. The output should be a useful
#' representation of the contained data.
#' @param store `dataStore` inheriting object
#' @param extent `SpatExtent` to filter on (optional)
#' @param tile_idx `integerlike` (optional) specific tile number(s) to read
#' @param fields `character` (optional) specific fields/columns to read
#' @param output `character` (default = "query"). Format to get values in:
#'
#'   * "query" - produces an arrow lazy query
#'   * "tibble" - materialized dplyr tibble
#'   * "terra" - materialized `SpatVector`
#'   * "sf" - materialized `sf` object
#'   * "duckdb" - (requires {duckdb} and {dbplyr}) produces a `tbl_dbi`
#'      lazy query. **Note:** should not be used in a parallelized
#'      context as duckdb handles parallelization internally.
#'   * "sedona" - (requires {sedonadb}) produces a `sedonadb_dataframe`.
#'      All pending Arrow-phase ops (`@ops`, `@crop`, `@window`) are
#'      applied via the Arrow pipeline before handing off. Pending affine
#'      transforms are applied via `ST_Affine` on the `geom` column.
#'      **Note:** `x_index`/`y_index` are not updated after a transform —
#'      use SedonaDB spatial functions for geometry-based filtering.
#' @param omit_internals `logical` (default `TRUE`). Whether to drop internal
#'   special columns (`row_index`, `source_id`) from materialized output
#'   (`"tibble"`, `"terra"`, `"sf"`). Set to `FALSE` to retain them, e.g.
#'   when downstream code needs these columns for joining.
#'   Has no effect for `"query"` or `"duckdb"` outputs (those always expose
#'   all columns).
#' @param ... additional params to pass (if any implemented)
NULL

# definitions ####

setMethod("storeRead", signature("ANY"), function(store, ...) {
    stop(sprintf("Reading not implemented for store type %s\n", class(store)))
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("fileStore"), function(store, ...) {
    if (.is_empty_fun(store@read_fun)) {
        stop("[storeRead] a specific 'read_fun' must be provided for `fileStore`\n",
             call. = FALSE)
    }
    
    .guard_store_written(store@path)
    .store_simple_read(store, ...)
})

#' @rdname storeRead
#' @param callback (optional) `function` where the first param
#'   should accept the {arrow} query. A function to apply to the
#'   query prior to output and after fields or other filters are
#'   applied.
#' 
#'   Mostly useful for outputs that require materialization.
#' @param duckdb_params named `list`. Params to pass
#'   to [duckdb::duckdb_register_arrow()] if `output = "duckdb"`.
#'   Key params:
#' 
#'   * `conn` - DBI connection to a duckdb instance.
#'   * `name` - `character` (optional) If not provided, a random ID
#'     for the registered table will be generated
#' @export
setMethod("storeRead", signature("queryableStore"), function(store,
    fields = NULL, 
    output = c("query", "tibble", "duckdb"), 
    callback = NULL,
    duckdb_params = list(),
    ...) {
    GiottoUtils::package_check("arrow")
    checkmate::assert_character(fields, null.ok = TRUE)
    checkmate::assert_function(callback, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble", "duckdb"))
    atab <- callNextMethod(store = store, ...)
    if (!is.null(fields)) {
        atab <- dplyr::select(atab, dplyr::all_of(fields))
    }
    if (!is.null(callback)) atab <- callback(atab)
    switch(output,
        "query" = atab,
        "tibble" = dplyr::collect(atab),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params)
    )
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetStore"), function(store,
    fields = NULL,
    output = c("query", "tibble", "duckdb", "sedona"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    checkmate::assert_character(fields, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble", "duckdb", "sedona"))
    if (output == "sedona") return(.pstore_to_sedona(store, fields = fields, ...))
    if (output == "duckdb") return(.pstore_to_duckdb(store,
        fields = fields, duckdb_params = duckdb_params, ...))
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output)

    # custom schema handling
    if (length(store@datatype) > 0L) {
        default_types <- as.list(.pstore_arrow_types(store))
        typelist <- lapply(default_types, .arrow_type_from_string)
        customlist <- lapply(store@datatype, .arrow_r_type_map)
        typelist <- modifyList(typelist, customlist)
        sc <- arrow::schema(typelist)

        formals(store@read_fun)$schema <- sc
    }

    atab <- callNextMethod(store, # queryableStore
        fields = lazy_fields,
        output = "query",
        callback = callback,
        ...
    )
  
    .pbase_storeread_processing(atab,
        store = store,
        fields = fields,
        output = output,
        duckdb_params = duckdb_params,
        omit_internals = omit_internals,
        ...
    )
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("unionParquetStore"), function(store,
    fields = NULL,
    output = c("query", "tibble", "duckdb", "sedona"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    checkmate::assert_character(fields, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble", "duckdb", "sedona"))
    if (output == "sedona") return(.pstore_to_sedona(store, fields = fields, ...))
    if (output == "duckdb") return(.pstore_to_duckdb(store,
        fields = fields, duckdb_params = duckdb_params, ...))
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output)

    # custom schema handling
    if (length(store@datatype) > 0L) {
        default_types <- as.list(store@params$arrow_types)
        typelist <- lapply(default_types, .arrow_type_from_string)
        customlist <- lapply(store@datatype, .arrow_r_type_map)
        typelist <- modifyList(typelist, customlist)
        sc <- arrow::schema(typelist)

        store@stores <- lapply(store@stores, function(substore) {
            formals(substore@read_fun)$schema <- sc
            substore
        })
    }

    atab <- arrow::open_dataset(lapply(store@stores, .store_simple_read))
    source_order <- vapply(store@stores, function(s) s@uid, FUN.VALUE = character(1L))

    if (!is.null(lazy_fields)) {
        atab <- dplyr::select(atab, dplyr::any_of(lazy_fields))
    }
    if (!is.null(callback)) atab <- callback(atab)

    .pbase_storeread_processing(atab,
        store = store,
        fields = fields,
        output = output,
        source_order = source_order,
        duckdb_params = duckdb_params,
        omit_internals = omit_internals,
        ...
    )
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("unionParquetGeomStore"), function(store,
    extent = NULL,
    fields = NULL,
    output = c("query", "tibble", "terra", "sf", "duckdb", "sedona"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb", "sedona"))
    if (output == "sedona") return(.pstore_to_sedona(store, fields = fields, extent = extent, ...))
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output)

    # extent checking
    e <- .pstore_active_extent(store) # either NULL or SpatExtent
    if (!is.null(extent)) {
        if (!is.null(e)) {
            e <- terra::intersect(e, ext(extent))
        } else {
            e <- ext(extent)
        }
    }

    upstream_callback <- NULL
    if (!is.null(e)) {
        upstream_callback <- function(atab) {
            .dplyr_crop(atab,
                sdimx = "x_index",
                sdimy = "y_index",
                extent = e,
                inclusive = TRUE
            )
        }
    }

    source_order <- vapply(store@stores, function(s) s@uid, FUN.VALUE = character(1L))

    atab <- callNextMethod(store, # unionParquetStore
        fields = lazy_fields,
        output = "query",
        callback = upstream_callback,
        ...
    )

    dropcols <- character(0L)
    if (!is.null(fields)) { # handle upstream special col injections
        dropcols <- setdiff(specialCols(store), fields)
    }

    if (!is.null(callback)) atab <- callback(atab)
    result <- switch(output,
        "query" = atab,
        "tibble" = .pstore_to_tibble(atab,
            dropcols = dropcols,
            arrangecols = c("source_id", "tile_index", "row_index"),
            source_order = source_order,
            omit_internals = omit_internals),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params),
        .pgstore_to_spatial(atab,
            output = output,
            dropcols = dropcols,
            arrangecols = c("source_id", "tile_index", "row_index"),
            source_order = source_order,
            crs = store@params$crs,
            use_xy_as_geom = isTRUE(store@params$use_xy_as_geom),
            omit_internals = omit_internals
        )
    )
    r_ops <- store@post_ops
    if (length(r_ops) > 0L && output %in% c("tibble", "terra", "sf")) {
        result <- .apply_post_ops(result, r_ops, output)
    }
    result
})

.pbase_storeread_processing <- function(atab, store,
    fields = NULL,
    output,
    source_order = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {

    # `spat_relate` ops are handled at this level (not via `.apply_op`) so we
    # can route the predicate through sedonadb and narrow the arrow query
    # with the surviving id rows, keeping the rest of the chain lazy on
    # the arrow side. The local `sr_cache` caches the id arrow Table per
    # op position so subsequent `spat_relate`s in the same chain don't
    # re-evaluate prior predicates; the cache lives only for this storeRead
    # invocation and is never written to `store@ops`.
    sr_cache <- list()
    for (i in seq_along(store@ops)) {
        op <- store@ops[[i]]
        if (identical(op$type, "spat_relate")) {
            ids_tab <- .spat_relate_narrow(store, i, sr_cache)
            sr_cache[[as.character(i)]] <- ids_tab
            id_cols <- names(ids_tab)
            atab <- dplyr::semi_join(atab, ids_tab, by = id_cols)
        } else {
            atab <- .ptabular_apply_op(atab, op)
        }
    }

    # Narrow back to caller-requested fields. The upstream projection may
    # have widened the schema with op-referenced cols (so filter exprs
    # could resolve); drop them now. Keep specialCols since materialization
    # paths (.pstore_to_tibble, .pgstore_to_spatial) need row_index for
    # ordering and other specials for geom rebuild; specials get cleaned
    # later via `dropcols` / `omit_internals`. `any_of` is defensive
    # against schema changes by distinct() or other ops.
    if (!is.null(fields)) {
        keep <- unique(c(fields, specialCols(store)))
        atab <- dplyr::select(atab, dplyr::any_of(keep))
    }

    dropcols <- character(0L)
    if (!is.null(fields)) {
        # remove undesired special inject cols
        dropcols <- c(dropcols, setdiff(specialCols(store), fields))
    }

    switch(output,
        "query" = atab,
        "tibble" = .pstore_to_tibble(atab,
            dropcols = dropcols,
            source_order = source_order,
            omit_internals = omit_internals),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params)
    )
}

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetGeomStore"), function(store,
    extent = NULL,
    fields = NULL,
    output = c("query", "tibble", "terra", "sf", "duckdb", "sedona"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb", "sedona"))
    if (output == "sedona") return(.pstore_to_sedona(store, fields = fields, extent = extent, ...))
    if (output == "duckdb") return(.pstore_to_duckdb(store,
        fields = fields, extent = extent,
        duckdb_params = duckdb_params, ...))
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output)

    # extent checking
    e <- .pstore_active_extent(store) # either NULL or SpatExtent
    if (!is.null(extent)) {
        if (!is.null(e)) {
            e <- terra::intersect(e, ext(extent))
        } else {
            e <- ext(extent)
        }
    }

    upstream_callback <- NULL
    if (!is.null(e)) {
        upstream_callback <- function(atab) {
            # extent filtering
            .dplyr_crop(atab,
                sdimx = "x_index",
                sdimy = "y_index",
                extent = e,
                inclusive = TRUE
            )
        }
    }

    atab <- callNextMethod(store, # parquetStore
        fields = lazy_fields,
        output = "query",
        callback = upstream_callback,
        ...
    )

    dropcols <- character(0L)
    if (!is.null(fields)) { # handle upstream special col injections
        dropcols <- setdiff(specialCols(store), fields)
    }

    if (!is.null(callback)) atab <- callback(atab)
    result <- switch(output,
        "query" = atab,
        "tibble" = .pstore_to_tibble(atab,
            dropcols = dropcols,
            omit_internals = omit_internals),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params),
        .pgstore_to_spatial(atab,
            output = output,
            dropcols = dropcols,
            crs = store@params$crs,
            use_xy_as_geom = isTRUE(store@params$use_xy_as_geom),
            omit_internals = omit_internals
        )
    )
    r_ops <- store@post_ops
    if (length(r_ops) > 0L && output %in% c("tibble", "terra", "sf")) {
        result <- .apply_post_ops(result, r_ops, output)
    }
    result
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetGeomTileStore"), function(store,
    extent = NULL,
    tile_idx = NULL,
    fields = NULL,
    output = c("query", "tibble", "terra", "sf", "duckdb", "sedona"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb", "sedona"))
    if (output == "sedona") return(.pstore_to_sedona(store, fields = fields, extent = extent, tile_idx = tile_idx, ...))
    if (output == "duckdb") return(.pstore_to_duckdb(store,
        fields = fields, extent = extent, tile_idx = tile_idx,
        duckdb_params = duckdb_params, ...))
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output)

    upstream_callback <- NULL
    effective_tile_idx <- if (!is.null(tile_idx)) {
        as.integer(tile_idx)
    } else if (length(store@tile_filter) > 0L) {
        store@tile_filter
    }
    if (!is.null(effective_tile_idx)) {
        upstream_callback <- function(atab) {
            atab <- dplyr::filter(atab, tile_index %in% effective_tile_idx)
        }
    }

    atab <- callNextMethod(store,
        extent = extent,
        fields = lazy_fields,
        output = "query",
        callback = upstream_callback,
        ...
    )

    dropcols <- c("source_id")
    if (!is.null(fields)) { # handle upstream special col injections
        dropcols <- setdiff(specialCols(store), fields)
    }

    if (!is.null(callback)) atab <- callback(atab)
    result <- switch(output,
        "query" = atab,
        "tibble" = .pstore_to_tibble(atab,
            dropcols = dropcols,
            arrangecols = c("source_id", "tile_index", "row_index"),
            omit_internals = omit_internals),
        "duckdb" = .arrow_to_duckdb(atab,
            duckdb_params = duckdb_params),
        .pgstore_to_spatial(atab,
            output = output,
            dropcols = dropcols,
            arrangecols = c("source_id", "tile_index", "row_index"),
            crs = store@params$crs,
            use_xy_as_geom = isTRUE(store@params$use_xy_as_geom),
            omit_internals = omit_internals
        )
    )
    r_ops <- store@post_ops
    if (length(r_ops) > 0L && output %in% c("tibble", "terra", "sf")) {
        result <- .apply_post_ops(result, r_ops, output)
    }
    result
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("h5ArrayStore"), function(store, ...) {
    HDF5Array::HDF5Array(
        filepath = store@path,
        name = store@params$name,
        ...
    )
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("tileDBArrayStore"), function(store, ...) {
    TileDBArray::TileDBArray(
        path = store@path,
        attr = store@params$name,
        ...
    )
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("bpcMatrixStore"), function(store, ...) {
    BPCells::open_matrix_dir(store@path, ...)
})


# as.data.frame ####

#' @export
setMethod("as.data.frame", "parquetBase", function(x, ...) {
    as.data.frame(storeRead(x, output = "tibble", ...))
})

#' @export
setMethod("as.data.frame", "parquetGeomBase",
    function(x, row.names = NULL, optional = FALSE, geom = NULL, ...) {
    if (is.null(geom)) {
        return(as.data.frame(storeRead(x, output = "tibble", ...)))
    }
    geom <- match.arg(geom, c("XY"))
    if (nzchar(x@geomtype) && x@geomtype != "points" &&
            !isTRUE(x@params$use_xy_as_geom)) {
        stop("geom = \"XY\" requires a point geometry store")
    }
    xy_cols <- c("x_index", "y_index")
    # When @fields is set by a prior select op, inject xy_cols directly so
    # .pstore_fields_requested() keeps them through its intersection.
    if (!is.null(x@fields)) {
        x@fields <- unique(c(x@fields, xy_cols))
        tbl <- storeRead(x, output = "tibble", ...)
    } else {
        tbl <- storeRead(x, output = "tibble",
            fields = c(colnames(x), xy_cols), ...)
    }
    as.data.frame(dplyr::rename(tbl, x = "x_index", y = "y_index"))
})

# as.vector ####

#' @export
setMethod("as.vector", "parquetBase", function(x, mode = "any") {
    cols <- colnames(x)
    if (length(cols) != 1L) stop(
        "as.vector() requires exactly one column selected; use x[, col] first",
        call. = FALSE
    )
    val <- storeRead(x, output = "query") |> dplyr::pull(cols)
    stats::setNames(list(val), cols)
})

# as.terra ####

#' @export
setMethod("as.terra", "parquetGeomBase", function(x, ...) {
    storeRead(x, output = "terra", ...)
})


# internals ####

.guard_store_written <- function(path) {
    bool <- file.exists(path) # may be more than one
    if (any(!bool)) {
        unwritten <- path[!bool]
        stop("[storeRead] file does not exist:\n  ", 
            paste(unwritten, collapse = "\n  "), call. = FALSE)
    }
}

# `atab` - lazy arrow query table
# `dropcols` - character vector of cols to drop after materialization (if present)
# `arrangecols` - character vector of cols to arrange on. Order matters.
.pstore_to_tibble <- function(atab,
    dropcols = character(0L),
    arrangecols = c("source_id", "row_index"),
    source_order = NULL,
    selection = integer(0L),
    omit_internals = TRUE) {
    # enforced drops
    if (isTRUE(omit_internals)) {
        dropcols <- unique(c("row_index", "source_id", dropcols))
    }
  
    if (!"row_index" %in% names(atab)) {
        warning("[storeRead][parquet->tibble] row_index missing\n",
            "  Materialized row order is indeterminate", call. = FALSE)
    } else {
        atab <- dplyr::arrange(atab, 
            dplyr::across(dplyr::any_of(arrangecols)))
    }
    data <- dplyr::collect(atab)
    
    # align materialized with expected ordering
    if (!is.null(source_order) && "source_id" %in% names(data)) {
        data$source_id <- factor(data$source_id, levels = source_order)
        data <- dplyr::arrange(data, dplyr::across(dplyr::any_of(arrangecols)))
    }
    # align materialized with user-specified ordering
    if (length(selection) > 0L && is.unsorted(selection, na.rm = TRUE)) {
        sel_order <- rank(selection, ties.method = "first")
        data <- data[sel_order, , drop = FALSE]
    }
  
    dplyr::select(data, -dplyr::any_of(dropcols))
}

.arrow_to_duckdb <- function(atab, 
    duckdb_params = list()) {
    checkmate::assert_list(duckdb_params)
    package_check("duckdb")
    package_check("dbplyr")
  
    a <- duckdb_params
    if (is.null(a$conn)) {
        stop("[storeRead][parquet->duckdb] param error:\n",
            "duckdb_params$conn: duckdb connection needed", 
            call. = FALSE)
    }
    a$arrow_scannable <- atab
    a$name <- a$name %||% .make_uid()
    do.call(duckdb::duckdb_register_arrow, a)
  
    dplyr::tbl(a$conn, a$name)
}

# should not be performed on the whole, only in chunks or tiles
.pgstore_to_spatial <- function(atab,
    output = c("terra", "sf"),
    dropcols = character(0L),
    crs = NULL,
    use_xy_as_geom = FALSE,
    omit_internals = TRUE,
    ...) { # passthrough to .pstore_to_tibble()
    output <- match.arg(output, choices = c("terra", "sf"))

    if (isTRUE(use_xy_as_geom)) {
        # centroid path: build point geometries from x_index/y_index
        dropcols <- unique(c("geom", dropcols))
        data <- .pstore_to_tibble(atab, dropcols = dropcols,
            omit_internals = omit_internals, ...)
        sv <- terra::vect(data, geom = c("x_index", "y_index"))
        data$x_index <- NULL
        data$y_index <- NULL
        if (!is.null(crs) && nzchar(crs)) terra::crs(sv) <- crs
        return(sv)
    }

    # enforced drops (never drop geom)
    dropcols <- setdiff(
        unique(c("x_index", "y_index", dropcols)),
        "geom"
    )

    if (!"geom" %in% names(atab)) {
        stop("[storeRead][parquet->spatial] geom col missing\n")
    }

    data <- .pstore_to_tibble(atab, dropcols = dropcols,
        omit_internals = omit_internals, ...)

    if (output == "sf") {
        sfdata <- sf::st_as_sf(data, sf_column_name = "geom")
        if (!is.null(crs) && nzchar(crs)) {
            sfdata <- sf::st_set_crs(sfdata, crs)
        }
        return(sfdata)
    }

    # "terra": use WKB directly -- avoids sf R-level geometry allocation overhead
    wkb <- as.list(data$geom)
    data$geom <- NULL
    sv <- terra::vect(wkb)
    if (!is.null(crs) && nzchar(crs)) {
        terra::crs(sv) <- crs
    }
    if (ncol(data) > 0L) {
        terra::values(sv) <- data
    }
    sv
}

# .pstore_to_sedona ####
#
# Translate a store's full state into a SedonaDB lazy dataframe backed by the
# original on-disk parquet files. SedonaDB reads the files directly via DataFusion,
# preserving GeoParquet metadata so geometry columns are natively recognised.
#
# Hive partition columns (`source_id`, `tile_index`) require special handling:
# the sedonadb R binding hard-codes `GeoParquetReadOptions::default()` and does
# not expose `table_partition_cols`, so DataFusion does not auto-promote hive
# directory segments to columns. Without them, `row_index` collides across tiles
# and source provenance is lost. We compensate at the R/SQL layer by registering
# each `tile_index=<n>/` directory as its own sub-view and reconstructing the
# partition cols via SQL literal injection. Per-tile SELECTs are UNION ALL'd
# into a single base view that downstream WHERE/projection/affine SQL targets.
#
# Pipeline:
#   1. Resolve tile specs (`.pstore_tile_specs`): one entry per
#      surviving `source_id=<uid>/tile_index=<n>/` directory. Honors
#      `tile_idx` arg and `@tile_filter` via file-level pruning — the pruned
#      tile dirs never reach DataFusion.
#   2. Register each tile dir as a sub-view; build per-tile
#      `SELECT *, '<uid>' AS source_id, <n> AS tile_index FROM <tile_view>`.
#      UNION ALL the per-tile SELECTs as the base view.
#   3. SQL WHERE from: @crop/@window AABB on x_index/y_index (parquet stats
#      pushdown), @ops filter exprs via .r_expr_to_sql().
#   4. SELECT DISTINCT / LIMIT for distinct/head ops; tail/sample/join warn
#      and skip.
#   5. If an affine transform is pending: wrap in outer SELECT with ST_Affine
#      (one LIMIT 0 schema probe to enumerate non-geom columns).
#
# x_index/y_index are NOT updated after a transform — use SedonaDB spatial
# predicates for geometry-based filtering post-transform.
.pstore_to_sedona <- function(store, fields = NULL, ...) {
    GiottoUtils::package_check("sedonadb",
        repository = "github:apache/sedona-db/r/sedonadb")

    extra_args <- list(...)
    extent_arg <- extra_args$extent
    tile_idx_arg <- extra_args$tile_idx

    # Resolve fields the same way the arrow/duckdb path does — combines the
    # caller's `fields` arg with `store@fields` (set by `[, j]`), then expands
    # to include geom-store specials required by the sedona output.
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output = "sedona")

    aff <- if (inherits(store, "parquetGeomBase")) .pgeom_pending_transform(store) else NULL

    # --- 1. Resolve tile specs + register per-tile sub-views ---
    specs <- .pstore_tile_specs(store, tile_idx_arg = tile_idx_arg)
    if (length(specs) == 0L) {
        stop("[storeRead][sedona] no tile directories match the requested filter\n",
            "  store path: ", toString(storePaths(store)), call. = FALSE)
    }
    # View names must be lowercase: DataFusion normalises unquoted identifiers
    # to lowercase at registration time, which would break double-quoted SQL refs.
    base_view_name <- tolower(paste0("gd_", .make_uid()))
    union_sqls <- vapply(seq_along(specs), function(i) {
        spec <- specs[[i]]
        tile_view_name <- sprintf("%s_t%d", base_view_name, i)
        sedonadb::sd_to_view(
            sedonadb::sd_read_parquet(spec$dir_path),
            tile_view_name, overwrite = TRUE
        )
        sprintf('SELECT *, %s FROM "%s"',
            .pstore_tile_literal_cols(spec), tile_view_name)
    }, FUN.VALUE = character(1L))
    base_sql <- paste(union_sqls, collapse = " UNION ALL ")
    sedonadb::sd_to_view(
        sedonadb::sd_sql(base_sql),
        base_view_name, overwrite = TRUE
    )

    # --- 2. Build inner SQL via shared helper ---
    # geom_sql_fn handles the SRID-aware ST_GeomFromText (DataFusion).
    sedona_geom_sql_fn <- function(op, store) {
        wkt_esc <- gsub("'", "''", op$y_wkt, fixed = TRUE)
        srid <- .spatrelate_store_srid(store)
        if (is.na(srid)) {
            sprintf("ST_GeomFromText('%s')", wkt_esc)
        } else {
            sprintf("ST_GeomFromText('%s', %d)", wkt_esc, srid)
        }
    }
    sedona_register_id_fn <- function(ids_tab) {
        name <- tolower(paste0("gd_idf_", .make_uid()))
        # `sd_to_view` accepts data.frame / nanoarrow, not arrow Table.
        sedonadb::sd_to_view(as.data.frame(ids_tab), name, overwrite = TRUE)
        name
    }
    inner_sql <- .pstore_sql_inner(store, base_view_name, extent_arg,
        lazy_fields,
        geom_sql_fn = sedona_geom_sql_fn,
        register_id_fn = sedona_register_id_fn,
        engine = "sedona")

    # --- 3. Wrap with ST_Affine if transform is pending ---
    sedona_probe_fn <- function(sql) {
        sedonadb::sd_collect(sedonadb::sd_sql(
            sprintf("SELECT * FROM (%s) AS _probe LIMIT 0", sql)
        ))
    }
    inner_sql <- .pstore_sql_affine_wrap(inner_sql, aff, sedona_probe_fn,
        engine = "sedona")

    sdf <- sedonadb::sd_sql(inner_sql)
    attr(sdf, "view_name") <- base_view_name
    sdf
}

#' Get the registered view name for a GiottoDisk SedonaDB dataframe
#'
#' Returns the DataFusion view name as a double-quoted SQL identifier, suitable
#' for embedding directly in [sedonadb::sd_sql()] queries. Only valid on
#' `sedonadb_dataframe` objects returned by `storeRead(output = "sedona")`.
#'
#' @param sdf A `sedonadb_dataframe` from `storeRead(output = "sedona")`.
#' @returns A length-1 character: the double-quoted view name, e.g. `'"gd_abc123"'`.
#' @export
sd_view_ref <- function(sdf) {
    nm <- attr(sdf, "view_name", exact = TRUE)
    if (is.null(nm)) stop(
        "no 'view_name' attribute found; ",
        "sdf must be produced by storeRead(output = \"sedona\")",
        call. = FALSE
    )
    sprintf('"%s"', nm)
}


# .pstore_to_duckdb ####
#
# Translate a store's full state into a DuckDB lazy `tbl_dbi` backed by the
# original on-disk parquet files via DuckDB's `read_parquet` scanner.
# Mirrors `.pstore_to_sedona`: per-tile UNION ALL with `source_id` /
# `tile_index` SQL literal injection, WHERE for @crop/@window AABB + @ops
# filter exprs, DISTINCT / LIMIT for distinct/head, spatial extension's
# ST_* for spat_relate, ST_Affine wrap for pending transforms.
#
# Connection lifecycle:
#   - If `duckdb_params$conn` is supplied, the user owns it; we just register
#     temp views on it and return a tbl_dbi pointing at them.
#   - Otherwise an ephemeral in-memory connection is created. dbplyr::tbl
#     holds a reference (via the tbl_dbi's $src$con), so the connection
#     stays alive as long as the returned tbl is referenced and is gc'd
#     with it.
#
# `id_filter` ops (injected by `.spat_relate_narrow` for multi-spat_relate
# chains) are translated to a correlated EXISTS subquery against a duckdb
# temp table registered from the cached id arrow Table via
# `duckdb_register_arrow`.
.pstore_to_duckdb <- function(store, fields = NULL, duckdb_params = list(),
    ...) {
    GiottoUtils::package_check("duckdb")
    GiottoUtils::package_check("dbplyr")
    GiottoUtils::package_check("DBI")
    checkmate::assert_list(duckdb_params)

    extra_args <- list(...)
    extent_arg <- extra_args$extent
    tile_idx_arg <- extra_args$tile_idx

    # Connection: user-supplied or ephemeral. tbl_dbi holds a reference to
    # conn, so ephemeral connections are kept alive as long as the returned
    # tbl is referenced.
    conn <- duckdb_params$conn
    if (is.null(conn)) {
        conn <- DBI::dbConnect(duckdb::duckdb())
    }

    # Spatial extension is required for any geom-base store (ST_*, ST_Affine).
    # `INSTALL` is idempotent + cached after first run; `LOAD` per-connection.
    if (inherits(store, "parquetGeomBase")) {
        DBI::dbExecute(conn, "INSTALL spatial; LOAD spatial;")
    }

    # SQL paths (sedona / duckdb-native) don't need op_referenced col widening
    # in lazy_fields — WHERE references FROM scope regardless of SELECT
    # projection. Pass "sedona" to the helper to share that branch.
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output = "sedona")

    aff <- if (inherits(store, "parquetGeomBase")) .pgeom_pending_transform(store) else NULL

    # --- 1. Resolve tile specs + register base view ---
    specs <- .pstore_tile_specs(store, tile_idx_arg = tile_idx_arg)
    if (length(specs) == 0L) {
        stop("[storeRead][duckdb] no tile directories match the requested filter\n",
            "  store path: ", toString(storePaths(store)), call. = FALSE)
    }
    base_view_name <- tolower(paste0("gd_dd_", .make_uid()))
    union_sqls <- vapply(seq_along(specs), function(i) {
        spec <- specs[[i]]
        # `read_parquet('dir/*.parquet')` with hive_partitioning disabled —
        # we inject source_id / tile_index as SQL literals (mirrors sedona)
        # and want to avoid duckdb auto-promoting them from the directory
        # layout (which would also coerce them to types that may not match).
        path_escaped <- gsub("'", "''", spec$dir_path, fixed = TRUE)
        sprintf("SELECT *, %s FROM read_parquet('%s*.parquet', hive_partitioning=false)",
            .pstore_tile_literal_cols(spec), path_escaped)
    }, FUN.VALUE = character(1L))
    base_sql <- paste(union_sqls, collapse = " UNION ALL ")
    DBI::dbExecute(conn,
        sprintf('CREATE OR REPLACE TEMP VIEW "%s" AS %s',
            base_view_name, base_sql))

    # --- 2. Build inner SQL via shared helper ---
    # DuckDB spatial's ST_GeomFromText takes (VARCHAR [, BOOLEAN]) — no
    # SRID arg, and no CRS enforcement (geometries are just GEOMETRY).
    duckdb_geom_sql_fn <- function(op, store) {
        wkt_esc <- gsub("'", "''", op$y_wkt, fixed = TRUE)
        sprintf("ST_GeomFromText('%s')", wkt_esc)
    }
    duckdb_register_id_fn <- function(ids_tab) {
        name <- tolower(paste0("gd_idf_", .make_uid()))
        duckdb::duckdb_register_arrow(conn, name, ids_tab)
        name
    }
    inner_sql <- .pstore_sql_inner(store, base_view_name, extent_arg,
        lazy_fields,
        geom_sql_fn = duckdb_geom_sql_fn,
        register_id_fn = duckdb_register_id_fn,
        engine = "duckdb")

    # --- 3. Wrap with ST_Affine if transform is pending ---
    duckdb_probe_fn <- function(sql) {
        DBI::dbGetQuery(conn,
            sprintf("SELECT * FROM (%s) AS _probe LIMIT 0", sql))
    }
    inner_sql <- .pstore_sql_affine_wrap(inner_sql, aff, duckdb_probe_fn,
        engine = "duckdb")

    # Register final view + return a lazy tbl. Final-view registration lets
    # the user reference it via plain `tbl(conn, name)` and gives dbplyr a
    # stable handle.
    final_view_name <- tolower(paste0("gd_dd_final_", .make_uid()))
    DBI::dbExecute(conn,
        sprintf('CREATE OR REPLACE TEMP VIEW "%s" AS %s',
            final_view_name, inner_sql))
    tbl <- dplyr::tbl(conn, final_view_name)
    attr(tbl, "view_name") <- base_view_name
    tbl
}


# .pestore_to_duckdb ####
#
# The expression-store twin of `.pstore_to_duckdb`, and deliberately a much
# smaller function.
#
# `.pstore_to_duckdb` has to compile `@ops` into SQL text because the tabular
# ops have no form both engines accept -- `spat_relate` needs ST_*, `id_filter`
# needs EXISTS, `crop`/`window` are AABB literals, and `subset()` captures
# arbitrary R expressions that `.r_expr_to_sql` has to walk. The cost of that
# is a second implementation of the op registry, which is why the SQL side
# drops `join` / `tail` / `sample` with a warning and flattens op ordering.
#
# A parquetExprStore has no such problem: its subset predicates and its op
# chain are already dplyr, and dplyr lowers to Acero and to dbplyr alike. So
# this function only swaps the CARRIER -- build a tbl_dbi over DuckDB's own
# `read_parquet` instead of an arrow Dataset, then run the SAME
# `.pe_apply_axis_pred()` / `.pe_apply_ops()` the arrow path runs. That is what
# makes `output = "query"` and `output = "duckdb"` return the same values by
# construction rather than by test.
#
# Connection lifecycle mirrors `.pstore_to_duckdb`: a user-supplied conn is
# theirs and is used as-is, otherwise an ephemeral one is created and kept
# alive by the returned tbl_dbi's `$src$con`.
.pestore_to_duckdb <- function(store, fields = NULL, callback = NULL,
    duckdb_params = list()) {
    GiottoUtils::package_check("duckdb")
    GiottoUtils::package_check("dbplyr")
    GiottoUtils::package_check("DBI")
    checkmate::assert_list(duckdb_params)

    # Ours goes through `.duckdb_connect()` so the memory-limit option
    # applies; a user's connection keeps whatever settings they gave it.
    conn <- duckdb_params$conn %||% .duckdb_connect()

    specs <- .pstore_tile_specs(store)
    .pestore_guard_specs(store, specs)

    # hive_partitioning=false + a literal `source_id`: the payload itself is
    # only (row_id, col_id, value), so source_id has to be injected, and
    # letting duckdb infer it from the directory name would also let it pick
    # the type. As a literal it is VARCHAR, matching what arrow's hive
    # partitioning hands back and what the op payloads join against.
    base_view_name <- tolower(paste0("gd_pe_", .make_uid()))
    arms <- vapply(specs, function(spec) {
        path_escaped <- gsub("'", "''", spec$dir_path, fixed = TRUE)
        sprintf("SELECT *, %s FROM read_parquet('%s*.parquet', hive_partitioning=false)",
            .pstore_tile_literal_cols(spec), path_escaped)
    }, FUN.VALUE = character(1L))
    DBI::dbExecute(conn, sprintf('CREATE OR REPLACE TEMP VIEW "%s" AS %s',
        base_view_name, paste(arms, collapse = " UNION ALL ")))

    x <- dplyr::tbl(conn, base_view_name)

    # --- the arrow pipeline, verbatim ---
    # A union carries its substores' subset state in one composite
    # source_id-aware expression; DuckDB pushes that disjunction through the
    # UNION ALL and constant-folds each arm's source_id away, leaving a bare
    # id range on each scan. No per-arm construction needed.
    if (inherits(store, "unionParquetExprStore")) {
        filt <- .union_substore_filter_expr(store@stores)
        if (!is.null(filt)) x <- dplyr::filter(x, !!filt)
    } else {
        x <- .pe_apply_axis_pred(x, .pe_axis_pred(store@cell_idx), "row_id")
        x <- .pe_apply_axis_pred(x, .pe_axis_pred(store@gene_idx), "col_id")
    }
    x <- .pe_apply_ops(x, store@ops)
    if (!is.null(fields)) x <- dplyr::select(x, dplyr::all_of(fields))
    if (!is.null(callback)) x <- callback(x)

    # Register the compiled query as a view so `dbplyr::remote_name()` on the
    # result is a handle a caller can hit with plain SQL, as it is for
    # `.pstore_to_duckdb`.
    final_view_name <- tolower(paste0("gd_pe_final_", .make_uid()))
    DBI::dbExecute(conn, sprintf('CREATE OR REPLACE TEMP VIEW "%s" AS %s',
        final_view_name, as.character(dbplyr::remote_query(x))))
    tbl <- dplyr::tbl(conn, final_view_name)
    attr(tbl, "view_name") <- base_view_name
    tbl
}

# `@uid` names the on-disk `source_id=` directory. A store minted fresh from a
# path gets a NEW uid, which points this scan at a directory that is not there
# -- the arrow path answers that with an empty result, so the failure is worth
# naming here rather than leaving as duckdb's "No files found".
.pestore_guard_specs <- function(store, specs) {
    if (length(specs) == 0L) {
        stop("[storeRead][pestore->duckdb] no source partitions resolved for ",
            toString(class(store)), call. = FALSE)
    }
    missing <- !vapply(specs, function(s) dir.exists(s$dir_path), logical(1L))
    if (any(missing)) {
        stop("[storeRead][pestore->duckdb] no parquet partition at:\n  ",
            paste(vapply(specs[missing], `[[`, character(1L), "dir_path"),
                collapse = "\n  "),
            "\n  `@uid` must match the on-disk `source_id=` directory.",
            call. = FALSE)
    }
    invisible(NULL)
}


# Resolve the per-tile directories that should back a SQL backend's base
# view. Used by both `.pstore_to_sedona` and `.pstore_to_duckdb`; the
# returned spec list is engine-neutral (just directory paths + partition
# metadata).
#
# Returns a list of specs, each:
#   list(uid, tile_index, dir_path, has_tile_index)
#
# - `uid`: source uid for the substore (becomes `source_id` literal).
# - `tile_index`: integer index parsed from `tile_index=<n>/` (or 0L for stores
#   without a tile_index partition layer; in that case `has_tile_index = FALSE`
#   and tile_index is unused).
# - `dir_path`: filesystem directory passed to the engine's parquet reader
#   (`sd_read_parquet` / `read_parquet`); always ends with `/` so the engine
#   lists files under it rather than treating it as a single file path.
# - `has_tile_index`: TRUE for `parquetGeomBase`-inheriting stores (flat geom +
#   tiled geom both write a `tile_index=<n>` partition); FALSE for flat tabular
#   `parquetStore` which only has the `source_id=<uid>` partition.
#
# Pruning rules:
#   - If `tile_idx_arg` is supplied: keep only matching `tile_index=<n>` dirs.
#   - Else if `store@tile_filter` is set (parquetGeomTileStore only): keep
#     those.
#   - Else: include all tile dirs.
#
# `unionParquetStore` substores are walked in order — one spec per substore
# tile directory.
#
# `parquetExprStore` shares the no-tile branch: it writes the same
# `source_id=<uid>/` partition and has no tile layer, so the spec it produces
# is the flat tabular one. Its union must be named explicitly here — unlike
# every other union in the package it does not extend `fileStore`, so falling
# through to `list(store)` would reach for a `@path` slot the class lacks.
.pstore_tile_specs <- function(store, tile_idx_arg = NULL) {
    substores <- if (inherits(store, c("unionParquetStore",
                                       "unionParquetExprStore"))) {
        store@stores
    } else {
        list(store)
    }
    has_tile_index <- inherits(store, "parquetGeomBase")

    # Only parquetGeomTileStore exposes @tile_filter — flat geom stores have a
    # single tile_index=0 and no filter.
    eff_tile_idx <- if (!is.null(tile_idx_arg)) {
        as.integer(tile_idx_arg)
    } else if (inherits(store, "parquetGeomTileStore") &&
               length(store@tile_filter) > 0L) {
        store@tile_filter
    } else {
        NULL
    }

    unlist(lapply(substores, function(s) {
        uid <- s@uid
        source_dir <- file.path(s@path, paste0("source_id=", uid))
        if (!has_tile_index) {
            return(list(list(
                uid = uid,
                tile_index = 0L,
                dir_path = paste0(source_dir, "/"),
                has_tile_index = FALSE
            )))
        }
        tile_dirs <- list.files(source_dir, pattern = "^tile_index=",
            full.names = TRUE)
        if (length(tile_dirs) == 0L) return(list())
        tile_idxs <- as.integer(sub("^tile_index=", "", basename(tile_dirs)))
        if (!is.null(eff_tile_idx)) {
            keep <- tile_idxs %in% eff_tile_idx
            tile_dirs <- tile_dirs[keep]
            tile_idxs <- tile_idxs[keep]
        }
        if (length(tile_dirs) == 0L) return(list())
        Map(function(n, p) list(
            uid = uid,
            tile_index = n,
            dir_path = paste0(p, "/"),
            has_tile_index = TRUE
        ), tile_idxs, tile_dirs)
    }), recursive = FALSE)
}

# Build the per-tile literal-column SQL fragment that injects
# `source_id` (and `tile_index` when applicable) into a per-tile SELECT.
# Shared by `.pstore_to_sedona` and `.pstore_to_duckdb`.
.pstore_tile_literal_cols <- function(spec) {
    uid_lit <- gsub("'", "''", spec$uid, fixed = TRUE)
    if (spec$has_tile_index) {
        # CAST AS INT keeps the literal at int32 so it matches the on-disk
        # hive partition column type. Without the cast, DataFusion infers
        # Int64 (default for integer literals); sd_collect then surfaces
        # that as R `numeric`, and arrow::as_arrow_table() promotes it to
        # float64 -- which then fails to semi_join against the int32
        # tile_index column from the parquet schema. `INT` is the
        # portable Int32 name across DataFusion / DuckDB.
        sprintf("'%s' AS source_id, CAST(%d AS INT) AS tile_index",
            uid_lit, spec$tile_index)
    } else {
        sprintf("'%s' AS source_id", uid_lit)
    }
}


# Build the inner SELECT SQL for a parquet-backed SQL backend (sedona,
# duckdb). Translates `@crop`/`@window` (AABB), `@ops` (filter / head /
# distinct / spat_relate / id_filter), and final projection
# (`lazy_fields` / DISTINCT / `*`) into a single SQL string against
# `base_view_name`.
#
# Engine-specific bits are passed as functions:
#   - `geom_sql_fn(op, store)` -> SQL for the WKT side of a spat_relate
#     predicate (sedona includes SRID; duckdb omits it).
#   - `register_id_fn(ids_tab)` -> registers a cached id arrow Table as an
#     engine-side view, returns the view name (sedona uses sd_to_view,
#     duckdb uses duckdb_register_arrow).
# `engine` is a short label embedded in warning messages.
.pstore_sql_inner <- function(store, base_view_name, extent_arg,
        lazy_fields, geom_sql_fn, register_id_fn, engine) {
    where_clauses <- character(0L)

    # Crop/window AABB — float range on centroid columns; parquet stats
    # pushdown. @crop/@window only exist on geom-base stores; flat
    # parquetStore skips.
    e <- if (inherits(store, "parquetGeomBase")) .pstore_active_extent(store) else NULL
    if (!is.null(extent_arg)) {
        e <- if (!is.null(e)) terra::intersect(e, ext(extent_arg)) else ext(extent_arg)
    }
    if (!is.null(e)) {
        ev <- .ext_to_num_vec(e)
        where_clauses <- c(where_clauses, sprintf(
            "x_index >= %.17g AND x_index <= %.17g AND y_index >= %.17g AND y_index <= %.17g",
            ev[[1L]], ev[[2L]], ev[[3L]], ev[[4L]]
        ))
    }

    distinct_cols <- NULL
    limit_n <- NULL
    for (op in store@ops) {
        switch(op$type,
            "filter"   = where_clauses <- c(where_clauses, .r_expr_to_sql(op$expr)),
            "head"     = limit_n <- op$n,
            "distinct" = distinct_cols <- op$cols,
            "spat_relate" = {
                if (!is.null(op$y_wkt)) {
                    sql_pred <- .sql_relation_fn(op$relation)
                    geom_sql <- geom_sql_fn(op, store)
                    where_clauses <- c(where_clauses, sprintf(
                        "%s(geom, %s)", sql_pred, geom_sql
                    ))
                } else {
                    warning(sprintf(
                        "[storeRead][%s] spat_relate with stored y is not yet implemented in the %s compile; the op is skipped",
                        engine, engine
                    ), call. = FALSE)
                }
            },
            "id_filter" = {
                # Internal op type — injected by `.spat_relate_narrow` when
                # a chain has multiple spat_relate ops so we don't re-run
                # earlier predicates. Engine registers the cached ids as a
                # view, returns its name; this helper builds the correlated
                # EXISTS subquery on the join keys.
                id_view_name <- register_id_fn(op$ids_tab)
                keys <- unname(op$by)
                join_conds <- vapply(keys, function(k) {
                    sprintf('"%s"."%s" = "%s"."%s"',
                        id_view_name, k, base_view_name, k)
                }, FUN.VALUE = character(1L))
                where_clauses <- c(where_clauses, sprintf(
                    'EXISTS (SELECT 1 FROM "%s" WHERE %s)',
                    id_view_name, paste(join_conds, collapse = " AND ")
                ))
            },
            "tail"     = ,
            "sample"   = ,
            "join"     = warning(sprintf(
                "[storeRead][%s] op '%s' cannot be expressed as SQL and is skipped",
                engine, op$type), call. = FALSE),
            warning(sprintf("[storeRead][%s] unknown op type '%s' skipped",
                engine, op$type), call. = FALSE)
        )
    }

    where_sql <- if (length(where_clauses) > 0L) {
        paste("WHERE", paste(where_clauses, collapse = " AND "))
    } else {
        ""
    }
    # Column projection — DISTINCT op wins (it specifies its own cols);
    # otherwise honor `lazy_fields` if the caller narrowed via `[, j]` or
    # `fields = ...`; otherwise project everything.
    select_sql <- if (!is.null(distinct_cols)) {
        paste("DISTINCT", paste(sprintf('"%s"', distinct_cols), collapse = ", "))
    } else if (!is.null(lazy_fields)) {
        paste(sprintf('"%s"', lazy_fields), collapse = ", ")
    } else {
        "*"
    }
    inner_sql <- trimws(sprintf('SELECT %s FROM "%s" %s',
        select_sql, base_view_name, where_sql))
    if (!is.null(limit_n)) inner_sql <- paste(inner_sql, "LIMIT", limit_n)
    inner_sql
}


# Wrap an inner SELECT in an outer SELECT that applies a pending affine
# transform via ST_Affine. No-op when `aff` is NULL. The engine-specific
# `schema_probe_fn(inner_sql)` returns a data.frame whose names enumerate
# the projection's columns (used to emit non-geom cols verbatim and
# replace geom with ST_Affine(geom, ...)).
#
# `engine` selects the ST_Affine argument convention:
#   * "duckdb" (PostGIS) -- x' = a*x + b*y + xoff, y' = d*x + e*y + yoff
#   * "sedona"           -- x' = a*x + d*y + xoff, y' = b*x + e*y + yoff
# SedonaDB's ST_Affine uses the transposed convention relative to PostGIS;
# verified empirically (sd_sql("ST_Affine(ST_Point(3622,-2142),
# 0.866,0.5,-0.5,0.866,0,0)") returns (4207.65, -43.97) rather than
# PostGIS's (2065.65, -3665.97)). Without this engine-aware swap, any
# pending @post_ops rotation on a backed polygon store via the sedona
# engine silently emits the wrong transform -- the geom never rotates,
# ST_Intersects against the rotated-frame query returns nothing.
.pstore_sql_affine_wrap <- function(inner_sql, aff, schema_probe_fn,
    engine = c("duckdb", "sedona")) {
    if (is.null(aff)) return(inner_sql)
    engine <- match.arg(engine)
    schema_df <- schema_probe_fn(inner_sql)
    other_cols <- setdiff(names(schema_df), "geom")
    m <- aff@affine
    # Post-multiply convention: [x, y, 1] %*% M
    #   x' = x*M[1,1] + y*M[2,1] + M[3,1],  y' = x*M[1,2] + y*M[2,2] + M[3,2]
    coefs <- if (engine == "sedona") {
        # transposed arg order: ST_Affine(geom, a, b, d, e, ...) with
        # x' = a*x + d*y, so we feed (M[1,1], M[1,2], M[2,1], M[2,2]).
        c(m[1L, 1L], m[1L, 2L], m[2L, 1L], m[2L, 2L], m[3L, 1L], m[3L, 2L])
    } else {
        # PostGIS / duckdb: ST_Affine(geom, a, b, d, e, ...) with
        # x' = a*x + b*y, so (M[1,1], M[2,1], M[1,2], M[2,2]).
        c(m[1L, 1L], m[2L, 1L], m[1L, 2L], m[2L, 2L], m[3L, 1L], m[3L, 2L])
    }
    affine_expr <- sprintf(
        "ST_Affine(geom, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g) AS geom",
        coefs[1L], coefs[2L], coefs[3L], coefs[4L], coefs[5L], coefs[6L]
    )
    outer_select <- if (length(other_cols) > 0L) {
        paste(c(paste(sprintf('"%s"', other_cols), collapse = ", "), affine_expr),
            collapse = ", ")
    } else {
        affine_expr
    }
    sprintf("SELECT %s FROM (%s) AS _t", outer_select, inner_sql)
}


# Translate an inlined R call expression to a SQL WHERE fragment.
# Handles operators produced by subset() after .inline_local_vars():
#   ==, !=, <, <=, >, >=, &/&&, |/||, !, (, %in%, is.na
# Half-plane filter exprs (a*x_index + b*y_index >= c) translate naturally
# since arithmetic operators map directly to SQL.
# Falls back to deparse() with a warning for unrecognised operators.
.r_expr_to_sql <- function(expr) {
    if (is.name(expr)) return(sprintf('"%s"', as.character(expr)))
    if (is.atomic(expr)) {
        if (is.character(expr)) {
            escaped <- gsub("'", "''", expr)
            if (length(escaped) == 1L) return(sprintf("'%s'", escaped))
            return(sprintf("('%s')", paste(escaped, collapse = "', '")))
        }
        if (is.logical(expr)) return(if (isTRUE(expr)) "TRUE" else "FALSE")
        if (length(expr) == 1L) return(as.character(expr))
        return(sprintf("(%s)", paste(expr, collapse = ", ")))
    }
    if (!is.call(expr)) {
        warning("[storeRead][sedona] cannot translate expression to SQL; using deparse()",
            call. = FALSE)
        return(deparse(expr))
    }
    fn   <- as.character(expr[[1L]])
    args <- as.list(expr[-1L])
    L <- function() .r_expr_to_sql(args[[1L]])
    R <- function() .r_expr_to_sql(args[[2L]])
    switch(fn,
        "=="   = sprintf("%s = %s",     L(), R()),
        "!="   = sprintf("%s != %s",    L(), R()),
        "<"    = sprintf("%s < %s",     L(), R()),
        "<="   = sprintf("%s <= %s",    L(), R()),
        ">"    = sprintf("%s > %s",     L(), R()),
        ">="   = sprintf("%s >= %s",    L(), R()),
        "+"    = sprintf("%s + %s",     L(), R()),
        "-"    = if (length(args) == 1L) sprintf("-%s", L())
                 else sprintf("%s - %s", L(), R()),
        "*"    = sprintf("%s * %s",     L(), R()),
        "/"    = sprintf("%s / %s",     L(), R()),
        "&"    = ,
        "&&"   = sprintf("(%s AND %s)", L(), R()),
        "|"    = ,
        "||"   = sprintf("(%s OR %s)",  L(), R()),
        "!"    = sprintf("NOT (%s)",    .r_expr_to_sql(args[[1L]])),
        "("    = sprintf("(%s)",        .r_expr_to_sql(args[[1L]])),
        "%in%" = {
            rhs <- args[[2L]]
            # inline c(...) literals arrive as call objects — evaluate them
            if (is.call(rhs)) rhs <- eval(rhs)
            .sql_in_clause(L(), rhs)
        },
        "is.na" = sprintf("%s IS NULL", .r_expr_to_sql(args[[1L]])),
        {
            warning(sprintf(
                "[storeRead][sedona] unsupported op '%s' in filter expression; using deparse()",
                fn), call. = FALSE)
            deparse(expr)
        }
    )
}

.pstore_fields_requested <- function(store, fields = NULL) {
    if (isTRUE(attr(fields, "lazy"))) return(fields) # already expanded, pass through
    if (is.null(fields) && is.null(store@fields)) return(NULL)
    if (!is.null(fields) && !is.null(store@fields)) {
        return(intersect(fields, store@fields))
    }
    fields %||% store@fields
}

.pstore_lazy_fields <- function(store, fields, output) {
    if (is.null(fields)) return(NULL)
    if (isTRUE(attr(fields, "lazy"))) return(fields)

    # Restrict user fields to upstream-available (on-disk) cols. Any
    # requested fields that come from queued join ops' y-stores will be
    # materialized post-join by `.apply_op` and narrowed back via the final
    # select in `.pbase_storeread_processing` -- they must NOT appear in
    # the upstream parquet projection.
    disk <- .pstore_disk_fields(store) %||% character(0L)
    lazy <- intersect(fields, disk)

    # Widen with cols referenced by queued ops -- arrow/dplyr's `select`
    # strips cols from the lazy schema, so filter predicates referencing
    # `[, j]`-dropped cols would fail at compile time. SQL engines (sedona)
    # reference cols via FROM regardless of SELECT projection, so no
    # widening is needed there. Final narrowing back to `fields` happens
    # in `.pbase_storeread_processing` after ops apply.
    if (output != "sedona") {
        op_refs <- .pstore_op_referenced_cols(store)
        if (length(op_refs) > 0L) {
            lazy <- unique(c(op_refs, lazy))
        }
    }

    # store requirements
    if (inherits(store, "parquetGeomBase")) {
        lazy <- unique(c("x_index", "y_index", "tile_index", lazy))
    }

    # output requirements
    if (!output %in% c("query", "duckdb", "sedona")) { # if a materialized format...
        lazy <- unique(c("source_id", "row_index", lazy))
    }
    if (output %in% c("terra", "sf")) {
        lazy <- unique(c("geom", lazy))
    }
    attr(lazy, "lazy") <- TRUE
    lazy
}

# Cols referenced by queued ops that downstream filter compile needs visible
# in the lazy schema. Intersected with disk_fields so synthetic / undeclared
# symbols (already inlined by `.inline_local_vars`) don't leak through, and
# so that y-side cols referenced by post-join filters don't get projected
# from the upstream (x-only) parquet -- those come in naturally via the
# join op's .apply_op handler.
.pstore_op_referenced_cols <- function(store) {
    if (length(store@ops) == 0L) return(character(0L))
    refs <- lapply(store@ops, function(op) {
        switch(op$type %||% "",
            "filter"   = all.vars(op$expr),
            "distinct" = op$cols,
            # Join keys must be visible in the upstream projection so
            # `.apply_op`'s dplyr::inner_join can find them on the x side.
            # `op$by` follows data.table syntax: names() = x cols, values
            # = y cols (or unnamed for same-name joins).
            "join"     = {
                nm <- names(op$by)
                if (is.null(nm)) op$by else nm
            },
            # spat_relate needs the `geom` column materialized so the arrow
            # path can build a SpatVector and call terra::relate.
            "spat_relate" = "geom",
            character(0L)
        )
    })
    refs <- unique(unlist(refs))
    disk <- .pstore_disk_fields(store)
    if (is.null(disk)) return(refs)
    intersect(refs, disk)
}

# Effective schema of a store at the point downstream code is about to add
# the next op. Equals on-disk cols plus any cols recursively brought in by
# already-queued join ops. y's specials (row_index, source_id, etc.) are
# dropped except for the join keys, mirroring `.apply_op`'s join compile.
#
# Used to decide what symbols in `subset()` predicates and `[, j]` selectors
# are legitimate column references vs local R variables to inline. Computed
# JIT from `@ops` rather than cached on the store, so direct `@ops`
# manipulation (a supported pattern per AGENTS.md) stays correct.
.pstore_effective_schema <- function(store) {
    base <- .pstore_disk_fields(store) %||% character(0L)
    for (op in store@ops) {
        if (identical(op$type, "join")) {
            y_eff <- .pstore_effective_schema(op$y)
            y_drop <- setdiff(specialCols(op$y), unname(op$by))
            base <- unique(c(base, setdiff(y_eff, y_drop)))
        }
    }
    base
}

# Emit a SQL IN clause, switching to a VALUES subquery above the threshold to
# allow DataFusion to use a hash join rather than a flat literal list scan.
.sql_in_clause <- function(col_sql, vals) {
    threshold <- getOption("giottodisk.sedona_in_subquery_threshold", 1000L)
    if (length(vals) > threshold) {
        vals_sql <- if (is.character(vals)) {
            paste(sprintf("('%s')", gsub("'", "''", vals)), collapse = ", ")
        } else {
            paste(sprintf("(%s)", vals), collapse = ", ")
        }
        sprintf("%s IN (SELECT column1 FROM (VALUES %s) AS _in)", col_sql, vals_sql)
    } else {
        vals_str <- if (is.character(vals)) {
            paste(sprintf("'%s'", gsub("'", "''", vals)), collapse = ", ")
        } else {
            paste(vals, collapse = ", ")
        }
        sprintf("%s IN (%s)", col_sql, vals_str)
    }
}
