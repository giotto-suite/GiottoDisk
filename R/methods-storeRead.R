
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
    if (output == "sedona") return(.pstore_to_sedona(store, ...))
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
    if (output == "sedona") return(.pstore_to_sedona(store, ...))
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
    if (output == "sedona") return(.pstore_to_sedona(store, extent = extent, ...))
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

    for (op in store@ops) {
        atab <- .do_op(atab, op)
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
    if (output == "sedona") return(.pstore_to_sedona(store, extent = extent, ...))
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
    if (output == "sedona") return(.pstore_to_sedona(store, extent = extent, tile_idx = tile_idx, ...))
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
# Pipeline:
#   1. sd_read_parquet(paths) + sd_to_view() — DataFusion reads files directly;
#      GeoParquet metadata intact, no WKB conversion needed
#   2. SQL WHERE from: @crop/@window AABB on x_index/y_index (parquet stats pushdown),
#      @tile_filter as tile_index IN (...), @ops filter exprs via .r_expr_to_sql()
#   3. SELECT DISTINCT / LIMIT for distinct/head ops; tail/sample/join warn and skip
#   4. If an affine transform is pending: wrap in outer SELECT with ST_Affine
#      (one LIMIT 0 schema probe to enumerate non-geom columns)
#
# x_index/y_index are NOT updated after a transform — use SedonaDB spatial
# predicates for geometry-based filtering post-transform.
.pstore_to_sedona <- function(store, ...) {
    GiottoUtils::package_check("sedonadb",
        repository = "github:apache/sedona-db/r/sedonadb")

    extra_args <- list(...)
    extent_arg <- extra_args$extent
    tile_idx_arg <- extra_args$tile_idx

    aff <- if (inherits(store, "parquetGeomBase")) .pgeom_pending_transform(store) else NULL

    # --- 1. Register base parquet paths ---
    # View name must be all lowercase: DataFusion normalises unquoted identifiers
    # to lowercase at registration time, which would break double-quoted SQL refs.
    view_name <- tolower(paste0("gd_", .make_uid()))
    sedonadb::sd_to_view(
        sedonadb::sd_read_parquet(.pstore_sedona_paths(store)),
        view_name, overwrite = TRUE
    )

    # --- 2. Build WHERE clause ---
    where_clauses <- character(0L)

    # Crop/window AABB — float range on centroid columns; parquet stats pushdown
    e <- .pstore_active_extent(store)
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

    # Tile partition filter
    eff_tile_idx <- if (!is.null(tile_idx_arg)) {
        as.integer(tile_idx_arg)
    } else if (inherits(store, "parquetGeomTileStore") && length(store@tile_filter) > 0L) {
        store@tile_filter
    } else {
        NULL
    }
    if (!is.null(eff_tile_idx)) {
        where_clauses <- c(where_clauses,
            sprintf("tile_index IN (%s)", paste(eff_tile_idx, collapse = ", ")))
    }

    # @ops — filter (including half-plane exprs), head, distinct translated to SQL;
    # tail/sample/join have no clean SQL equivalent and are skipped with a warning
    distinct_cols <- NULL
    limit_n <- NULL
    for (op in store@ops) {
        switch(op$type,
            "filter"   = where_clauses <- c(where_clauses, .r_expr_to_sql(op$expr)),
            "head"     = limit_n <- op$n,
            "distinct" = distinct_cols <- op$cols,
            "tail"     = ,
            "sample"   = ,
            "join"     = warning(sprintf(
                "[storeRead][sedona] op '%s' cannot be expressed as SQL and is skipped",
                op$type), call. = FALSE),
            warning(sprintf("[storeRead][sedona] unknown op type '%s' skipped", op$type),
                call. = FALSE)
        )
    }

    # --- 3. Build inner SQL ---
    where_sql <- if (length(where_clauses) > 0L) {
        paste("WHERE", paste(where_clauses, collapse = " AND "))
    } else {
        ""
    }
    select_sql <- if (!is.null(distinct_cols)) {
        paste("DISTINCT", paste(sprintf('"%s"', distinct_cols), collapse = ", "))
    } else {
        "*"
    }
    inner_sql <- trimws(sprintf('SELECT %s FROM "%s" %s', select_sql, view_name, where_sql))
    if (!is.null(limit_n)) inner_sql <- paste(inner_sql, "LIMIT", limit_n)

    # --- 4. Wrap with ST_Affine if transform is pending ---
    # Post-multiply convention: [x, y, 1] %*% M
    #   x' = x*M[1,1] + y*M[2,1] + M[3,1],  y' = x*M[1,2] + y*M[2,2] + M[3,2]
    # ST_Affine(geom, a, b, d, e, xoff, yoff):
    #   x' = a*x + b*y + xoff,  y' = d*x + e*y + yoff
    if (!is.null(aff)) {
        schema_df <- sedonadb::sd_collect(
            sedonadb::sd_sql(sprintf("SELECT * FROM (%s) AS _probe LIMIT 0", inner_sql))
        )
        other_cols <- setdiff(names(schema_df), "geom")
        m <- aff@affine
        affine_expr <- sprintf(
            "ST_Affine(geom, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g) AS geom",
            m[1L, 1L], m[2L, 1L], m[1L, 2L], m[2L, 2L], m[3L, 1L], m[3L, 2L]
        )
        outer_select <- if (length(other_cols) > 0L) {
            paste(c(paste(sprintf('"%s"', other_cols), collapse = ", "), affine_expr),
                collapse = ", ")
        } else {
            affine_expr
        }
        inner_sql <- sprintf("SELECT %s FROM (%s) AS _t", outer_select, inner_sql)
    }

    sdf <- sedonadb::sd_sql(inner_sql)
    attr(sdf, "view_name") <- view_name
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

# Collect all parquet file paths for a store for sd_read_parquet().
# For directory-based stores (hive-partitioned tile stores), recursively
# finds all .parquet files under each base path.
.pstore_sedona_paths <- function(store) {
    base_paths <- storePaths(store)
    unlist(lapply(base_paths, function(p) {
        if (dir.exists(p)) {
            list.files(p, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
        } else {
            p
        }
    }), use.names = FALSE)
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
            vals <- if (is.character(rhs)) {
                paste(sprintf("'%s'", gsub("'", "''", rhs)), collapse = ", ")
            } else {
                paste(rhs, collapse = ", ")
            }
            sprintf("%s IN (%s)", L(), vals)
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
    lazy <- fields

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
