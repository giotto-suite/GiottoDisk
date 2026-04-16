
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
    output = c("query", "tibble", "duckdb"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    checkmate::assert_character(fields, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble", "duckdb"))
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
    output = c("query", "tibble", "duckdb"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    checkmate::assert_character(fields, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble", "duckdb"))
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
    output = c("query", "tibble", "terra", "sf", "duckdb"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb"))
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
    switch(output,
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
    output = c("query", "tibble", "terra", "sf", "duckdb"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb"))
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
    switch(output,
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
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetGeomTileStore"), function(store,
    extent = NULL,
    tile_idx = NULL,
    fields = NULL,
    output = c("query", "tibble", "terra", "sf", "duckdb"),
    callback = NULL,
    duckdb_params = list(),
    omit_internals = TRUE,
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb"))
    fields <- .pstore_fields_requested(store, fields)
    lazy_fields <- .pstore_lazy_fields(store, fields, output)

    upstream_callback <- NULL
    if (!is.null(tile_idx)) {
        tile_idx <- as.integer(tile_idx)
        upstream_callback <- function(atab) {
            atab <- dplyr::filter(atab, tile_index %in% tile_idx)
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
    switch(output,
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

    # "terra": use WKB directly — avoids sf R-level geometry allocation overhead
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
    if (!output %in% c("query", "duckdb")) { # if a materialized format...
        lazy <- unique(c("source_id", "row_index", lazy))
    }
    if (output %in% c("terra", "sf")) {
        lazy <- unique(c("geom", lazy))
    }
    attr(lazy, "lazy") <- TRUE
    lazy
}
