
#' @name storeRead
#' @title Read a `dataStore`
#' @description
#' Read from a `dataStore` inheriting object. The output should be a useful
#' representation of the contained data.
#' @param store `dataStore` inheriting object
#' @param extent `SpatExtent` to filter on (optional)
#' @param tile `integerlike` (optional) specific tile number(s) to read
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
    store@read_fun(store@path)
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
        atab <- dplyr::select(atab, arrow::all_of(fields))
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
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "duckdb"))
  
    lazy_fields <- fields
    if (output == "tibble" && !is.null(fields)) {
        # special inject row_index col for row ordering
        lazy_fields <- unique(c("source_id", "row_index", lazy_fields))
    }
    atab <- callNextMethod(store, # queryableStore
        fields = lazy_fields,
        output = "query",
        callback = callback,
        ...
    )
    
    dropcols <- character(0L)
    if (!is.null(fields)) {
        # remove undesired special inject cols
        dropcols <- c(dropcols, setdiff(specialCols(store), fields))
    }
  
    switch(output,
        "query" = atab,
        "tibble" = .pstore_to_tibble(atab, dropcols = dropcols),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params)
    )
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetGeomStore"), function(store,
    extent = NULL,
    fields = NULL,
    output = c("query", "tibble", "terra", "sf", "duckdb"), 
    callback = NULL,
    duckdb_params = list(),
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb"))
    
    # extent checking
    if (length(store@extent) == 0L) {
        stop("[storeRead] extent info not found\n", call. = FALSE)
    }
    e <- ext(store@extent)
    if (!is.null(extent)) {
        e <- terra::intersect(e, ext(extent))
    }
    if (is.null(e)) {
        stop("[storeRead] No geometries within requested extent\n",
             call. = FALSE)
    }

    # inject special cols if needed
    lazy_fields <- fields
    if (!is.null(fields)) {
        # needed for extent filtering in upstream_callback
        lazy_fields <- unique(c("x_index", "y_index", lazy_fields))
        if (!output %in% c("query", "duckdb")) { # row ordering
            lazy_fields <- unique(c("source_id", "row_index", lazy_fields))
        }
        if (output %in% c("terra", "sf")) {
            # needed for geometries
            lazy_fields <- unique(c("geom", lazy_fields))
        }
    }
  
    upstream_callback <- function(atab) {
        # extent filtering
        atab <- .dplyr_crop(atab,
            sdimx = "x_index",
            sdimy = "y_index",
            extent = e,
            inclusive = TRUE
        )
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
        "tibble" = .pstore_to_tibble(atab, dropcols = dropcols),
        "duckdb" = .arrow_to_duckdb(atab, duckdb_params = duckdb_params),
        .pgstore_to_spatial(atab, 
            output = output, 
            dropcols = dropcols,
            crs = store@params$crs
        )
    )
})

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("parquetGeomTileStore"), function(store,
    extent = NULL, 
    tile = NULL,
    fields = NULL,
    output = c("query", "tibble", "terra", "sf", "duckdb"),
    callback = NULL,
    duckdb_params = list(),
    ...) {
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf", "duckdb"))
  
    lazy_fields <- fields
    upstream_callback <- NULL
    if (!is.null(fields)) { # inject special cols if needed
        lazy_fields <- unique(c("source_id", "tile_index", lazy_fields))
    }
    if (!is.null(tile)) {
        # this callback is conditional
        tile <- as.integer(tile)
        upstream_callback <- function(atab) {
            # tile filtering
            atab <- dplyr::filter(atab, tile_index %in% tile)
        }
    }

    atab <- callNextMethod(store,
        extent = extent,
        fields = lazy_fields,
        output = "query",
        callback = upstream_callback,
        ...
    )
  
    dropcols <- c("tile_index", "source_id")
    if (!is.null(fields)) { # handle upstream special col injections
        dropcols <- setdiff(specialCols(store), fields)
    }

    if (!is.null(callback)) atab <- callback(atab)
    switch(output,
        "query" = atab,
        "tibble" = .pstore_to_tibble(atab, 
            dropcols = dropcols, 
            arrangecols = c("source_id", "tile_index", "row_index")),
        "duckdb" = .arrow_to_duckdb(atab, 
            duckdb_params = duckdb_params),
        .pgstore_to_spatial(atab, 
            output = output,
            dropcols = dropcols,
            arrangecols = c("source_id", "tile_index", "row_index"),
            crs = store@params$crs
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


# internals ####

.test_store_written <- function(path) {
    all(file.exists(path))
}

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
    arrangecols = c("source_id", "row_index")) {
    # enforced drops
    dropcols <- unique(c("row_index", "source_id", dropcols))
  
    if (!"row_index" %in% names(atab)) {
        warning("[storeRead][parquet->tibble] row_index missing\n",
            "  Materialized row order is indeterminate", call. = FALSE)
    } else {
        atab <- dplyr::arrange(atab, 
            dplyr::across(dplyr::any_of(arrangecols)))
    }
    atab |>
        dplyr::collect() |>
        dplyr::select(-dplyr::any_of(dropcols))
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
    arrangecols = c("source_id", "row_index"),
    crs = NULL) {
    output <- match.arg(output, choices = c("terra", "sf"))
    # enforced drops (never drop geom)
    dropcols <- setdiff(
        unique(c("x_index", "y_index", dropcols)), 
        "geom"
    )

    if (!"geom" %in% names(atab)) {
        stop("[storeRead][parquet->spatial] geom col missing\n")
    }
  
    # collect + consume indices -> sf readin
    sfdata <- .pstore_to_tibble(atab, 
        dropcols = dropcols,
        arrangecols = arrangecols) |>
        sf::st_as_sf(sf_column_name = "geom"
    )
    if (!is.null(crs) && nzchar(crs)) {
        sfdata <- sf::st_set_crs(sfdata, crs)
    }
    switch(output,
        "sf" = sfdata,
        "terra" = terra::vect(sfdata)
    )
}
