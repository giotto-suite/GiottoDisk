

# definitions ####

setMethod("storeRead", signature("ANY"), function(store, ...) {
    stop(sprintf("Reading not implemented for store type %s\n", class(store)))
})

setMethod("storeRead", signature("fileStore"), function(store, ...) {
    if (.is_empty_fun(store@read_fun)) {
        stop("[storeRead] a specific 'read_fun' must be provided for `fileStore`\n",
             call. = FALSE)
    }
    if (!file.exists(store@path)) {
        stop("[storeRead] file does not exist\n", call. = FALSE)
    }
    store@read_fun(store@path)
})

setMethod("storeRead", signature("parquetStore"), function(store, fields = NULL, output = c("query", "tibble"), ...) {
    GiottoUtils::package_check("arrow")
    checkmate::assert_character(fields, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble"))
    atab <- callNextMethod(store = store, ...)
    res <- dplyr::arrange(atab, row_index)
    if (!is.null(fields)) {
        getcols <- unique(c("row_index", fields))
        atab <- dplyr::select(atab, arrow::all_of(getcols))
    }
    switch(output,
        "query" = atab,
        "tibble" = dplyr::collect(atab)[, fields]
    )
})

setMethod("storeRead", signature("parquetGeomStore"), function(store,
    extent = NULL, fields = NULL,
    output = c("query", "tibble", "terra", "sf"), ...) {
    GiottoUtils::package_check("arrow")
    checkmate::assert_character(fields, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf"))
    atab <- callNextMethod(store, ...) # parquetStore

    # extent filtering
    if (length(store@extent) == 0L) {
        stop("[storeRead] extent info not found\n", call. = FALSE)
    }
    e_final <- ext(store@extent)
    if (!is.null(extent)) {
        e_final <- terra::intersect(e, ext(extent))
    }
    if (is.null(e_final)) {
        stop("[storeRead] No geometries within requested extent\n",
             call. = FALSE)
    }

    atab <- .dplyr_crop(atab,
        sdimx = "x_index",
        sdimy = "y_index",
        extent = e_final,
        inclusive = TRUE
    )

    # fields filtering
    if (!is.null(fields)) {
        getcols <- unique(c("row_index", fields))
        if (output %in% c("terra", "sf")) getcols <- unique(c("geom", getcols))
        atab <- dplyr::select(atab, arrow::all_of(getcols))
    }
    switch(output,
        "query" = atab,
        "tibble" = {
            data <- dplyr::collect(atab)
            if (!is.null(fields)) data <- data[, fields]
            data
        },
        "sf" = .parquet_format_to_spatial(atab, output = "sf", fields = fields),
        "terra" = .parquet_format_to_spatial(atab, output = "terra", fields = fields)
    )
})

setMethod("storeRead", signature("parquetGeomTileStore"), function(store,
    extent = NULL, tile = NULL, fields = NULL,
    output = c("query", "tibble", "terra", "sf"), ...) {
    GiottoUtils::package_check("arrow")
    checkmate::assert_character(fields, null.ok = TRUE)
    output <- match.arg(output, choices = c("query", "tibble", "terra", "sf"))

    atab <- callNextMethod(store, ...)

    # extent filtering
    if (length(store@extent) == 0L) {
        stop("[storeRead] extent info not found\n", call. = FALSE)
    }
    e_final <- ext(store@extent)
    if (!is.null(extent)) {
        e_final <- terra::intersect(e, ext(extent))
    }
    if (is.null(e_final)) {
        stop("[storeRead] No geometries within requested extent\n",
             call. = FALSE)
    }

    if (!is.null(tile)) {
        tile <- as.integer(tile)
        atab <- dplyr::filter(atab, tile_index %in% tile)
    }
    atab <- .dplyr_crop(atab,
        sdimx = "x_index",
        sdimy = "y_index",
        extent = e_final,
        inclusive = TRUE
    )

    # fields filtering
    if (!is.null(fields)) {
        getcols <- unique(c("row_index", fields))
        if (output %in% c("terra", "sf")) getcols <- unique(c("geom", getcols))
        atab <- dplyr::select(atab, arrow::all_of(getcols))
    }

    # spatial class coercion
    if (output %in% c("terra", "sf")) {

    }
    switch(output,
        "query" = atab,
        "tibble" = {
            data <- dplyr::collect(atab)
            if (!is.null(fields)) data <- data[, fields]
            data
        },
        "sf" = .parquet_format_to_spatial(atab, output = "sf", fields = fields),
        "terra" = .parquet_format_to_spatial(atab, output = "terra", fields = fields)
    )
})

setMethod("storeRead", signature("h5ArrayStore"), function(store, ...) {
    HDF5Array::HDF5Array(
        filepath = store@path,
        name = store@name,
        ...
    )
})

setMethod("storeRead", signature("tileDBMatrixStore"), function(store, ...) {
   TileDBArray::TileDBArray(path = store@path, attr = store@name)
})


# internals ####

# should not be performed on the whole, only in chunks or tiles
.parquet_format_to_spatial <- function(data, output = c("terra", "sf"), fields = NULL) {
    output <- match.arg(output, choices = c("terra", "sf"))
    # sf readin
    if (!is.null(fields)) {
        fields <- unique(c(fields, "geom"))
        data <- data %>% dplyr::select(dplyr::any_of(fields))
    }
    sfdata <- dplyr::collect(data) %>%
        sf::st_as_sf(sf_column_name = "geom")
    switch(output,
        "sf" = sfdata,
        "terra" = terra::vect(sfdata)
    )
}
