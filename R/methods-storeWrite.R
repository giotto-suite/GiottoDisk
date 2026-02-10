
# docs ####

#' @name storeWrite
#' @title Write to a `dataStore`
#' @description
#' Write data to a `dataStore` inheriting class (usually a subclass of
#' [fileStore-class]). This documentation covers the basic expected API.
#' For more functional examples, see documentation for specific store
#' implementations.
#'
#' * [storeWrite-parquetStore]
#' * [storeWrite-parquetGeomStore]
#' * [storeWrite-parquetGeomTileStore]
#'
#' @param store `dataStore` inheriting class
#' @param data data to write
#' @param ... additional params to pass
#' @family storeWrite methods
NULL




.guard_lazy_access <- function(x) {
    if (!inherits(x, c("FileSystemDataset", "tbl_lazy", "arrow_dplyr_query"))) {
        stop("[storeWrite] 'data' should be a `FileSystemDataset`, `arrow_dplyr_query`, or `tbl_lazy`")
    }
}

.guard_no_data_rows <- function(x) {
    if (nrow(x) == 0L) stop("[storeWrite] nrow = 0. No values to write.")
}



# * general ####

setMethod("storeWrite", signature("ANY", "ANY"), function(store, data, ...) {
    stop(sprintf("Writing not implemented for store type %s\n", class(store)),
         call. = FALSE)
})

# from fileStore ETL chaining
# e.g. data coercible to arrow FileSystemDataset -> <fileStore> -> read to FSD -> <parquetStore> write
#' @rdname storeWrite
#' @export
setMethod("storeWrite", signature("fileStore", "fileStore"), function(store, data, ...) {
    storeWrite(store, storeRead(data))
})

# * parquetStore ####
# in-memory writes
#' @name storeWrite-parquetStore
#' @title Write to a Parquet Storage Spec
#' @description
#' Write tabular data to a parquet. Enforces the presence of an integer
#' `row_index` column since parquet does not intrinsically encode row order.
#'
#' The `ANY` signature method is provided for writing lazily accessed data in
#' batches with possible pre-processing via `callback`.
#' @inheritParams storeWrite
#' @param callback `function` (optional). Function to apply to `data` before
#' writing. The first param of this function should accept the `data`.
#' @param row_offset `numeric`, will be coerced to `integer` (default = 0L).
#' Offset to apply to row number indexing. For example `row_offset = 0L` means
#' the first `row_index` value in this file will be `1`.
#' @param ... additional params to pass to [arrow::write_dataset()]
#' @family storeWrite methods
#' @returns A `parquetStore` inheriting class
#' @export
setMethod("storeWrite", signature("parquetStore", "data.frame"),
    function(store, data,
        callback = NULL,
        row_offset = 0L,
        ...) {
    GiottoUtils::package_check("arrow")
    .guard_no_data_rows(data)
    checkmate::assert_function(callback, null.ok = TRUE)
    row_offset <- as.integer(row_offset)
    if (!is.null(callback)) {
        data <- callback(data)
    }
    has_row_index <- "row_index" %in% colnames(data)
    if (!has_row_index) {
        data <- .dt_set_row_index(data, offset = row_offset, col = "row_index")
    }
    store@fields <- colnames(data) # record fields info
    arrow::write_dataset(dataset = data, path = store@path, format = "parquet", ...)
    store
})

# batched writes


#' @rdname storeWrite-parquetStore
#' @export
setMethod("storeWrite", signature("parquetStore", "ANY"), function(store, data,
    callback = NULL,
    row_offset = 0L,
    ...) {
    GiottoUtils::package_check("arrow")
    checkmate::assert_function(callback, null.ok = TRUE)
    row_offset <- as.integer(row_offset)
    .guard_lazy_access(data)
    .guard_no_data_rows(data)

    # if no callback fun, and row_index exists, directly write out
    if (is.null(callback) && "row_index" %in% colnames(data)) {
        store@fields <- colnames(data)
        arrow::write_dataset(dataset = data, path = store@path, format = "parquet", ...)
        return(store)
    }

    # else, if there is something to be batch processed...
    # perform any callbacks and add a row_index if needed
    data <- .arrow_map_batches(data, FUN = callback)
    # if row_index still missing...
    if (!"row_index" %in% colnames(data)) {
        data <- .arrow_add_row_index(data,
            col = "row_index", offset = row_offset
        )
    }
    arrow::write_dataset(dataset = data, path = store@path, format = "parquet", ...)
    store
})

# * parquetGeomStore ####
#' @name storeWrite-parquetGeomStore
#' @title Write to a Parquet Geometry Storage
#' @description
#' Write geometry data to a `parquetGeomStore`. Each row is a single geometry.
#' Enforces the presence of `x_index` and `y_index` cols which are the centroid
#' (as calculated by the center of the geometry's spatial extent) and the
#' `row_index` which preserves row ordering.
#'
#' May be written either directly from terra `SpatVector` or `data.frame` inputs
#' that are first parsed via GiottoClass `createGiottoPoints` or
#' `createGiottoPolygon`, depending on the `type` param.
#' @inheritParams storeWrite
#' @inheritParams storeWrite-parquetStore
#' @param `type` `character`. One of `point` or `polygon` depending on the
#' geometries to create from the table data.
#' @param geom_param `list` (optional). Additional params to pass to
#' [GiottoClass::createGiottoPoints] or [GiottoClass::createGiottoPolygon],
#' depending on `type`.
#' @param ... addtional params to pass to `parquetStore` method and
#' [arrow::write_dataset()]
#' @family storeWrite methods
#' @export
setMethod("storeWrite", signature("parquetGeomStore", "SpatVector"), function(store, data, row_offset = 0, ...) {
    if (nrow(data) == 0L) return(NULL)
    GiottoUtils::package_check("arrow")
    store@extent <- .ext_to_num_vec(ext(data))
    data <- .terra_to_parquet_format(data, row_offset = row_offset)
    store_write_next <- methods::getMethod("storeWrite",
        signature("parquetStore", "data.frame")
    )
    store <- store_write_next(store, data, ...) # call a lower method
    store
})

# parse geometries via {GiottoClass} then pass to further conversion and writes
#' @rdname storeWrite-parquetGeomStore
#' @param type `character`. Either `"point"` or `"polygon"`. What type of
#'   geometry to write. Determines whether [GiottoClass::createGiottoPoints()]
#'   or [GiottoClass::createGiottoPolygon()] is used to parse the values.
#' @export
setMethod("storeWrite", signature("parquetGeomStore", "data.frame"),
    function(store, data,
        type = c("point", "polygon"),
        geom_param = list(verbose = FALSE),
        row_offset = 0,
        ...) {
    if (nrow(data) == 0L) return(NULL)
    GiottoUtils::package_check("arrow")
    checkmate::assert_list(geom_param)
    type <- match.arg(type, choices = c("point", "polygon"))
    fun <- switch(type,
        "point" = GiottoClass::createGiottoPoints,
        "polygon" = GiottoClass::createGiottoPolygon
    )
    # coerce to giotto representation -> SpatVector
    data <- do.call(fun, args = c(list(data), geom_param))
    store <- storeWrite(store, data[], row_offset = row_offset, ...)
    store
})

# * parquetGeomTileStore ####
# this one cannot use data inputs other than fileStore since it is parallelized
#' @name storeWrite-parquetGeomTileStore
#' @title Write to a Parquet Geometry Tiled Storage
#' @description
#' Write geometry data in a spatially tiled manner.
#' @inheritParams storeWrite
#' @inheritParams storeWrite-parquetStore
#' @inheritParams storeWrite-parquetGeomStore
#' @param n_tiles `numeric`. Minimum number of tiles to write across
#' @param tile tilework `tilePlan` (optional)
#' @param sdimx,sdimy `character`. Names of columns containing x and y values,
#'   respectively.
#' @family storeWrite methods
#' @export
setMethod("storeWrite", signature("parquetGeomTileStore", "fileStore"),
    function(store, data,
             n_tiles = 100,
             tile = NULL,
             sdimx = "x", sdimy = "y",
             poly_id = "id",
             type = c("point", "polygon"),
             dry_run = FALSE,
             geom_param = list(verbose = FALSE),
             verbose = NULL,
             ...) {
        GiottoUtils::package_check("arrow")
        a <- storeRead(data)
        if (!inherits(a, c("FileSystemDataset", "tbl_lazy", "arrow_dplyr_query"))) {
            stop("[storeWrite] storeRead(data) should be a `FileSystemDataset` `arrow_dplyr_query`, or `tbl_lazy`")
        }
        if (nrow(a) == 0L) return(NULL)
        type <- match.arg(type, choices = c("point", "polygon"))
        tiles <- store@tiles
        # offsets are needed, so this can't be done by child processes
        envelope <- switch(type,
            "point" = FALSE,
            "polygon" = TRUE
        )
        tiles <- .annotate_tileiterator(tiles,
            data = a,
            n_tiles = n_tiles,
            sdimx = sdimx,
            sdimy = sdimy,
            poly_id = poly_id,
            envelope = envelope
        )
        tile_indices <- tiles$tile[tiles$n_records > 0L]
        vmsg(.v = verbose, sprintf("nonzero tiles: %d out of %d", length(tile_indices), n_tiles))

        if (isTRUE(dry_run)) {
            message("[storeWrite] dry run: planned tiles")
            plot(tiles, values = "n_records", main = "geoms per planned tile")
            return(tiles)
        }
        store@tiles <- tiles
        store@extent <- .ext_to_num_vec(ext(tiles))

        if (!is.null(tile)) tile_indices <- as.integer(tile) # write a specific tile
        written_stores <- lapply_flex(tile_indices, function(tile_i) {
            .pgts_write_tile(
                store = store,
                data = data,
                tile_i = tile_i,
                sdimx = sdimx,
                sdimy = sdimy,
                poly_id = poly_id,
                type = type,
                geom_param = geom_param
            )
        },
        future.seed = TRUE,
        # future.packages = c("terra", "arrow", "data.table"),
        future.globals = list(
            store = store,
            data = data,
            sdimx = sdimx,
            sdimy = sdimy,
            poly_id = poly_id,
            type = type,
            geom_param = geom_param
        ))
        written_stores <- written_stores[!vapply(written_stores, is.null, FUN.VALUE = logical(1L))]
        store@fields <- written_stores[[1]]@fields
        store
})

# this one cannot use data inputs other than fileStore since it is parallelized
#' @rdname storeWrite-parquetGeomTileStore
#' @export
setMethod("storeWrite", signature("parquetGeomTileStore", "fileStore"),
    function(store, data,
             n_tiles = 100,
             tile = NULL,
             sdimx = "x", sdimy = "y",
             poly_id = "id",
             type = c("point", "polygon"),
             dry_run = FALSE,
             geom_param = list(verbose = FALSE),
             verbose = NULL,
             ...) {
        GiottoUtils::package_check("arrow")
        a <- storeRead(data)
        # only permit lazy reading
        if (!inherits(a, c("FileSystemDataset", "tbl_lazy", "arrow_dplyr_query"))) {
            stop("[storeWrite] storeRead(data) should be a `FileSystemDataset` `arrow_dplyr_query`, or `tbl_lazy`")
        }
        if (nrow(a) == 0L) return(NULL)
        type <- match.arg(type, choices = c("point", "polygon"))
        tiles <- store@tiles
        # offsets are needed, so this can't be done by child processes
        envelope <- switch(type,
            "point" = FALSE,
            "polygon" = TRUE
        )
        tiles <- .annotate_tileiterator(tiles,
            data = a,
            n_tiles = n_tiles,
            sdimx = sdimx,
            sdimy = sdimy,
            poly_id = poly_id,
            envelope = envelope
        )
        tile_indices <- tiles$tile[tiles$n_records > 0L]
        vmsg(.v = verbose, sprintf("nonzero tiles: %d out of %d", length(tile_indices), n_tiles))

        if (isTRUE(dry_run)) {
            message("[storeWrite] dry run: planned tiles")
            plot(tiles, values = "n_records", main = "geoms per planned tile")
            return(tiles)
        }
        store@tiles <- tiles
        store@extent <- .ext_to_num_vec(ext(tiles))

        if (!is.null(tile)) tile_indices <- as.integer(tile) # write a specific tile
        written_stores <- lapply_flex(tile_indices, function(tile_i) {
            .pgts_write_tile(
                store = store,
                data = data,
                tile_i = tile_i,
                sdimx = sdimx,
                sdimy = sdimy,
                poly_id = poly_id,
                type = type,
                geom_param = geom_param
            )
        },
        future.seed = TRUE,
        # future.packages = c("terra", "arrow", "data.table"),
        future.globals = list(
            store = store,
            data = data,
            sdimx = sdimx,
            sdimy = sdimy,
            poly_id = poly_id,
            type = type,
            geom_param = geom_param
        ))
        written_stores <- written_stores[!vapply(written_stores, is.null, FUN.VALUE = logical(1L))]
        store@fields <- written_stores[[1]]@fields
        store
})

# * h5ArrayStore ####
#' @rdname storeWrite
#' @export
setMethod("storeWrite", signature("h5ArrayStore", "memoryMatrix"),
    function(store, data, ...) {
        HDF5Array::writeHDF5Array(
            x = data,
            filepath = store@path,
            name = store@name,
            ...
        )
        store
    })

# * tileDBArrayStore ####
#' @rdname storeWrite
#' @export
setMethod("storeWrite", signature("tileDBArrayStore", "memoryMatrix"),
    function(store, data, ...) {
        p <- store@path
        if (!dir.exists(p)) dir.create(p, recursive = TRUE)
        TileDBArray::writeTileDBArray(data,
            path = p,
            ...
        )
        store
    })



# internals ####

.pgts_write_tile <- function(store, data, tile_i,
    sdimx = "x", sdimy = "y", poly_id = "id", type = c("point", "polygon"),
    geom_param = list(verbose = FALSE), ...) {
    if (!inherits(data, "fileStore")) {
        stop("'data' must inherit from `fileStore`\n", call. = FALSE)
    }
    type <- match.arg(type, choices = c("point", "polygon"))
    tile_i = as.integer(tile_i)
    envelope <- switch(type,
        "point" = FALSE,
        "polygon" = TRUE
    )

    # generate arrow file pointer
    a <- storeRead(data)
    # filter down to tile and arrange by id if needed
    tile_data <- .tile_crop(
        tiles = store@tiles,
        data = a,
        i = tile_i,
        sdimx = sdimx,
        sdimy = sdimy,
        group_col = poly_id,
        envelope = envelope
    )
    if (type == "polygon") {
        tile_data <- dplyr::arrange(tile_data, !!as.name(poly_id))
    }
    # pull into memory
    mem_data <- tile_data %>% dplyr::collect()
    tile_store <- storeCreate(
        path = .tilepath(store@path, idx = tile_i),
        type = "parquetGeom"
    )
    # pass to parquetGeomStore, data.frame
    written_store <- storeWrite(tile_store, mem_data,
        type = type,
        geom_param = geom_param,
        row_offset = tiles$row_offset[tile_i],
        ...
    )
    written_store
}

.terra_to_parquet_format <- function(x, row_offset = 0) {
    checkmate::assert_class(x, "SpatVector")
    wkb <- terra::geom(x, wkb = TRUE)
    ctrs <- XY(centroids(x))
    if (!is.matrix(ctrs)) ctrs <- t(as.matrix(ctrs))

    data <- data.frame(
        row_index = seq_len(nrow(ctrs)) + row_offset,
        x_index = ctrs[,1],
        y_index = ctrs[,2]
    )
    data$geom <- wkb # raw needs to be added separately
    data <- cbind(data, terra::values(x))
    data
}

# hive partitioning style
.hive_part_col <- function(col, idx) {
    checkmate::assert_character(col)
    idx <- as.integer(idx)
    sprintf("%s=%03d", col, idx)
}

.hive_part_path <- function(dir, cols, indices) {
    for (i in seq_along(cols)) {
        dir <- file.path(dir, .hive_part_col(cols[[i]], indices[[i]]))
    }
    dir
}

.tilepath <- function(dir, idx) {
    vapply(idx, function(idx_i) {
        .hive_part_path(dir = dir, cols = "tile_index", indices = idx_i)
    }, character(1L))
}


