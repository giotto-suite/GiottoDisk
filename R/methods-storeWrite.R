
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

.lazy_classes <- c("FileSystemDataset", "tbl_lazy", "arrow_dplyr_query")

.guard_lazy_access <- function(x) {
    if (!inherits(x, .lazy_classes)) {
        stop("[storeWrite] 'data' should be one of: ", .lazy_classes)
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
#'   writing. The first param of this function should accept the `data`.
#' @param row_offset `numeric`, will be coerced to `integer` (default = 0L).
#'   Offset to apply to row number indexing. For example `row_offset = 0L` means
#'   the first `row_index` value in this file will be `1`.
#' @param uid_partition `logical` When `TRUE`, tags the store uid as a column
#'   called `"source_id"` using filepath hive partitioning rules.
#' @param .arrow_meta `named list` of `character` strings to add as parquet
#'   metadata
#' @param ... additional params to pass to [arrow::write_dataset()]
#' @family storeWrite methods
#' @returns A `parquetStore` inheriting class
#' @export
setMethod("storeWrite", signature("parquetStore", "data.frame"),
    function(store, data,
        callback = NULL,
        row_offset = 0L,
        uid_partition = TRUE,
        .arrow_meta = NULL,
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
    # never include source_id col in write
    data <- dplyr::select(data, -dplyr::any_of("source_id"))
    
    # apply any parquet metadata
    if (!is.null(.arrow_meta)) {
        checkmate::assert_list(.arrow_meta)
        data <- arrow::arrow_table(data)
        data$metadata <- c(data$metadata, .arrow_meta)
    }
      
    # apply uid partitioning
    if (isTRUE(uid_partition)) {
        path <- .idpath(store@path, store@uid)
    } else {
        path <- store@path
    }
    arrow::write_dataset(dataset = data, path = path, format = "parquet", ...)
    store
    })

# batched writes

# ANY signature intended for writing lazy queries
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
        # never include source_id col in write
        data <- dplyr::select(data, -dplyr::any_of("source_id"))
        path <- .idpath(store@path, store@uid)
        arrow::write_dataset(dataset = data, path = path, format = "parquet", ...)
        return(store)
    }

    # else, if there is something to be batch processed...
    # perform any callbacks and add a row_index if needed
    if (!is.null(callback)) {
        data <- .arrow_map_batches(data, FUN = callback)
    }
    # if row_index still missing...
    if (!"row_index" %in% colnames(data)) {
        data <- .arrow_add_row_index(data,
            col = "row_index", offset = row_offset
        )
    }
    # never include source_id col in write
    data <- dplyr::select(data, -dplyr::any_of("source_id"))
    path <- .idpath(store@path, store@uid)
    arrow::write_dataset(dataset = data, path = path, format = "parquet", ...)
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
#' @param type `character`. One of `"points"` or `"polygons"` depending on the
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
    store@geomtype <- terra::geomtype(data)
    store@params$crs <- terra::crs(data)
    # convert to data.frame with column 'geom' as WKB
    data <- .terra_to_parquet_format(data, row_offset = row_offset)
    # get schema to apply metadata, append metadata for geom col
    geo_meta <- list(geo = .geoparquet_metadata(
            geom_col = "geom",
            geomtype = store@geomtype,
            crs = store@params$crs,
            extent = store@extent
    ))
    store_write_next <- methods::getMethod("storeWrite",
        signature("parquetStore", "data.frame")
    )
    store <- store_write_next(store, data, .arrow_meta = geo_meta, ...) # call a lower method
    store
})

# parse geometries via {GiottoClass} then pass to further conversion and writes
#' @rdname storeWrite-parquetGeomStore
#' @param type `character`. Either `"points"` or `"polygons"`. What type of
#'   geometry to write. Determines whether [GiottoClass::createGiottoPoints()]
#'   or [GiottoClass::createGiottoPolygon()] is used to parse the values.
#' @export
setMethod("storeWrite", signature("parquetGeomStore", "data.frame"),
    function(store, data,
        type = c("points", "polygons"),
        geom_param = list(verbose = FALSE),
        row_offset = 0,
        ...) {
    if (nrow(data) == 0L) return(NULL)
    GiottoUtils::package_check("arrow")
    checkmate::assert_list(geom_param)
    type <- match.arg(type, choices = c("points", "polygons"))
    fun <- switch(type,
        "points" = GiottoClass::createGiottoPoints,
        "polygons" = GiottoClass::createGiottoPolygon
    )
    # coerce to giotto representation -> SpatVector
    data <- do.call(fun, args = c(list(data), geom_param))
    store <- storeWrite(store, data[], row_offset = row_offset, ...)
    store
})

#' @describeIn storeWrite-parquetGeomStore reads in as `tibble`, then passes
#' to `data.frame` method.
#' @export
setMethod("storeWrite", signature("parquetGeomStore", "parquetStore"), function(store, data, ...) {
    storeWrite(store, storeRead(data, output = "tibble"), ...)
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
#' @param i `integerlike` specific tile(s) to write. Leave as NULL (default)
#'   to write all.
#' @param sdimx,sdimy `character`. Names of columns containing x and y values,
#'   respectively.
#' @param poly_id `character`. If `type = "polygons"`, name of column containing
#'   IDs that group vertices
#' @param contiguous `logical` (default = TRUE) When `TRUE`, tile inclusivity
#' rules (see [getBoundedData]) prevent duplication of data at tile boundaries
#' (assuming no padding).
#' @param dry_run `logical` (default = FALSE). When `TRUE`, stops after planning
#' the tiles to write and plots the tiles to write.
#' @param write_param named `list` (optional). Additional params to pass to the
#' writing function.
#' @param ... additional params to pass to `tileApply`
#' @family storeWrite methods
#' @export
setMethod("storeWrite", signature("parquetGeomTileStore", "queryableStore"),
    function(store, data,
             n_tiles = 100,
             i = NULL,
             sdimx, sdimy, poly_id,
             type = c("points", "polygons"),
             contiguous = TRUE,
             dry_run = FALSE,
             geom_param = list(verbose = FALSE),
             write_param = list(),
             verbose = NULL,
             ...) {
        type <- match.arg(type, choices = c("points", "polygons"))
        poly_stores <- c("parquetStore")
        if (type == "polygons" && !inherits(data, poly_stores)) {
            msg <- c("[storeWrite] not a recognized type to write to 'parquetGeomTileStore'.\n",
                " Consider writing to 'parquetStore' inheriting class first.")
            stop(msg, call. = FALSE)
        }
        GiottoUtils::package_check("arrow")
        a <- storeRead(data)
        # only permit lazy reading
        .guard_lazy_access(a)
        if (nrow(a) == 0L) return(NULL)
        tiles <- store@tiles
        # row offsets are needed -- can't be calculated by child processes
        envelope <- switch(type,
            "points" = FALSE,
            "polygons" = TRUE
        )

        # 1. setup tiles with requested n
        # 2. add $n_records metadata/tile
        tiles <- .tile_annotate(tiles,
            data = a,
            n_tiles = n_tiles,
            sdimx = sdimx,
            sdimy = sdimy,
            poly_id = poly_id,
            envelope = envelope
        )
        # get tile indices with non-zero counts
        tile_indices <- tiles$tile[tiles$n_records > 0L]
        vmsg(.v = verbose, sprintf("nonzero tiles: %d out of %d", length(tile_indices), n_tiles))
        # select non-zero + requested tiles
        # coerces to tiles to `tileSelection`
        if (!is.null(i)) {
            tile_sel <- tiles[i = intersect(tile_indices, as.integer(i)), drop = FALSE]
        } else {
            tile_sel <- tiles[i = tile_indices, drop = FALSE]
        }
        # 'dry_run' stops here. Show planning with a plot
        if (isTRUE(dry_run)) {
            message("[storeWrite] dry run: planned tiles")
            plot(tile_sel, values = "n_records", main = "geoms per planned tile")
            return(tile_sel)
        }
        # update store slots
        store@tiles <- tile_sel@tp
        store@extent <- .ext_to_num_vec(ext(tile_sel@tp))

        # write tile stores
        written_stores <- tileApply(
            x = data,
            tiles = tile_sel,
            contiguous = contiguous, # passes to getTile via ...
            get_params_x = list( # passes to getBoundedData
                sdimx = sdimx,
                sdimy = sdimy,
                group_col = poly_id,
                envelope = envelope
            ),
            FUN = function(tile_atab, .I) {
                # pull into memory
                mem_data <- tile_atab |>
                    dplyr::arrange(source_id, row_index) |>
                    dplyr::collect() |>
                    dplyr::select(-dplyr::any_of(c("row_index", "source_id")))
                tile_store <- storeCreate(
                    path = .tilepath(.idpath(store@path, store@uid), idx = .I),
                    type = "parquetGeom",
                    uid = store@uid # inherit parent uid
                )
                # pass to parquetGeomStore, data.frame
                sw_args <- c(
                    list(tile_store, mem_data, 
                        type = type,
                        geom_param = geom_param,
                        uid_partition = FALSE # would double tag otherwise
                    ), 
                    write_param
                )
                do.call(storeWrite, sw_args)
            },
            ...
        )
        written_stores <- written_stores[!vapply(written_stores, is.null, FUN.VALUE = logical(1L))]
        store@fields <- written_stores[[1]]@fields
        store@params$crs <- written_stores[[1]]@params$crs
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
            name = store@params$name,
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

# * bpcMatrixStore ####
#' @rdname storeWrite
#' @export
setMethod("storeWrite", signature("bpcMatrixStore", "memoryMatrix"),
    function(store, data, ...) {
        p <- store@path
        BPCells::write_matrix_dir(data, dir = p, ...)
        store
    })


# internals ####

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

.idpath <- function(dir, uid) {
    file.path(dir, paste0("source_id=", uid))
}

.tilepath <- function(dir, idx) {
    vapply(idx, function(idx_i) {
        .hive_part_path(dir = dir, cols = "tile_index", indices = idx_i)
    }, character(1L))
}


