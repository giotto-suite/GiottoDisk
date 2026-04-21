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
    stop(
        sprintf("Writing not implemented for store type %s\n", class(store)),
        call. = FALSE
    )
})

# from fileStore ETL chaining
# e.g. data coercible to arrow FileSystemDataset -> <fileStore> -> read to FSD -> <parquetStore> write
#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("fileStore", "fileStore"),
    function(store, data, ...) {
        vmsg(.is_debug = TRUE, "[storeWrite] fileStore,fileStore")
        storeWrite(store, storeRead(data))
    }
)

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
setMethod(
    "storeWrite",
    signature("parquetStore", "data.frame"),
    function(
        store,
        data,
        callback = NULL,
        row_offset = 0L,
        uid_partition = TRUE,
        tile_idx = NULL,
        .arrow_meta = NULL,
        ...
    ) {
        vmsg(.is_debug = TRUE, "[storeWrite] parquetStore,data.frame")
        GiottoUtils::package_check("arrow")
        .guard_no_data_rows(data)
        checkmate::assert_function(callback, null.ok = TRUE)
        row_offset <- as.integer(row_offset)
        if (!is.null(callback)) {
            data <- callback(data)
        }
        has_row_index <- "row_index" %in% colnames(data)
        if (!has_row_index) {
            data <- .dt_set_row_index(
                data,
                offset = row_offset,
                col = "row_index"
            )
        }
        # never include source_id col in write
        data <- dplyr::select(data, -dplyr::any_of("source_id"))

        # apply any parquet metadata
        if (!is.null(.arrow_meta)) {
            checkmate::assert_list(.arrow_meta)
            data <- arrow::arrow_table(data)
            data$metadata <- c(data$metadata, .arrow_meta)
        }
      
        .write_parquet(store, data, uid_partition = uid_partition, tile_idx = tile_idx, ...)
        initialize(store)
    }
)

# batched writes

# ANY signature intended for writing lazy queries
#' @rdname storeWrite-parquetStore
#' @export
setMethod(
    "storeWrite",
    signature("parquetStore", "ANY"),
    function(store, data, 
        callback = NULL,
        row_offset = 0L,
        uid_partition = TRUE,
        ...) {
        vmsg(.is_debug = TRUE, "[storeWrite] parquetStore,ANY")
        GiottoUtils::package_check("arrow")
        checkmate::assert_function(callback, null.ok = TRUE)
        row_offset <- as.integer(row_offset)
        .guard_lazy_access(data)
        .guard_no_data_rows(data)
      
        # if no callback fun, and row_index exists, directly write out
        if (is.null(callback) && "row_index" %in% colnames(data)) {
            # never include source_id col in write
            data <- dplyr::select(data, -dplyr::any_of("source_id"))
            .write_parquet(store, data, uid_partition = uid_partition, ...)
            return(initialize(store))
        }

        # else, if there is something to be batch processed...
        # perform any callbacks and add a row_index if needed
        if (!is.null(callback)) {
            data <- .arrow_map_batches(data, FUN = callback)
        }
        # if row_index still missing...
        if (!"row_index" %in% colnames(data)) {
            data <- .arrow_add_row_index(
                data,
                col = "row_index",
                offset = row_offset
            )
        }

        # never include source_id col in write
        data <- dplyr::select(data, -dplyr::any_of("source_id"))
        .write_parquet(store, data, uid_partition = uid_partition, ...)
        initialize(store)
    }
)

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
#' @param meta (optional) `data.frame-like` that contains attributes info in the
#'   same ordering as the spatvector in `data`. If not provided and attributes
#'   are in the spatvector, they will be automatically extracted with
#'   [terra::values()]
#' @param tile_idx `integer(1)`. Tile index to write under as a hive partition
#'   (`tile_index=<n>/`). Defaults to `0L` so flat geom stores share a uniform
#'   `(tile_index, row_index)` join key with tiled stores without needing
#'   `dplyr::coalesce`. Pass `NULL` to write without a tile subdir (not
#'   recommended for `parquetGeomStore`).
#' @param ... addtional params to pass to `parquetStore` method and
#' [arrow::write_dataset()]
#' @family storeWrite methods
#' @export
setMethod(
    "storeWrite",
    signature("parquetGeomStore", "SpatVector"),
    function(store, data, meta = NULL, row_offset = 0, tile_idx = 0L, .arrow_meta = NULL, ...) {
        if (nrow(data) == 0L) {
            return(NULL)
        }
        vmsg(.is_debug = TRUE, "[storeWrite] parquetGeomStore,SpatVector")
        GiottoUtils::package_check("arrow")
        store@params$disk_extent <- .ext_to_num_vec(ext(data))
        store@geomtype <- terra::geomtype(data)
        store@params$crs <- terra::crs(data)
        if (store@geomtype == "polygons") {
            .pgeom_max_poly_radius(store) <- sqrt(
                suppressWarnings(max(terra::expanse(data), na.rm = TRUE)) / pi
            )
        }
        # get meta if not already separate
        if (is.null(meta)) {
            if (ncol(data) > 0L) {
                meta <- terra::values(data)
            } else {
                meta <- data.frame()
            }
        }
        # convert to data.frame with column 'geom' as WKB
        data <- .terra_to_parquet_format(data, 
            meta = meta, 
            row_offset = row_offset
        )
        # add geoparquet metadata if needed
        .arrow_meta <- .arrow_meta_add_geoparquet(store, .arrow_meta)

        store_write_next <- methods::getMethod(
            "storeWrite",
            signature("parquetStore", "data.frame")
        )
        store <- store_write_next(store, data,
            tile_idx = tile_idx,
            .arrow_meta = .arrow_meta,
            ...
        )
        store
    }
)

# parse geometries via {GiottoClass} then pass to further conversion and writes
#' @rdname storeWrite-parquetGeomStore
#' @param type `character`. Either `"points"` or `"polygons"`. What type of
#'   geometry to write.
#' @param id_col `character`. Column name in `data` that holds the geometry
#'   identifier. For polygons this groups vertices into individual geometries;
#'   for points this is the feature name. Regardless of the input name, the
#'   column is written to disk as `poly_ID` (polygons) or `feat_ID` (points)
#'   so downstream consumers (e.g. [calculateOverlap()]) can assume standardized
#'   names.
#' @param sdimx `character`. Column name containing x spatial coordinates.
#' @param sdimy `character`. Column name containing y spatial coordinates.
#' @param part_col `character` (optional, polygons only). Column name for
#'   multi-part polygon grouping. See [GiottoClass::createGiottoPolygon].
#' @export
setMethod(
    "storeWrite",
    signature("parquetGeomStore", "data.frame"),
    function(
        store,
        data,
        type = c("points", "polygons"),
        id_col,
        sdimx,
        sdimy,
        part_col = NULL,
        row_offset = 0,
        ...
    ) {
        vmsg(.is_debug = TRUE, "[storeWrite] parquetGeomStore,data.frame")
        if (nrow(data) == 0L) {
            return(NULL)
        }
        GiottoUtils::package_check("arrow")
        geom_param <- list(
            id_col = id_col,
            sdimx = sdimx,
            sdimy = sdimy,
            verbose = FALSE
        )
        type <- match.arg(type, choices = c("points", "polygons"))
        if (type == "polygons") {
            geom_param$part_col <- part_col
        }
        fun <- switch(type,
            "points" = .df_to_terra_pts,
            "polygons" = .df_to_terra_poly
        )
        # coerce to giotto representation -> SpatVector
        out <- do.call(fun, args = c(list(data), geom_param))
        # standardize ID column name on disk regardless of user-supplied id_col
        standard_id <- if (type == "polygons") "poly_ID" else "feat_ID"
        if (id_col != standard_id && id_col %in% names(out$meta)) {
            data.table::setnames(out$meta, id_col, standard_id)
        }
        store <- storeWrite(store,
            data = out$geom,
            meta = out$meta,
            row_offset = row_offset,
            ...
        )
        store
    }
)

setMethod(
    "storeWrite",
    signature("parquetStore", "parquetStore"),
    function(store, data, ...) {
        vmsg(.is_debug = TRUE, "[storeWrite] parquetStore,parquetStore")
        store@params <- data@params
        # remove items that need to be regenerated
        .pstore_disk_fields(store) <- NULL

        arrange_cols <- c("source_id", "tile_index", "row_index")
        drop_cols <- c("source_id", "tile_index", "row_index")

        atab <- storeRead(data, output = "query") |>
            dplyr::arrange(dplyr::across(dplyr::any_of(arrange_cols))) |>
            dplyr::select(-dplyr::any_of(drop_cols))
        # dispatch to parquetStore/ANY which regenerates the row_index
        storeWrite(store, atab, ...)
    }
)

setMethod(
    "storeWrite",
    signature("parquetGeomStore", "parquetGeomStore"),
    function(store, data, ...) {
        vmsg(.is_debug = TRUE, "[storeWrite] parquetGeomStore,parquetGeomStore")
        # carry spatial metadata from source
        store@geomtype <- data@geomtype
        callNextMethod(store, data, ...)
    }
)

# #' @describeIn storeWrite-parquetGeomStore reads in as `tibble`, then passes
# #' to `data.frame` method.
# #' @export
# setMethod("storeWrite", signature("parquetGeomStore", "parquetStore"), function(store, data, ...) {
#     vmsg(.is_debug = TRUE, "[storeWrite] parquetGeomStore,parquetStore")
#     storeWrite(store, storeRead(data, output = "tibble"), ...)
# })

# * parquetGeomTileStore ####
# this one cannot use data inputs other than fileStore since it is parallelized
#' @name storeWrite-parquetGeomTileStore
#' @title Write to a Parquet Geometry Tiled Storage
#' @description
#' Write geometry data to hive-partitoned spatially tiled parquets. The partition
#' column is `"tile_index"`.
#' 
#' ## Accepted `data` inputs
#' 
#' * `queryableStore`-inheriting - Any dplyr lazy queryable store. Requires
#'   selection of `type` param. 
#'   * `"points"` - written directly.
#'   * `"polygons"` - must already be written to `"parquetStore"`-inheriting
#'     class. Computes geometries from vertices which requires row ordering.
#'   
#' * `parquetGeomStore`-inheriting - Simple write to a new tiled store.
#' @inheritParams storeWrite
#' @inheritParams storeWrite-parquetStore
#' @inheritParams storeWrite-parquetGeomStore
#' @param threshold `numeric` or `NULL`. Maximum number of records per tile.
#'   `NULL` (default) auto-selects based on total row count and `type`:
#'   * `"points"`: 1M/tile up to 1B rows; clamped ramp `n/1000` (1M--100M) above.
#'   * `"polygons"` (vertex rows): 500k/tile up to 500M rows;
#'     clamped ramp `n/1000` (500k--5M) above that.
#' @param tiles `freeTilePlan`, `tilePlan`, or `NULL`. When a `freeTilePlan`
#'   is provided (e.g. the output of a `dry_run`), it is used directly without
#'   replanning. When a `tilePlan` is provided, it is used as the seed grid for
#'   [tilework::quadtreePlan()]. When `NULL`, a default seed grid is built from
#'   the data extent's aspect ratio.
#' @param i `integerlike` specific tile(s) to write. Leave as NULL (default)
#'   to write all.
#' @param sdimx,sdimy `character`. Names of columns containing x and y spatial
#'   coordinate values, respectively. Used for spatial tile planning and
#'   geometry construction.
#' @param group_col `character`. If `type = "polygons"`, name of column used
#'   for spatial tile planning via envelope centroids. Defaults to `id_col`.
#' @param contiguous `logical` (default = TRUE) When `TRUE`, tile inclusivity
#' rules (see [getBoundedData]) prevent duplication of data at tile boundaries
#' (assuming no padding).
#' @param dry_run `logical` (default = FALSE). When `TRUE`, runs tile planning
#'   and plots the adaptive tile layout with record counts, then returns the
#'   `freeTilePlan` invisibly without writing for inspection. This can then be
#'   passed to `tiles` to skip tile planning.
#' @param write_param named `list` (optional). Additional params to pass to the
#' writing function.
#' @param ... additional params to pass to `tileApply`
#' @family storeWrite methods
NULL

#' @name storeWrite-parquetGeomTileStore
#' @export
setMethod(
    "storeWrite",
    signature("parquetGeomTileStore", "queryableStore"),
    function(
        store,
        data,
        threshold = NULL,
        tiles = NULL,
        i = NULL,
        id_col,
        sdimx,
        sdimy,
        type = c("points", "polygons"),
        contiguous = TRUE,
        dry_run = FALSE,
        write_param = list(),
        verbose = NULL,
        ...
    ) {
        vmsg(
            .is_debug = TRUE,
            "[storeWrite] parquetGeomTileStore,queryableStore"
        )
        type <- match.arg(type, choices = c("points", "polygons"))
        if (type == "polygons") {
            # polygons vertices need exact row ordering. A correct `row_index` is needed.
            msg <- c(
                "[storeWrite] not a recognized type to write to 'parquetGeomTileStore'.\n",
                " Consider writing to 'parquetStore' inheriting class first."
            )
            stop(msg, call. = FALSE)
        }
        geom_param <- list(
            id_col = id_col,
            sdimx = sdimx,
            sdimy = sdimy
        )

        a <- storeRead(data)
        .guard_lazy_access(a) # only permit lazy reading
        if (nrow(a) == 0L) {
            return(NULL)
        }

        n <- .dplyr_nrow(a)
        threshold <- threshold %||% .auto_threshold(n, type = type)
        if (is.null(tiles)) {
            tiles <- tilework::tilePlan("spatial")
            data_ext <- .dplyr_ext(a, sdimx = sdimx, sdimy = sdimy)
            terra::ext(tiles) <- data_ext
            erange <- range(data_ext)
            length(tiles) <- round(max(erange) / min(erange)) * 4L
        }
        plans <- .tile_plan_quadtree(store, data,
            threshold = threshold,
            tiles = tiles,
            i = i,
            dry_run = dry_run,
            verbose = verbose,
            qtree_args = list(
                sdimx = sdimx,
                sdimy = sdimy,
                contiguous = contiguous
            )
        )
        if (isTRUE(dry_run)) return(invisible(plans$tile_sel@tp))
        store <- plans$store
        tile_sel <- plans$tile_sel
        store@geomtype <- type

        # write tile stores
        written_stores <- tileApply(
            x = data,
            tiles = tile_sel,
            contiguous = contiguous, # passes to getTile via ...
            get_params_x = list(
                # passes to getBoundedData
                sdimx = sdimx,
                sdimy = sdimy,
                envelope = FALSE
            ),
            FUN = function(tile_atab, .I) {
                # pull into memory
                mem_data <- tile_atab |>
                    dplyr::arrange(dplyr::across(dplyr::any_of(c("source_id", "row_index")))) |>
                    dplyr::collect() |>
                    dplyr::select(-dplyr::any_of(c("row_index", "source_id")))
                tile_store <- storeCreate(
                    path = .tilepath(.idpath(store@path, store@uid), idx = .I),
                    type = "parquetGeom",
                    uid = store@uid # inherit parent uid
                )
                # pass to parquetGeomStore, data.frame
                sw_args <- c(
                    list(
                        tile_store,
                        mem_data,
                        type = type,
                        uid_partition = FALSE, # would double tag otherwise
                        .arrow_meta = .arrow_meta_add_geoparquet(store)
                    ),
                    geom_param,
                    write_param
                )
                do.call(storeWrite, sw_args)
            },
            ...
        )
        written_stores <- written_stores[
            !vapply(written_stores, is.null, FUN.VALUE = logical(1L))
        ]
        store@params$crs <- written_stores[[1]]@params$crs
        if (type == "polygons") {
            poly_radii <- vapply(written_stores, function(s) {
                .pgeom_max_poly_radius(s) %||% NA_real_
            }, numeric(1L))
            if (any(!is.na(poly_radii))) {
                .pgeom_max_poly_radius(store) <- max(poly_radii, na.rm = TRUE)
            }
        }
        initialize(store)
    }
)

#' @name storeWrite-parquetGeomTileStore
#' @export
setMethod(
    "storeWrite",
    signature("parquetGeomTileStore", "parquetStore"),
    function(
        store,
        data,
        threshold = NULL,
        tiles = NULL,
        i = NULL,
        id_col,
        sdimx,
        sdimy,
        part_col = NULL,
        group_col = id_col,
        type = c("points", "polygons"),
        contiguous = TRUE,
        dry_run = FALSE,
        write_param = list(),
        verbose = NULL,
        ...
    ) {
        vmsg(
            .is_debug = TRUE,
            "[storeWrite] parquetGeomTileStore,parquetStore"
        )
        type <- match.arg(type, choices = c("points", "polygons"))
        geom_param <- list(
            id_col = id_col,
            sdimx = sdimx,
            sdimy = sdimy,
            part_col = part_col
        )

        a <- storeRead(data)
        .guard_lazy_access(a)
        if (nrow(a) == 0L) {
            return(NULL)
        }
        envelope <- type != "points"

        n <- .dplyr_nrow(a)
        threshold <- threshold %||% .auto_threshold(n, type = type)
        if (is.null(tiles)) {
            tiles <- tilework::tilePlan("spatial")
            data_ext <- .dplyr_ext(a, sdimx = sdimx, sdimy = sdimy)
            terra::ext(tiles) <- data_ext
            erange <- range(data_ext)
            length(tiles) <- round(max(erange) / min(erange)) * 4L
        }
        plans <- .tile_plan_quadtree(store, data,
            threshold = threshold,
            tiles = tiles,
            i = i,
            dry_run = dry_run,
            verbose = verbose,
            qtree_args = list(
                sdimx = sdimx,
                sdimy = sdimy,
                contiguous = contiguous
            )
        )
        if (isTRUE(dry_run)) return(invisible(NULL))
        store <- plans$store
        tile_sel <- plans$tile_sel
        store@geomtype <- type

        # write tile stores
        written_stores <- tileApply(
            x = data,
            tiles = tile_sel,
            contiguous = contiguous, # passes to getTile via ...
            get_params_x = list(
                # passes to getBoundedData
                sdimx = sdimx,
                sdimy = sdimy,
                group_col = group_col,
                envelope = envelope
            ),
            FUN = function(tile_atab, .I) {
                # pull into memory
                mem_data <- tile_atab |>
                    dplyr::arrange(dplyr::across(dplyr::any_of(c("source_id", "row_index")))) |>
                    dplyr::collect() |>
                    dplyr::select(-dplyr::any_of(c("row_index", "source_id")))
                tile_store <- storeCreate(
                    path = .tilepath(.idpath(store@path, store@uid), idx = .I),
                    type = "parquetGeom",
                    uid = store@uid # inherit parent uid
                )
                # pass to parquetGeomStore, data.frame
                sw_args <- c(
                    list(
                        tile_store,
                        mem_data,
                        type = type,
                        uid_partition = FALSE, # would double tag otherwise
                        .arrow_meta = .arrow_meta_add_geoparquet(store)
                    ),
                    geom_param,
                    write_param
                )
                do.call(storeWrite, sw_args)
            },
            ...
        )
        written_stores <- written_stores[
            !vapply(written_stores, is.null, FUN.VALUE = logical(1L))
        ]
        store@params$crs <- written_stores[[1]]@params$crs
        if (type == "polygons") {
            poly_radii <- vapply(written_stores, function(s) {
                .pgeom_max_poly_radius(s) %||% NA_real_
            }, numeric(1L))
            if (any(!is.na(poly_radii))) {
                .pgeom_max_poly_radius(store) <- max(poly_radii, na.rm = TRUE)
            }
        }
        initialize(store)
    }
)

#' @name storeWrite-parquetGeomTileStore
#' @export
setMethod(
    "storeWrite",
    signature("parquetGeomTileStore", "parquetGeomStore"),
    function(
        store,
        data,
        threshold = NULL,
        tiles = NULL,
        i = NULL,
        contiguous = TRUE,
        dry_run = FALSE,
        write_param = list(),
        verbose = NULL,
        ...
    ) {
        vmsg(
            .is_debug = TRUE,
            "[storeWrite] parquetGeomTileStore,parquetGeomStore"
        )

        a <- storeRead(data)
        .guard_lazy_access(a) # only permit lazy reading
        if (nrow(a) == 0L) {
            return(NULL)
        }

        threshold <- threshold %||% .auto_threshold(.dplyr_nrow(a), type = data@geomtype)
        plans <- .tile_plan_quadtree(store, data,
            threshold = threshold,
            tiles = tiles,
            i = i,
            dry_run = dry_run,
            verbose = verbose,
            qtree_args = list(
                sdimx = "x_index",
                sdimy = "y_index",
                contiguous = contiguous
            )
        )
        if (isTRUE(dry_run)) return(invisible(NULL))
        store <- plans$store
        tile_sel <- plans$tile_sel
        store@geomtype <- data@geomtype
      
        # write tile stores
        written_stores <- tileApply(
            x = data,
            tiles = tile_sel,
            contiguous = contiguous, # passes to getTile via ...
            get_params_x = list(
                # passes to getBoundedData
                sdimx = "x_index",
                sdimy = "y_index",
                envelope = FALSE
            ),
            FUN = function(tile_atab, .I) {
                tile_atab <- tile_atab |>
                    dplyr::arrange(dplyr::across(dplyr::any_of(c("source_id", "row_index")))) |>
                    dplyr::select(-dplyr::any_of(c("row_index", "source_id", "tile_index")))
                tile_store <- storeCreate(
                    path = .tilepath(.idpath(store@path, store@uid), idx = .I),
                    type = "parquetGeom",
                    uid = store@uid # inherit parent uid
                )
                # pass to parquetGeomStore, ANY (lazy)
                sw_args <- c(
                    list(
                        tile_store,
                        tile_atab,
                        uid_partition = FALSE # would double tag otherwise
                    ),
                    write_param
                )
                do.call(storeWrite, sw_args)
            },
            ...
        )
        written_stores <- written_stores[
            !vapply(written_stores, is.null, FUN.VALUE = logical(1L))
        ]
        store@params$crs <- data@params$crs
        if (!is.null(.pgeom_max_poly_radius(data))) {
            .pgeom_max_poly_radius(store) <- .pgeom_max_poly_radius(data)
        }
        initialize(store)
    }
)

# * h5ArrayStore ####
#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("h5ArrayStore", "memoryMatrix"),
    function(store, data, ...) {
        HDF5Array::writeHDF5Array(
            x = data,
            filepath = store@path,
            name = store@params$name,
            ...
        )
        store
    }
)


# * tileDBArrayStore ####
#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("tileDBArrayStore", "memoryMatrix"),
    function(store, data, ...) {
        p <- store@path
        if (!dir.exists(p)) {
            dir.create(p, recursive = TRUE)
        }
        TileDBArray::writeTileDBArray(data, path = p, ...)
        store
    }
)

# * bpcMatrixStore ####
#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("bpcMatrixStore", "memoryMatrix"),
    function(store, data, ...) {
        p <- store@path
        BPCells::write_matrix_dir(data, dir = p, ...)
        store
    }
)


# internals ####

.auto_threshold <- function(n, type = c("points", "polygons")) {
    type <- match.arg(type)
    if (type == "points") {
        if (n < 1e6) return(Inf)
        if (n < 1e9) return(1e6)
        return(min(max(n / 1000, 1e6), 1e8))
    } else {
        if (n < 5e5) return(Inf)
        if (n < 5e8) return(5e5)
        return(min(max(n / 1000, 5e5), 5e6))
    }
}

.tile_plan_quadtree <- function(store, data, threshold, tiles, i, dry_run, verbose, qtree_args) {
    if (inherits(tiles, "freeTilePlan")) {
        fp <- tiles # use directly -- already a finished plan (e.g. from dry_run)
    } else {
        fp <- do.call(tilework::quadtreePlan, c(
            list(data, tiles = tiles, threshold = threshold),
            qtree_args
        ))
    }

    # filter to non-empty tiles, then intersect with requested i
    nonempty <- which(fp@metadata$n_records > 0L)
    tile_indices <- if (!is.null(i)) intersect(nonempty, as.integer(i)) else nonempty

    vmsg(.v = verbose,
        sprintf("nonzero tiles: %d out of %d", length(tile_indices), length(fp)))

    tile_sel <- fp[i = tile_indices, drop = FALSE]

    if (isTRUE(dry_run)) {
        message("[storeWrite] dry run: planned tiles")
        plot(tile_sel, values = "n_records", main = "geoms per planned tile")
        return(invisible(tile_sel))
    }

    b <- fp@bounds
    disk_ext <- as.numeric(c(min(b[, 1L]), max(b[, 2L]), min(b[, 3L]), max(b[, 4L])))
    store@tiles <- fp
    .pstore_disk_extent(store) <- disk_ext
    list(tile_sel = tile_sel, store = store)
}

.write_parquet <- function(store, data, uid_partition = TRUE, tile_idx = NULL, ...) {
    if (isTRUE(uid_partition)) {
        path <- .idpath(store@path, store@uid)
    } else {
        path <- store@path
    }
    if (!is.null(tile_idx)) {
        path <- file.path(path, .hive_part_col("tile_index", tile_idx))
    }
    arrow::write_dataset(
        dataset = data,
        path = path,
        format = "parquet",
        ...
    )
}

.terra_to_parquet_format <- function(x, meta, row_offset = 0L) {
    checkmate::assert_class(x, "SpatVector")
    wkb <- terra::geom(x, wkb = TRUE)
    ctrs <- XY(centroids(x))
    if (!is.matrix(ctrs)) {
        ctrs <- t(as.matrix(ctrs))
    }

    data <- data.frame(
        row_index = as.integer(seq_len(nrow(ctrs)) + row_offset),
        x_index = ctrs[, 1],
        y_index = ctrs[, 2]
    )
    data$geom <- wkb # raw needs to be added separately
    if (ncol(meta) > 0L) {
        data <- cbind(data, meta)
    }
    data
}

# generate the geoparquet metadata for a geom store if needed
.arrow_meta_add_geoparquet <- function(store, .arrow_meta = NULL, extent = NULL) {
    checkmate::assert_class(store, "parquetGeomStore")
    checkmate::assert_list(.arrow_meta, null.ok = TRUE)
    
    if (is.null(.arrow_meta)) {
        .arrow_meta <- list()
    }
    # append to existing metadata if any/modify if extent specified
    if (is.null(.arrow_meta$geo) || !is.null(extent)) {
        extent <- extent %||% store@params$disk_extent
        geo_meta <- list(
            geo = .geoparquet_metadata(
                geom_col = "geom",
                geomtype = store@geomtype,
                crs = store@params$crs %||% "",
                extent = extent
            )
        )
        .arrow_meta <- c(.arrow_meta, geo_meta)
    }
    .arrow_meta
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
    vapply(
        idx,
        function(idx_i) {
            .hive_part_path(dir = dir, cols = "tile_index", indices = idx_i)
        },
        character(1L)
    )
}
