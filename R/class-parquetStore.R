#' @include class-dataStore.R
NULL

#' @name parquetStore-class
#' @title Parquet Store
#' @description
#' S4 Class extending [fileStore-class] for indexed storage of tabular data
#' in Apache Parquet format. `parquetStores` provide delayed/query-based
#' access to data rather than loading into memory. They are not intended for
#' in-place updates or edits. Create with [parquetStore()].
#'
#' Stores represent a query interface to on-disk data and do not contain the
#' actual table data as state.
#'
#' @slot path character. Local file path or directory, or remote URI
#'   (s3://, gs://, az://). For remote paths, authentication is handled via
#'   environment variables (AWS_ACCESS_KEY_ID, etc.) or credential files.
#'   Can point to a single file or a directory with optional hive-style
#'   partitioning.
#' @slot fields character. Cached column names from the parquet dataset.
#' @slot read_fun function. Preset to `arrow::open_dataset()` for standard
#'   parquet access. Can be customized for edge cases (see [fileStore-class]).
#' @slot extent numeric(4). (`parquetGeomStore`-inheriting) Spatial extent as
#'   xmin, xmax, ymin, ymax of contained geometries.
#' @slot geomtype `character` (`parquetGeomStore`-inheriting) Type of geometry
#'   contained (i.e. polygons/points)
#' @slot selection `integer`. Tracks row indexing/selection. At the lazy
#'   query level, these are filtered for. When materialized, ordering is also
#'   obeyed.
#' @slot tiles tileIterator. (parquetGeomTileStore only) \{tilework\} object
#'   defining which tile(s) in a tile plan this store is responsible for.
#'
#' @family store types
NULL

# cols: row_index, id, ...
#' @rdname parquetStore-class
setClass("parquetStore",
    contains = "queryableStore",
    slots = list(
        path = "character",
        fields = "character",
        selection = "integer"
    )
)

# cols: row_index, x_index, y_index, geom, ...
#' @rdname parquetStore-class
setClass("parquetGeomStore",
    contains = "parquetStore",
    slots = list(
        extent = "numeric",
        geomtype = "character"  
    )
)

# cols: row_index, x_index, y_index, tile_index, geom, ...
#' @rdname parquetStore-class
setClass("parquetGeomTileStore",
    contains = "parquetGeomStore",
    slots = list(tiles = "tilePlan"),
    prototype = list(
        tiles = tilework::tilePlan("spatial")
    )
)

# constructors ####

#' @name parquetStore
#' @title Create a Parquet Store
#' @description
#' Create a parquet-backed storage.
#' @param path `character`. Disk path to file or hive storage directory
#' @family store constructors
#' @seealso [store][parquetStore-class]
#' @export
parquetStore <- function(path = tempfile(), ...) {
    new("parquetStore", path = path, ...)
}

#' @rdname parquetStore
#' @export
parquetGeomStore <- function(path = tempfile(), ...) {
    new("parquetGeomStore", path = path, ...)
}

#' @rdname parquetStore
#' @export
parquetGeomTileStore <- function(path = tempfile(), ...) {
    new("parquetGeomTileStore", path = path, ...)
}

# initialize ####

setMethod("initialize", signature("parquetStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    # register default read_fun
    .Object@read_fun <- function(x) {                                       
        if (length(x) == 1L) {
            arrow::open_dataset(sources = x)
        } else {
            arrow::open_dataset(sources = lapply(x, arrow::open_dataset))
        }
    }
    if (!.store_exists(.Object)) return(.Object) # skip if not ready to read
    .Object@fields <- .get_fields(.Object) # cache colnames
    if (!"row_index" %in% .Object@fields) {
        stop("[initialize] no 'row_index' column", call. = FALSE)
    }
    .Object
})

setMethod("initialize", signature("parquetGeomStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (!.store_exists(.Object)) return(.Object) # skip if not ready to read
    # fields should be cached already in `parquetStore` method
    check_cols <- c("x_index", "y_index", "geom")
    missing_cols <- check_cols[!check_cols %in% .Object@fields]
    if (length(missing_cols) > 0L) {
        warning(sprintf("[initialize] '%s' column(s) should be present",
            toString(missing_cols)), call. = FALSE)
    }
    if (length(.Object@extent) == 0L &&
        all(c("x_index", "y_index") %in% .Object@fields)) {
        data <- storeRead(.Object)
        .Object@extent <- .ext_to_num_vec(.dplyr_ext(data,
            sdimx = "x_index",
            sdimy = "y_index"
        ))
    }
    .Object
})

# internals ####

# cached fields check
# used for repeated checks across cascading initialize
.get_fields <- function(x) {
    if (length(x@fields) > 0L) return(x@fields)
    colnames(x)
}