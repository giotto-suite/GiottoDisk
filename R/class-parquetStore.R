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
#' `[, j]` indexing is supported. `[i, ]` row indexing is not -- use
#' filter/crop-based access instead.
#' 
#' `nrow()` and `dim()` calls return as `numeric` instead of `integer` to
#' support row counts exceeding R's int32 limit (~2.1 billion).
#'
#' @slot path character. Local file path or directory, or remote URI
#'   (s3://, gs://, az://). For remote paths, authentication is handled via
#'   environment variables (AWS_ACCESS_KEY_ID, etc.) or credential files.
#'   Can point to a single file or a directory with optional hive-style
#'   partitioning.
#' @slot fields character. Cached column names from the parquet dataset.
#' @slot read_fun function. Preset to `arrow::open_dataset()` for standard
#'   parquet access. Can be customized for edge cases (see [fileStore-class]).
#' @slot window numeric(4). (`parquetGeomStore`-inheriting) Spatial window as
#'   xmin, xmax, ymin, ymax.
#' @slot geomtype `character` (`parquetGeomStore`-inheriting) Type of geometry
#'   contained (i.e. polygons/points)
#' @slot tiles tileIterator. (parquetGeomTileStore only) \{tilework\} object
#'   defining which tile(s) in a tile plan this store is responsible for.
#' @slot datatype named `list` of `character` (optional). Datatype to read
#'   a column as. These are applied to the arrow schema at dataset opening.
#' @family store types
#' @examples
#' # creation, writing, and existence tests
#' ps <- parquetStore()
#' storeExists(ps)
#' ps <- storeWrite(ps, mtcars)
#' storeExists(ps)
#' 
#' # object shape
#' nrow(ps) # numeric
#' ncol(ps) # integer
#' dim(ps) # numeric
#' 
#' # different outputs
#' storeRead(ps)
#' storeRead(ps, output = "tibble")
#' 
#' # filtering ops
#' ps_1 <- subset(ps, subset = cyl > gear)
#' ps_1
#' storeRead(ps_1, output = "tibble")
#' head(ps_1, 2) |> storeRead(output = 'tibble')
#' # sample size is approximate
#' rowSample(ps_1, 4) |> storeRead(output = "tibble")
NULL

setClass("parquetBase",
    contains = "VIRTUAL",
    slots = list(
        fields = "ANY",
        datatype = "list",
        ops = "list"
    )
)

setClass("parquetGeomBase",
    contains = c("VIRTUAL", "parquetBase"),
    slots = list(
        crop = "numeric",
        window = "numeric",
        geomtype = "character",
        post_ops = "list"
    )
)

# cols: row_index, id, ...
#' @rdname parquetStore-class
setClass("parquetStore",
    contains = c("queryableStore", "parquetBase")
)

# cols: row_index, x_index, y_index, geom, ...
#' @rdname parquetStore-class
setClass("parquetGeomStore",
    contains = c("parquetStore",  "parquetGeomBase")
)

# cols: row_index, x_index, y_index, tile_index, geom, ...
#' @rdname parquetStore-class
setClass("parquetGeomTileStore",
    contains = "parquetGeomStore",
    slots = list(
        tiles = "tilePlan",
        tile_filter = "integer"
    ),
    prototype = list(
        tiles = tilework::tilePlan("spatial"),
        tile_filter = integer(0L)
    )
)

setClass("unionParquetStore",
    contains = "parquetBase",
    slots = list(
        stores = "list",
        params = "list"
    )
)

setClass("unionParquetGeomStore",
    contains = c("unionParquetStore", "parquetGeomBase")
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
parquetStore <- function(path = .dump_tempfile(), ...) {
    new("parquetStore", path = path, ...)
}

#' @rdname parquetStore
#' @export
parquetGeomStore <- function(path = .dump_tempfile(), ...) {
    new("parquetGeomStore", path = path, ...)
}

#' @rdname parquetStore
#' @export
parquetGeomTileStore <- function(path = .dump_tempfile(), ...) {
    new("parquetGeomTileStore", path = path, ...)
}

# initialize ####

setMethod("initialize", signature("parquetStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    # register default read_fun
    .Object@read_fun <- function(x, schema = NULL) {
        arrow::open_dataset(sources = x, schema = schema)
    }
    if (!storeExists(.Object)) return(.Object) # skip if not ready to read
  
    # init cached values (disk state)
    .Object <- .pstore_cache_values(.Object)
  
    # checks
    atypes <- .Object@params$arrow_types
    is_int64 <- atypes == "int64"
    if (any(is_int64)) {
        warning(sprintf("[parquetStore] int64 is not well supported:%s",
            paste0("\n  ", toString(names(atypes[is_int64])))))
    }
  
    if (!"row_index" %in% .pstore_disk_fields(.Object)) {
        stop("[initialize] no 'row_index' column", call. = FALSE)
    }
    .Object
})

setMethod("initialize", signature("parquetGeomStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (!storeExists(.Object)) return(.Object) # skip if not ready to read
    # disk_fields should be cached already in `parquetStore` method
    check_cols <- c("x_index", "y_index", "geom")
    missing_cols <- check_cols[!check_cols %in% .pstore_disk_fields(.Object)]
    if (length(missing_cols) > 0L) {
        warning(sprintf("[initialize] '%s' column(s) should be present",
            toString(missing_cols)), call. = FALSE)
    }
    .Object
})

setMethod("initialize", signature("unionParquetStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@stores) == 0L) return(.Object)
    
    store_exists <- storeExists(.Object, all = FALSE)
    if (!all(store_exists)) {
        paths <- vapply(.Object@stores, function(x) x@path, FUN.VALUE = character(1L))
        uids <- vapply(.Object@stores, function(x) x@uid, FUN.VALUE = character(1L))
        unwritten <- sprintf("[%s]@ %s", uids[!store_exists], paths[!store_exists]) |>
            paste(collapse = "\n")
        stop("[unionParquetStore] not all substores exist:\n", unwritten, call. = FALSE)
    }
  
    .pstore_disk_fields(.Object) <- lapply(.Object@stores,
        .pstore_disk_fields) |>
        unlist() |>
        unique()
  
    atypes <- lapply(.Object@stores,
        .pstore_arrow_types) |>
        unlist()
    anames <- unique(names(atypes))
    .pstore_arrow_types(.Object) <- atypes[anames]
    .Object
})

# internals ####

# central helper for cacheing intialization steps to reduce
# file touches. All values should only require a simple read
# at most, without additional processing
.pstore_cache_values <- function(x) {
    need_fields <- is.null(.pstore_disk_fields(x))
    need_ddtypes <- is.null(.pstore_disk_datatypes(x)) ||
        is.null(.pstore_arrow_types(x))
  
    need_cache <- any(c(need_fields, need_ddtypes))
    if (!need_cache) return(x) # early exit
    
    a <- .store_simple_read(x)
    
    if (need_fields) {
        .pstore_disk_fields(x) <- names(a)
    }
    if (need_ddtypes) {
        # determine R datatypes by pulling in empty table
        preview <- dplyr::collect(head(a, 0L))
        .pstore_disk_datatypes(x) <- lapply(preview, class)
        sc <- arrow::schema(a)
        .pstore_arrow_types(x) <- stats::setNames(vapply(sc$fields,
            function(f) f$type$ToString(), 
            FUN.VALUE = character(1L)
        ), names(sc))
    }
    x
}

# return SpatExtent if a limiting window or crop is set
# return NULL otherwise
.pstore_active_extent <- function(x) {
    if (length(x@window) > 0L) return(ext(x@window))
    if (length(x@crop) > 0L) return(ext(x@crop))
    NULL
}

# cached fields check - output should be disk state
# used for repeated checks across cascading initialize
.pstore_get_fields <- function(x) {
    disk_fields <- .pstore_disk_fields(x) # check cached
    if (!is.null(disk_fields)) return(disk_fields)
  
    names(.store_simple_read(x)) # fallback (direct read from disk)
}

##  parquetStore-specific get/set ####
.pstore_disk_fields <- function(x) x@params$disk_fields
`.pstore_disk_fields<-` <- function(x, value) {
    x@params$disk_fields <- value
    x
}
.pstore_disk_datatypes <- function(x) x@params$disk_datatypes
`.pstore_disk_datatypes<-` <- function(x, value) {
    x@params$disk_datatypes <- value
    x
}
.pstore_arrow_types <- function(x) x@params$arrow_types
`.pstore_arrow_types<-` <- function(x, value) {
    x@params$arrow_types <- value
    x
}
.pstore_disk_extent <- function(x) x@params$disk_extent
`.pstore_disk_extent<-` <- function(x, value) {
    x@params$disk_extent <- value
    x
}

##  parquetGeomStore-specific get/set ####
.pgeom_max_poly_radius <- function(x) x@params$max_poly_radius
`.pgeom_max_poly_radius<-` <- function(x, value) {
    x@params$max_poly_radius <- value
    x
}