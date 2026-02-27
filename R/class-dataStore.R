#' @include utils.R
#' @include pkg_imports.R
NULL

# docs ####

#' @name dataStore-class
#' @aliases store
#' @title Data Storage
#' @description
#' S4 class define a storage of data. The base class `dataStore` is VIRTUAL
#' @family store types
NULL

#' @name fileStore-class
#' @title File Store
#' @description
#' S4 class for disk-backed data storage. Defines a file system `path` and
#' a `read_fun` method for accessing the data. Create with [fileStore()].
#'
#' `fileStore` serves dual purposes:
#' * base class for disk-backed stores
#' * wild card customizable class allowing ad-hoc compatibility with
#'   [storeRead()] for non-standard file formats.
#'
#' @section Extending Store Classes:
#' Specific stores with particular **file formats** (parquet, HDF5, etc.)
#' extend this class. These specialized stores provide preset `read_fun`
#' implementations in their `initialize()` methods, but the ability to
#' customize remains available for edge cases.
#'
#' When extending, register a lightweight `read_fun` to *access* the data,
#' while specific optimization on reading/filtering should be in the
#' `storeRead()` method.
#'
#' @section Customizing Access:
#' The `read_fun` is primarily useful for custom file formats:
#'
#' ```r
#' # Wildcard usage - custom formats
#' rds_store <- fileStore(path = "data.rds", read_fun = readRDS)
#' fst_store <- fileStore(path = "data.fst", read_fun = fst::read_fst)
#' ```
#'
#' Specialized stores (parquetStore, h5ArrayStore, etc.) have preset readers
#' but can be overridden for edge cases if needed.
#'
#' @slot path `character`. File system path (local or remote URI)
#' @slot uid `character` (optional). unique ID for artifact tracking
#' @slot params `list` Additional params such as remote names that are 
#' relevant to the store type may be stored here as a named list.
#' @slot read_fun `function`. Function to access data where `read_fun(path)`
#'   returns the data in a useful format. Can be customized for
#'   compression or other access requirements. Do not embed credentials in
#'   `read_fun()`
#' @family store types
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
#' @slot extent numeric(4). (parquetGeomStore only) Spatial extent as
#'   xmin, xmax, ymin, ymax of contained geometries.
#' @slot tiles tileIterator. (parquetGeomTileStore only) \{tilework\} object
#'   defining which tile(s) in a tile plan this store is responsible for.
#'
#' @family store types
NULL







# definitions ####

setOldClass("matrix")
setOldClass("data.frame")
setOldClass("data.table")
setClassUnion("memoryMatrix", c("matrix", "Matrix"))
setClassUnion("memoryStore", c("memoryMatrix", "data.frame", "data.table"))

#' @rdname dataStore-class
setClass("dataStore", contains = "VIRTUAL")
#' @rdname fileStore-class
setClass("fileStore",
    contains = "dataStore",
    slots = list(
        path = "character",
        uid = "character",
        read_fun = "function",
        params = "list"
    )
)

#' @rdname fileStore-class
#' @section `queryableStore`:
#' Subclass of `fileStore` that is accessible via dplyr query semantics.
#' Flag that the `fileStore` is queryable by coercing to this class
#' ```r
#' as(x, "queryableStore")
#' ```
setClass("queryableStore",
    contains = "fileStore"
)

# * parquet ####

# cols: row_index, id, ...
#' @rdname parquetStore-class
setClass("parquetStore",
    contains = "queryableStore",
    slots = list(
        path = "character",
        fields = "character"
    )
)

# cols: row_index, x_index, y_index, geom, id, ...
#' @rdname parquetStore-class
setClass("parquetGeomStore",
    contains = "parquetStore",
    slots = list(extent = "numeric")
)

# cols: row_index, x_index, y_index, tile_index, geom, id, ...
#' @rdname parquetStore-class
setClass("parquetGeomTileStore",
    contains = "parquetGeomStore",
    slots = list(tiles = "tilePlan"),
    prototype = list(
        tiles = tilework::tilePlan("spatial")
    )
)

# * coercion ####

setAs("fileStore", "queryableStore", function(from) {
    qs <- new("queryableStore",
        path = from@path,
        uid = from@uid,
        read_fun = from@read_fun,
        params = from@params
    )
    access <- storeRead(qs)
    .guard_lazy_access(access)
    qs
})

# * delayed ####

setClass("h5ArrayStore",
    contains = "fileStore"
)

setClass("tileDBArrayStore",
    contains = "fileStore"
)

setClass("bpcMatrixStore",
    contains = "fileStore"
)

# constructors ####

#' @name fileStore
#' @title Create a File Store
#' @description
#' Create a `fileStore` to represent a non-standard file on disk. A `read_fun`
#' function must be provided where the first param of the function can read
#' the `path` to load in the data in a useful format.
#' @param path `character`. Disk path to file
#' @param read_fun `function`. Used to read in the data based on the path.
#' @family store constructors
#' @export
fileStore <- function(path = tempfile(), read_fun, ...) {
    if (missing(read_fun)) stop(call. = FALSE,
        "[fileStore] `read_fun` must be provided.\n"
    )
    checkmate::assert_function(read_fun)
    fs <- new("fileStore", path = path, read_fun = read_fun, ...)
    if (.is_empty_fun(fs@read_fun)) {
        stop("[fileStore] `fileStore` requires a `read_fun`\n",
             call. = FALSE)
    }
    fs
}

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

#' @name arrayStore
#' @title Array Storage
#' @param path storage directory
#' @param name `character` Remote naming used by some storage types,
#'   for example HDF5Array and TileDB.
#' @export
h5ArrayStore <- function(
    path = tempfile(), name = HDF5Array::getHDF5DumpName(), ...) {
    package_check("HDF5Array", repository = "Bioc")
    new("h5ArrayStore",
        path = path, 
        params = list("name" = name), 
        ...
    )
}

#' @rdname arrayStore
#' @export
tileDBArrayStore <- function(
        path = file.path(tempdir(), .make_uid()),
        name = TileDBArray::getTileDBAttr(),
        ...
    ) {
    package_check("TileDBArray", repository = "Bioc")
    new("tileDBArrayStore", 
        path = path,
        params = list("name" = name),
        ...
    )
}

#' @rdname arrayStore
#' @export
bpcMatrixStore <- function(
    path = file.path(tempdir(), .make_uid()),
    ...
) {
    package_check(
        pkg_name = "BPCells", 
        repository = "github:bnprks/BPCells/r"
    )
    new("bpcMatrixStore",
        path = path, 
        ...
    )
}

#' @name storeCreate
#' @title Create a Store
#' @description
#' Hub function for creating a concrete store class.
#' @param path `character`. Disk path to file or hive storage directory
#' @param type `character`. Type of store to create. Currently one of:
#' 
#' * `"file"`
#' * `"parquet"`
#' * `"parquetGeom"`
#' * `"parquetGeomTile"`
#' * `"h5"`
#' * `"bpcells"`
#' * `"tiledb"`
#' @family store constructors
#' @seealso [store]
#' @export
storeCreate <- function(path = tempfile(), type = "parquet", ...) {
    type <- match.arg(tolower(type),
        choices = c("parquet", "parquetgeom", "parquetgeomtile", "file", "h5", "bpcells", "tiledb")
    )
    store <- switch(type,
        "parquet" = parquetStore(path = path, ...),
        "parquetgeom" = parquetGeomStore(path = path, ...),
        "parquetgeomtile" = parquetGeomTileStore(path = path, ...),
        "file" = fileStore(path = path, ...),
        "h5" = h5ArrayStore(path = path, ...),
        "bpcells" = bpcMatrixStore(path = path, ...),
        "tiledb" = tileDBArrayStore(path = path, ...)
    )
    store
}




# initialize ####

# cached fields check
# used for repeated checks across cascading initialize
.get_fields <- function(x) {
    if (length(x@fields) > 0L) return(x@fields)
    colnames(x)
}

# test if a store has been written and 'exists'
.store_exists <- function(x) {
    # no path provided
    if (length(x@path) == 0L || is.na(x@path) || x@path == "") {
        return(FALSE)
    }

    # test if path is remote URI
    if (grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", x@path)) {
        # Expensive: actually try to open it
        if (is.null(x@read_fun)) return(TRUE)
        remote_status <- tryCatch({
            dataset <- x@read_fun(x@path)
            !is.null(dataset)
        }, error = function(e) {
            FALSE  # Could be auth error OR missing file
        })
        return(remote_status)
    }

    file.exists(x@path)
}

setMethod("initialize", signature("parquetStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    # register default read_fun
    .Object@read_fun <- function(x) arrow::open_dataset(sources = x)
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
