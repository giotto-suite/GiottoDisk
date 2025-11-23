#' @include utils.R
#' @include pkg_imports.R
NULL

# docs ####

#' @name store
#' @title Data Storage
#' @description
#' S4 class define a storage of data. The base class `dataStore` is VIRTUAL
#' @slot name description
#' @family store types
NULL

#' @name fileStore-class
#' @title File Store
#' @description
#' S4 class defining a file location and read spec. `fileStore` is a
#' general purpose customizable class that can be used for non-standard file
#' formats to read them in with [storeRead()]. A function to read the file
#' just has to be supplied to the `read_fun` param of the `fileStore()`
#' constructor.
#'
#' Specific file format stores can extend this class.
#' @param path disk path
#' @param read_fun `function`. Reader function for the file format. The first
#' param should be the `path` given as input, such that `read_fun(path)`
#' reads in the data in a useful format.
#' @family store_types
NULL

#' @name parquetStore-class
#' @title Parquet Store
#' @description
#' S4 Class extending [fileStore] with specific expected columns for
#' indexed storage of tabular data. `parquetStores` are not intended to be
#' updated or edited. They only query and do delayed ops. They can be generated
#' using constructors described in [parquetStore].
#'
#' No table is contained so that the idea of this class storing state cannot be
#' entertained.
#' @slot path `character`. Disk path to file or hive storage directory
#' @slot fields `character`. Cached column names.
#' @slot read_fun `function`. Reader function for the file format.
#' @slot extent `numeric(4)`. xmin, xmax, ymin, ymax extent of data contained.
#' @slot tiles a \{GiottoTile\} `tileIterator`. Provides understanding of
#' which tile(s) in a tile plan this store instance is responsible for.
#' @family store types
NULL







# definitions ####

setOldClass("matrix")
setOldClass("data.frame")
setOldClass("data.table")
setClassUnion("memoryMatrix", c("matrix", "Matrix"))
setClassUnion("memoryStore", c("memoryMatrix", "data.frame", "data.table"))

#' @rdname store
setClass("dataStore", contains = "VIRTUAL")
#' @rdname store
setClass("fileStore",
    contains = "dataStore",
    slots = list(
        path = "character",
        read_fun = "function"
    )
)

# * parquet ####

# cols: row_index, id, ...
#' @rdname parquetStore
setClass("parquetStore",
    contains = "fileStore",
    slots = list(
        path = "character",
        fields = "character"
    )
)

# cols: row_index, x_index, y_index, geom, id, ...
#' @rdname parquetStore
setClass("parquetGeomStore",
    contains = "parquetStore",
    slots = list(extent = "numeric")
)

# cols: row_index, x_index, y_index, tile_index, geom, id, ...
#' @rdname parquetStore
setClass("parquetGeomTileStore",
    contains = "parquetGeomStore",
    slots = list(tiles = "tileIterator")
)


# * delayed ####

setClass("h5ArrayStore",
    contains = "fileStore",
    slots = list(
        name = "character"
    )
)

setClass("tileDBArrayStore",
    contains = "fileStore",
    slots = list(
        name = "character"
    )
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
#' @family Store constructors
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
#' @family Store constructors
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
#' @export
h5ArrayStore <- function(
    path = tempfile(), name = HDF5Array::getHDF5DumpName(), ...) {
    package_check("HDF5Array", repository = "Bioc")
    new("h5ArrayStore", path = path, name = name, ...)
}

#' @rdname arrayStore
#' @export
tileDBArrayStore <- function(
        path = file.path(tempdir(), .random_id()),
        name = TileDBArray::getTileDBAttr(),
        ...
    ) {
    package_check("TileDBArray", repository = "Bioc")
    new("tileDBArrayStore", path = path, name = name, ...)
}


#' @name storeCreate
#' @title Create a Store
#' @description
#' Hub function for creating a concrete store class.
#' @param path `character`. Disk path to file or hive storage directory
#' @param type `character`. Type of store to create. Currently one of
#' `"parquet"`, `"parquetGeom"`, `"parquetGeomTile"`, or `"file"`.
#' @family Store constructors
#' @seealso [store]
#' @export
storeCreate <- function(path = tempfile(), type = "parquet", ...) {
    type <- match.arg(type,
        choices = c("parquet", "parquetGeom", "parquetGeomTile", "file")
    )
    store <- switch(type,
        "parquet" = parquetStore(path = path, ...),
        "parquetGeom" = parquetGeomStore(path = path, ...),
        "parquetGeomTile" = parquetGeomTileStore(path = path, ...),
        "file" = fileStore(path = path, ...)
    )
    store
}




# initialize ####

setMethod("initialize", signature("parquetStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    read_fun <- function(x) arrow::open_dataset(sources = x)
    .Object@read_fun <- read_fun
    exists <- checkmate::test_file_exists(.Object@path)
    if (!exists) return(.Object)
    .Object@fields <- colnames(.Object)
    if (!"row_index" %in% .Object@fields) {
        warning("[initialize] 'row_index' column should be present",
                call. = FALSE)
    }
    .Object
})

setMethod("initialize", signature("parquetGeomStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    exists <- checkmate::test_file_exists(.Object@path)
    if (!exists) return(.Object)
    .Object@fields <- colnames(.Object)
    check_cols <- c("x_index", "y_index", "geom")
    missing_cols <- check_cols[!check_cols %in% .Object@fields]
    if (length(missing_cols) > 0L) {
        warning(sprintf("[initialize] '%s' column(s) should be present",
                        paste(collapse = "', '", missing_cols)),
                call. = FALSE)
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
