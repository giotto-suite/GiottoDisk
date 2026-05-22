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
#' @slot uid `character` automatically generated unique ID for artifact
#'   tracking. See [artifact_uid]
#' @slot params `list` Additional params such as remote names that are 
#' relevant to the store type may be stored here as a named list.
#' @slot read_fun `function`. Function to access data where `read_fun(path)`
#'   returns the data in a useful format. Can be customized for
#'   compression or other access requirements. Do not embed credentials in
#'   `read_fun()`
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
fileStore <- function(path = .dump_tempfile(), read_fun, ...) {
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

#' @name arrayStore
#' @title Array Storage
#' @param path storage directory
#' @param name `character` Remote naming used by some storage types,
#'   for example HDF5Array and TileDB.
#' @export
h5ArrayStore <- function(
    path = .dump_tempfile(), name = HDF5Array::getHDF5DumpName(), ...) {
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
#' * `"parquetExpr"`
#' * `"parquetEdge"`
#' * `"h5"`
#' * `"bpcells"`
#' * `"tiledb"`
#' @family store constructors
#' @seealso [store]
#' @export
storeCreate <- function(path = .dump_tempfile(), type = "parquet", ...) {
  .store_type_aliases <- c(
      "parquet"         = "parquetstore",
      "parquetgeom"     = "parquetgeomstore",
      "parquetgeomtile" = "parquetgeomtilestore",
      "parquetexpr"     = "parquetexprstore",
      "parquetedge"     = "parquetedgestore",
      "file"            = "filestore",
      "h5"              = "h5arraystore",
      "bpcells"         = "bpcmatrixstore",
      "tiledb"          = "tiledbarraystore"
  )
  type <- tolower(type)
  type <- .store_type_aliases[type] %na% type

    type <- match.arg(type,
        choices = c(
            "parquetstore",
            "parquetgeomstore",
            "parquetgeomtilestore",
            "parquetexprstore",
            "parquetedgestore",
            "filestore",
            "h5arraystore",
            "bpcmatrixstore",
            "tiledbarraystore"
        )
    )
    store <- switch(type,
        "parquetstore" = parquetStore(path = path, ...),
        "parquetgeomstore" = parquetGeomStore(path = path, ...),
        "parquetgeomtilestore" = parquetGeomTileStore(path = path, ...),
        "parquetexprstore" = parquetExprStore(path = path, ...),
        "parquetedgestore" = parquetEdgeStore(path = path, ...),
        "filestore" = fileStore(path = path, ...),
        "h5arraystore" = h5ArrayStore(path = path, ...),
        "bpcmatrixstore" = bpcMatrixStore(path = path, ...),
        "tiledbarraystore" = tileDBArrayStore(path = path, ...)
    )
    store
}

# initialize ####

setMethod("initialize", signature("fileStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    if (length(.Object@uid) == 0L) .Object@uid <- .make_uid()
    .Object
})

# internals ####

# simple store read without additional processing
# useful for determining on-disk state.
.store_simple_read <- function(x, ...) {  
    x@read_fun(x@path, ...)
}
  