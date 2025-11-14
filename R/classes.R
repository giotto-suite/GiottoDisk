#' @include pkg_imports.R
NULL

setClass("dataStore", contains = "VIRTUAL")
setClass("memoryStore", contains = "dataStore")

setClass("fileStore",
    contains = "dataStore",
    slots = list(
        path = "character",
        read_fun = "function"
    )
)

# cols: row_index, id, ...
setClass("parquetStore",
    contains = "fileStore",
    slots = list(
        path = "character",
        fields = "character"
    )
)

# cols: row_index, x_index, y_index, geom, id, ...
setClass("parquetGeomStore", contains = "parquetStore", slots = list(extent = "numeric"))

# cols: row_index, x_index, y_index, tile_index, geom, id, ...
setClass("parquetGeomTileStore",
         contains = "parquetGeomStore", slots = list(tiles = "tileIterator"))

diskGiotto <- setClass(
    Class = "diskGiotto",
    contains = "giotto",
    slots = list(
        src = "character"
    ),
    prototype = list(
        src = NA_character_
    )
)

# flag that it's a proxy subobject?
# setClass(Class = "proxy", contains = "VIRTUAL", slots = list(
#
# ))

#'
#' tileMap <- setClass(
#'     "tileMap",
#'     slots = list(
#'         tiles = "tilePlan",
#'         affine = "affine2d"
#'     )
#' )
#'
#'
#'
#'
#' nullfun <- \() NULL
#'
#' #' @name parquetStore
#' #' @title Parquet Data Storage
#' #' @description
#' #' Classes for representing parquet storage.
#' #' @slot schema store the call needed to (re)generate an {arrow} schema for the
#' #' data it references. See `?arrow::Schema`. Generate the call using the active
#' #' binding on the Schema R6 object: `$code(namespace = TRUE)`
#' setClass(
#'     "parquetStore",
#'     contains = "dataStore",
#'     slots = list(
#'         schema = "call"
#'     ),
#'     prototype = list(
#'         schema = call("nullfun")
#'     )
#' )
#'
#' #' @rdname parquetStore
#' setClass(
#'     "parquetGeomStore",
#'     contains = "parquetStore"
#' )
#'
#' setClass(
#'     "parquetGeomTileStore",
#'     contains = "parquetGeomStore",
#'     slots = list(
#'         tiles = "tileIterator"
#'     )
#' )
#'
#'
#' # processing parameters ####
#' #' @name dataBatch
#' #' @slot store dataStore object to read in values from
#' #' @slot file output file or directory (can be none if lazy)
#' setClass(
#'     "dataBatch",
#'     contains = "VIRTUAL",
#'     slots = list(
#'         store = "dataStore",
#'         file = "character"
#'     )
#' )
#'
#' # default batching by arrow FileSystemDataset interface
#' setClass(
#'     "nativeBatch",
#'     contains = "dataBatch",
#'     slots = list(
#'         idx = "numeric"
#'     )
#' )
#'
#' # selected (manual) rows per batch
#' setClass(
#'     "selectBatch",
#'     contains = "dataBatch",
#'     slots = list(
#'         idx = "numeric"
#'     )
#' )
#'
#' # spatially tiled batching
#' setClass(
#'     "tileBatch",
#'     contains = "dataBatch",
#'     slots = list(
#'         tiles = "tileIterator",
#'         idx = "numeric"
#'     )
#' )
#'
#'
#'
#'
#'
#' setClass(
#'     "computeResult",
#'     slots = list(
#'         type = "character",  # "lazy", "eager", "spatial"
#'         operations = "list",  # List of operations in order
#'         source = "ANY"       # Original data source
#'     ),
#'     prototype = list(
#'         type = "lazy",
#'         operations = list()
#'     )
#' )
#'
#' # Operation class to store each step
#' setClass(
#'     "spatialOperation",
#'     slots = list(
#'         fn = "function",
#'         params = "list",
#'         name = "character"
#'     )
#' )
#'
#'
#'
#' setClass("ParquetSpatVector", contains = "SpatVector", slots = c(
#'     file = "character"
#' ))




