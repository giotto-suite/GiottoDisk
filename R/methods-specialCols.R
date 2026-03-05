#' @name specialCols
#' @title Special Store Columns
#' @description
#' Documents columns that are special to particular store implementations.
#' These columns are enforced on disk and available when accessing the
#' data as a query. These columns may be automatically consumed depending on
#' the materialization format.
#' @param store `dataStore`-inheriting class
#' @returns `character`
NULL

#' @rdname specialCols
#' @export
setMethod("specialCols", signature("dataStore"), function(store) character(0L))
#' @rdname specialCols
#' @export
setMethod("specialCols", signature("parquetStore"), function(store) {
    c(callNextMethod(store), "row_index", "source_id")
})
#' @rdname specialCols
#' @export
setMethod("specialCols", signature("parquetGeomStore"), function(store) {
    c(callNextMethod(store), "x_index", "y_index", "geom")
})
#' @rdname specialCols
#' @export
setMethod("specialCols", signature("parquetGeomTileStore"), function(store) {
    c(callNextMethod(store), "tile_index")
})
