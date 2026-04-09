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
setMethod("specialCols", signature("ANY"), function(store) character(0L))
#' @rdname specialCols
#' @export
setMethod("specialCols", signature("parquetBase"), function(store) {
    c(callNextMethod(store), "row_index", "source_id")
})
#' @rdname specialCols
#' @export
setMethod("specialCols", signature("parquetGeomBase"), function(store) {
    c(callNextMethod(store), "x_index", "y_index", "geom", "tile_index")
})
#' @rdname specialCols
#' @export
setMethod("specialCols", signature("unionParquetStore"), function(store) {
    unique(unlist(lapply(store@stores, specialCols)))
})
