setMethod("show", signature("fileStore"), function(object) {
    cat(sprintf("<%s>\n", class(object)))
    cat(sprintf("path: %s\n", 
        paste0(str_abbreviate(object@path), collapse = "\n      ")
    ))
    if (!.test_store_written(object@path)) {
        cat("<empty>")
        return(invisible())
    }
})

setMethod("show", signature("parquetStore"), function(object) {
    callNextMethod(object)
    if (!.test_store_written(object@path)) return(invisible())
    cat(sprintf("columns: %s\n", paste(object@fields, collapse = ", ")))
    cat(sprintf("nrows: %d\n", as.integer(nrow(object))))
})

setMethod("show", signature("parquetGeomStore"), function(object) {
    callNextMethod(object)
    if (!.test_store_written(object@path)) return(invisible())
    if (length(object@extent) > 0L) {
        cat(sprintf("extent: %s (xmin, xmax, ymin, ymax)\n",
                    paste(collapse = ", ", round(object@extent, digits = 3))))
    } else {
        cat("extent: none settable: check 'x_index' and 'y_index' are columns\n")
    }
})

setMethod("show", signature("parquetGeomTileStore"), function(object) {
    callNextMethod(object)
    if (!.test_store_written(object@path)) return(invisible())
    cat(sprintf("tiles: %d\n", length(object@tiles)))
})

setMethod("show", signature("h5ArrayStore"), function(object) {
    callNextMethod(object)
    if (!.test_store_written(object@path)) return(invisible())
    cat(sprintf("name: \"%s\"\n", object@params$name))
})

setMethod("show", signature("tileDBArrayStore"), function(object) {
    callNextMethod(object)
    if (!.test_store_written(object@path)) return(invisible())
    cat(sprintf("name: \"%s\"\n", object@params$name))
})
