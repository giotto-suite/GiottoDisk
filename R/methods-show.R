setMethod("show", signature("gDirSource"), function(object) {
    cat(sprintf("<%s>\n", class(object)))

    if (file.exists(object@path)) {
        object <- object@read()
        cat("stores:\n")
        print_list(object@catalog$stores)
        cat("\n")
        cat("artifacts:", length(object@catalog$artifacts), "\n")
        cat("versions:", length(object@catalog$versions), "\n")
    } else {
        cat("giottodir.json not written yet. Use `@write()`\n")
    }

    cat("\n")
    cat("* `@read(path)` to read from a Giotto directory json\n")
    cat("* `@write(path)` to write a Giotto directory json\n")
})

setMethod("show", signature("fileStore"), function(object) {
    cat(sprintf("<%s>\n", class(object)))
    cat(sprintf("path: %s\n", str_abbreviate(object@path)))
    if (!file.exists(object@path)) {
        cat("<empty>")
        return(invisible())
    }
})

setMethod("show", signature("parquetStore"), function(object) {
    callNextMethod(object)
    if (!file.exists(object@path)) return(invisible())
    cat(sprintf("columns: %s\n", paste(object@fields, collapse = ", ")))
    cat(sprintf("nrows: %d\n", as.integer(nrow(object))))
})

setMethod("show", signature("parquetGeomStore"), function(object) {
    callNextMethod(object)
    if (!file.exists(object@path)) return(invisible())
    if (length(object@extent) > 0L) {
        cat(sprintf("extent: %s (xmin, xmax, ymin, ymax)\n",
                    paste(collapse = ", ", round(object@extent, digits = 3))))
    } else {
        cat("extent: none settable: check 'x_index' and 'y_index' are columns\n")
    }
})

setMethod("show", signature("parquetGeomTileStore"), function(object) {
    callNextMethod(object)
    if (!file.exists(object@path)) return(invisible())
    cat(sprintf("tiles: %d\n", length(object@tiles)))
})

setMethod("show", signature("h5ArrayStore"), function(object) {
    callNextMethod(object)
    if (!file.exists(object@path)) return(invisible())
    cat(sprintf("name: \"%s\"\n", object@name))
})

setMethod("show", signature("tileDBArrayStore"), function(object) {
    callNextMethod(object)
    if (!file.exists(object@path)) return(invisible())
    cat(sprintf("name: \"%s\"\n", object@name))
})
