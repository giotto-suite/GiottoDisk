
setMethod("initialize", "tileMap", \(.Object, ...) {
    tiles <- .Object@tiles
    tiles <- initialize(tiles)

    if (length(tiles) == 0) return(.Object) # return early if no tiles
    #
})

setMethod("initialize", "ParquetSpatVector", function(.Object, ...) {
    .Object <- structure(terra::vect(), class = c("ParquetSpatVector", "SpatVector", "S4"))
    a <- list(...)
    attr(.Object, "file") <- if (is.null(a$file)) NA_character_
    else a$file
    .Object
})



