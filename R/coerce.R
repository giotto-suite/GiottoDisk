setAs("parquetStore", "parquetTileStore", \(from) {
    attr(from, tile) <- new("tileIterator")
    attr(from, "class") <- "parquetTileStore"

    initialize(from)
})

# R geospatial conversions ####

setMethod("as.sf", "parquetGeomStore", \(x, ...) {

})

setMethod("as.terra", "parquetGeomStore", \(x, ...) {

})

setMethod("as.data.frame", "parquetStore", \(x, ...) {

})

