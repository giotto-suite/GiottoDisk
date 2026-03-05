
# definitions ####

# * dim ####
setMethod("dim", signature("fileStore"), function(x) {
    dim(storeRead(x))
})

# * nrow ####
setMethod("nrow", signature("fileStore"), function(x) {
    dims <- dim(x)
    if (is.null(dims)) return(0L)
    dims[1]
})

setMethod("nrow", signature("parquetStore"), function(x) {
    atab <- storeRead(x)
    .dplyr_nrow(atab)
})

# * ncol ####
setMethod("ncol", signature("fileStore"), function(x) {
    dims <- dim(x)
    if (is.null(dims)) return(0L)
    dims[2]
})

# * colnames ####
setMethod("colnames", signature("parquetStore"), function(x) {
    atab <- storeRead(x)
    names(atab)
})

# * crop ####
setMethod("crop", signature("parquetGeomStore", "ANY"), function(x, y, ...) {
    e <- ext(x)
    e <- terra::intersect(e, ext(y))
    if (is.null(e)) stop("[crop] no geometries within requested extent\n", call. = FALSE)
    x@extent <- .ext_to_num_vec(e)
    x
})

# * ext ####
setMethod("ext", signature("parquetGeomStore"), function(x, ...) {
    if (length(x@extent) == 0L) {
        stop("no extent set\n", call. = FALSE)
    }
    ext(x@extent)
})

setMethod("geomtype", signature("parquetGeomStore"), function(x) {
    x@geomtype
})

# * rbind ####
setMethod("rbind2", signature("parquetStore", "parquetStore"), function(x, y, ...) {
    x@uid <- c(x@uid, y@uid)
    x@path <- c(x@path, y@path)
    if (anyDuplicated(x@path)) {                                 
        warning("[rbind] duplicate store paths detected — same data will be interleaved")             
  }
    x
})

setMethod("rbind2", signature("parquetStore", "parquetGeomStore"), function(x, y, ...) {
    rbind2(y, x, ...)
})

setMethod("rbind2", signature("parquetGeomStore", "parquetStore"), function(x, y, ...) {
    if (!inherits(y, "parquetGeomStore")) {
        stop("[rbind] parquetGeomStore can only be rbinded to other parquetGeomStore-inheriting")
    }
    x <- callNextMethod(x, y, ...)
    x@extent <- .ext_to_num_vec(terra::union(ext(x@extent), ext(y@extent)))
    x
})
