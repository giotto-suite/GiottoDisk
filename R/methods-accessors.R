
# definitions ####

# * nrow ####
setMethod("nrow", signature("parquetStore"), function(x) {
    atab <- storeRead(x)
    .dplyr_nrow(atab)
})

# * colnames ####
setMethod("colnames", signature("parquetStore"), function(x) {
    atab <- storeRead(x)
    colnames(atab)
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
