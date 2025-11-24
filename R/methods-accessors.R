
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

# * $ ####

setMethod("$", "gDirSource", function(x, name) {
    file.path(x@path, x@catalog$stores[[name]])
})

#' @keywords internal
#' @export
.DollarNames.gDirSource <- function(x, pattern) names(x@catalog$stores)

# * [ ####

#' @keywords internal
#' @export
setMethod("[", c(x = "gDirSource", i = "character", j = "missing", drop = "missing"), function(x, i, j, ..., drop) {
    file.path(x@path, x@catalog$stores[[i]])
})
