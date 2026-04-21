# Spatial manipulation verbs for parquetGeomBase
#
# affine(parquetGeomBase, affine2d) is the primitive: it sets the anchor,
# composes with any existing pending transform, and writes a single
# "transform" entry into @post_ops.
#
# Each higher-level verb retrieves the current pending affine2d (or creates
# a fresh identity), applies the corresponding operation to the affine2d,
# and writes the result back via .pgeom_set_transform() directly -- bypassing
# the anchor / composition logic in affine() that is only needed when the
# caller supplies an affine2d with an unknown anchor.
#
# Verbs that need a centroid (spin, rescale, shear) live-scan the intrinsic
# extent to anchor the affine2d before delegating to the affine2d method.
# Verbs with no centroid dependency (spatShift, flip) skip the scan.



# helpers ####

# return the pending affine2d from @post_ops, or NULL
.pgeom_pending_transform <- function(x) {
    ops <- Filter(function(op) op$type == "transform", x@post_ops)
    if (length(ops) == 0L) return(NULL)
    ops[[1L]]$affine2d
}

# retrieve current pending affine2d or initialise a fresh identity
.pgeom_get_transform <- function(x) {
    aff <- .pgeom_pending_transform(x)
    if (is.null(aff)) GiottoClass::affine()  # identity affine2d
    else aff
}

# write an affine2d back as the sole "transform" entry in @post_ops
.pgeom_set_transform <- function(x, aff) {
    x@post_ops <- c(
        Filter(function(op) op$type != "transform", x@post_ops),
        list(list(type = "transform", affine2d = aff))
    )
    x
}



# affine ####

#' @name affine
#' @title Lazily record an affine transform on a parquetGeomBase store
#' @description
#' Record an affine transformation on a `parquetGeomBase`-inheriting store. The
#' transform is stored as a `"transform"` entry in `@post_ops` and applied to output
#' coordinates at materialization time (`storeRead(output = "tibble"/"terra"/"sf")`).
#' It is a no-op for the Arrow query phase.
#'
#' Subsequent `affine()` calls compose with any existing pending transform -- the
#' chain collapses to a single `affine2d` entry. Spatial filters
#' ([crop()]/[`window<-`]) applied after an `affine()` call back-project their
#' query extents to intrinsic (on-disk) space for accurate Arrow filtering.
#'
#' `output = "query"` and `output = "duckdb"` are unaffected.
#' @param x `parquetGeomBase`-inheriting store
#' @param y `affine2d` transform object (from GiottoClass)
#' @param ... additional arguments (ignored)
#' @returns `x` with the transform recorded in `@post_ops`
#' @export
setMethod("affine", signature("parquetGeomBase", "affine2d"),
    function(x, y, ...) {
    # Set anchor only if y still carries the prototype default.
    # Callers that need an accurate anchor (spin, rescale) set the extent
    # on their affine2d externally before passing it here.
    if (isTRUE(all.equal(y@anchor, c(-180, 180, -90, 90)))) {
        ext(y) <- .pgeom_ext_intrinsic(x)
    }

    existing <- .pgeom_pending_transform(x)
    if (!is.null(existing)) {
        # Compose: apply y on top of existing. affine(affine2d, matrix) uses
        # the 3x3 homogeneous matrix (including translation) for composition.
        y <- affine(existing, y@affine)
    }

    .pgeom_set_transform(x, y)
})



# spin ####

#' @name spin
#' @title Rotate a parquetGeomBase store
#' @description
#' Lazily record a rotation on a `parquetGeomBase`-inheriting store. The
#' rotation is composed with any existing pending transform and applied at
#' read time. The centroid anchor is live-scanned so that any preceding
#' crops or filters are reflected in the rotation pivot.
#' @param x `parquetGeomBase`-inheriting store
#' @param angle rotation angle in degrees
#' @param x0,y0 pivot coordinates. If `NULL` (default), the centroid of
#'   the current data extent is used.
#' @param ... additional arguments (ignored)
#' @returns `x` with the rotation recorded in `@ops`
#' @export
setMethod("spin", signature("parquetGeomBase"), function(x, angle, x0 = NULL, y0 = NULL, ...) {
    aff <- .pgeom_get_transform(x)
    ext(aff) <- .pgeom_ext_intrinsic(x)  # live-scan for accurate centroid
    a <- list(x = aff, angle = angle)
    if (!is.null(x0)) a$x0 <- x0
    if (!is.null(y0)) a$y0 <- y0
    aff <- do.call(spin, args = a)
    .pgeom_set_transform(x, aff)
})



# rescale ####

#' @name rescale
#' @title Scale a parquetGeomBase store
#' @description
#' Lazily record a scaling transform on a `parquetGeomBase`-inheriting store.
#' @param x `parquetGeomBase`-inheriting store
#' @param fx,fy scale factors for x and y axes. `fy` defaults to `fx`.
#' @param x0,y0 pivot coordinates. If omitted, the centroid of the current
#'   data extent is used.
#' @param ... additional arguments (ignored)
#' @returns `x` with the scaling recorded in `@ops`
#' @export
setMethod("rescale", signature("parquetGeomBase"), function(x, fx = 1, fy = fx, x0, y0, ...) {
    aff <- .pgeom_get_transform(x)
    ext(aff) <- .pgeom_ext_intrinsic(x)
    a <- list(x = aff, fx = fx, fy = fy)
    if (!missing(x0)) a$x0 <- x0
    if (!missing(y0)) a$y0 <- y0
    aff <- do.call(rescale, args = a)
    .pgeom_set_transform(x, aff)
})



# shear ####

#' @name shear
#' @title Shear a parquetGeomBase store
#' @description
#' Lazily record a shear transform on a `parquetGeomBase`-inheriting store.
#' @param x `parquetGeomBase`-inheriting store
#' @param fx,fy shear factors for x and y axes (default 0 = no shear)
#' @param x0,y0 shear centre coordinates. If omitted, the centroid of the
#'   current data extent is used.
#' @param ... additional arguments (ignored)
#' @returns `x` with the shear recorded in `@ops`
#' @export
setMethod("shear", signature("parquetGeomBase"), function(x, fx = 0, fy = 0, x0, y0, ...) {
    aff <- .pgeom_get_transform(x)
    ext(aff) <- .pgeom_ext_intrinsic(x)
    a <- list(x = aff, fx = fx, fy = fy)
    if (!missing(x0)) a$x0 <- x0
    if (!missing(y0)) a$y0 <- y0
    aff <- do.call(shear, args = a)
    .pgeom_set_transform(x, aff)
})



# spatShift ####

#' @name spatShift
#' @title Translate a parquetGeomBase store
#' @description
#' Lazily record a translation on a `parquetGeomBase`-inheriting store.
#' @param x `parquetGeomBase`-inheriting store
#' @param dx,dy translation distances along x and y axes (default 0)
#' @param ... additional arguments (ignored)
#' @returns `x` with the translation recorded in `@ops`
#' @export
setMethod("spatShift", signature("parquetGeomBase"), function(x, dx = 0, dy = 0, ...) {
    aff <- .pgeom_get_transform(x)
    aff <- spatShift(aff, dx = dx, dy = dy)
    .pgeom_set_transform(x, aff)
})



# t ####

#' @name t
#' @title Transpose a parquetGeomBase store
#' @description
#' Lazily record a transpose (swap x and y coordinates) on a
#' `parquetGeomBase`-inheriting store.
#' @param x `parquetGeomBase`-inheriting store
#' @returns `x` with the transpose recorded in `@post_ops`
#' @export
setMethod("t", signature("parquetGeomBase"), function(x) {
    aff <- .pgeom_get_transform(x)
    aff <- t(aff)
    .pgeom_set_transform(x, aff)
})



# flip ####

#' @name flip
#' @title Flip a parquetGeomBase store
#' @description
#' Lazily record a reflection on a `parquetGeomBase`-inheriting store.
#' @param x `parquetGeomBase`-inheriting store
#' @param direction `"vertical"` (default, reflects over a horizontal axis)
#'   or `"horizontal"` (reflects over a vertical axis)
#' @param x0,y0 coordinate of the reflection axis (default 0). For a vertical
#'   flip, `y0` is the y-position of the axis; for horizontal, `x0`.
#' @param ... additional arguments (ignored)
#' @returns `x` with the reflection recorded in `@ops`
#' @export
setMethod("flip", signature("parquetGeomBase"), function(x, direction = "vertical", x0 = 0, y0 = 0, ...) {
    aff <- .pgeom_get_transform(x)
    aff <- flip(aff, direction = direction, x0 = x0, y0 = y0)
    .pgeom_set_transform(x, aff)
})
