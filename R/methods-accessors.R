# internal union class for dispatch
setClassUnion(".index", c("numeric", "logical", "character"))

# definitions ####

# * [ ####
setMethod("[", signature("parquetBase", ".index", ".index", "missing"), 
    function(x, i, j, ..., drop) {
    .guard_pstore_i_index(i)
})

setMethod("[", signature("parquetBase", ".index", "missing", "missing"), 
    function(x, i, j, ..., drop) {
    .guard_pstore_i_index(i)
})

setMethod("[", signature("parquetBase", "missing", "character", "missing"),
    function(x, i, j, ..., drop) {
    # Universe is the effective schema (on-disk + cols brought in by queued
    # join ops), not just disk_fields -- otherwise narrowing to a y-side col
    # after a join would error here. The final post-ops select in
    # `.pbase_storeread_processing` enforces this narrowing on the
    # materialized result.
    sdiff <- setdiff(j, .pstore_effective_schema(x))
    if (length(sdiff) > 0L) {
        stop(sprintf("[parquetStore] cols %s do not exist on disk",
            toString(sdiff)), call. = FALSE)
    }
    x@fields <- j
    x
})

setMethod("[", signature("parquetBase", "missing", "numeric", "missing"),
    function(x, i, j, ..., drop) {
    checkmate::assert_integerish(j)
    fields <- colnames(x)
    too_many <- any(j > length(fields))
    if (too_many) stop("[parquetStore] not that many columns")
    get_fields <- fields[j]
    x[, j = get_fields, ...]
})

setMethod("[", signature("parquetBase", "missing", "logical", "missing"),
    function(x, i, j, ..., drop) {
    fields <- colnames(x)
    len_f <- length(fields)
    if (length(j) < len_f) j <- rep(j, length.out = len_f)
    get_fields <- fields[j]
    x[, j = get_fields, ...]
})

setMethod("[", signature("parquetBase", "parquetBase", "missing"),
    function(x, i, j, ..., drop) {
    dots <- list(...)
    on <- dots$on
    # data.table convention: default `nomatch = NA` preserves x rows with NA
    # fill for unmatched y; explicit `nomatch = NULL` drops unmatched rows
    # (inner). NULL is lost when stored in a list -- recode as "inner"/"left"
    # for the queued op.
    nomatch <- if ("nomatch" %in% names(dots)) {
        nm <- dots$nomatch
        if (is.null(nm)) {
            "inner"
        } else if (length(nm) == 1L && is.na(nm)) {
            "left"
        } else {
            stop("[parquetStore] `nomatch` must be `NULL` (inner) or `NA` (left)",
                call. = FALSE)
        }
    } else {
        "left"  # data.table default: preserve x with NA fill on miss
    }
    checkmate::assert_character(on, min.len = 1L,
        .var.name = "on", null.ok = FALSE)
    x@ops <- c(x@ops, list(list(type = "join", y = i, by = on, nomatch = nomatch)))
    x
})

## internals ####
.guard_pstore_i_index <- function(i) {
    fmt <- "[parquetStore] `i` %s row indexing not supported"
    if (is.logical(i)) {
        fmt <- paste0(fmt, "\n If trying to downsample, use `rowSample()` instead.")
    }
    stop(sprintf(fmt, class(i)), call. = FALSE)
}

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

setMethod("nrow", signature("parquetBase"), function(x) {
    .dplyr_nrow(storeRead(x, output = "query"))
})

setMethod("ncol", signature("parquetBase"), function(x) {
    length(colnames(x))
})

setMethod("dim", signature("parquetBase"), function(x) {
    c(nrow(x), ncol(x))
})

# * ncol ####
setMethod("ncol", signature("fileStore"), function(x) {
    dims <- dim(x)
    if (is.null(dims)) return(0L)
    dims[2]
})

# * colnames ####
#' @export
setMethod("colnames", signature("parquetBase"), function(x) {
    # Effective schema: on-disk cols + cols recursively brought in by
    # queued join ops. Lets `colnames()` reflect what will be available at
    # the point downstream code is about to read or add another op.
    fields <- .pstore_effective_schema(x)
    if (!is.null(x@fields)) fields <- intersect(fields, x@fields)
    fields <- setdiff(fields, specialCols(x))
    fields
})

#' @export
setMethod("colnames", signature("unionParquetStore"), function(x) {
    fields <- lapply(x@stores, colnames) |>
        unlist() |>
        unique()
    if (!is.null(x@fields)) fields <- intersect(fields, x@fields)
    fields <- setdiff(fields, specialCols(x))
    fields
})

# * ext ####
#' @name ext
#' @title Spatial Extent
#' @aliases bbox
#' @description
#' Spatially mapped bounds.
#' `exact = TRUE` (default): scans coordinates with all spatial filter ops
#' applied; pending transform is projected during the scan.
#' `exact = FALSE`: fast estimate from metadata bounds (`@crop` > `disk_extent`)
#' intersected with `@window` and projected through any pending transform.
#' Axis-aligned transforms are exact; rotation/shear gives a conservative AABB.
#' Row-level ops (`subset`, `head`, etc.) are never reflected in either mode.
#' @param x object to use
#' @param exact `logical(1)`. If `TRUE` (default), scans for a true extent.
#'   If `FALSE`, returns a fast metadata-based estimate without scanning.
#' @param ... additional params to pass (not used)
#' @export
setMethod("ext", signature("parquetGeomBase"), function(x, exact = TRUE, ...) {
    aff <- .pgeom_pending_transform(x)
    if (!exact) return(.pgeom_ext_estimate(x, aff))
    q <- storeRead(x, output = "query")
    if (is.null(aff)) .dplyr_ext(q) else .dplyr_ext_affine(q, aff)
})

# Always returns extent in intrinsic (on-disk) x_index/y_index space,
# regardless of any pending "transform" op.  Used internally by crop(),
# window<-, and affine() where intrinsic bounds are required.
.pgeom_ext_intrinsic <- function(x) {
    .dplyr_ext(storeRead(x, output = "query"))
}

setMethod("geomtype", signature("parquetGeomBase"), function(x) {
    x@geomtype
})
