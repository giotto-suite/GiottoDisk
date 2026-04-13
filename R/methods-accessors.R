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
    sdiff <- setdiff(j, .pstore_disk_fields(x))
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
    # nomatch = NULL (inner join, data.table convention) is indistinguishable
    # from "not passed" via dots$nomatch — check names explicitly.
    # Store as "inner"/"left" string to survive list storage (NULL is lost).
    nomatch <- if ("nomatch" %in% names(dots)) {
        if (is.null(dots$nomatch)) "inner" else "left"
    } else {
        "inner"  # default: inner join — unmatched rows indicate a data problem
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
    fields <- .pstore_disk_fields(x) # get on-disk state
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
#' Spatially mapped bounds
#' @param x object to use
#' @param ... additional params to pass (not used)
#' @export
setMethod("ext", signature("parquetGeomBase"), function(x, ...) {
    .dplyr_ext(storeRead(x, output = "query"))
})

setMethod("geomtype", signature("parquetGeomBase"), function(x) {
    x@geomtype
})
