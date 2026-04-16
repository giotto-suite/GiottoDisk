#' @name storeExists
#' @title Store Existence
#' @description
#' Test if a store has been written and is ready to be used.
#' @param x `store` object
#' @param ... additional params to pass
#' @returns `logical`
NULL

#' @rdname storeExists
#' @export
setMethod("storeExists", signature("fileStore"), function(x, ...) {
    if (.test_empty_path(x@path)) return(FALSE)
  
    # determine if path is remote URI
    is_uri <- isTRUE(grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", x@path))
    if (is_uri) return(.test_store_exists_remote(x, x@path))
  
    # test local existence
    file.exists(x@path)
})

#' @rdname storeExists
#' @param all `logical` return single `all()` value instead of a vector of `logical`
#' for each substore
#' @export
setMethod("storeExists", signature("unionParquetStore"), function(x, all = TRUE, ...) {
    res <- vapply(x@stores, storeExists, FUN.VALUE = logical(1L))
    if (all) res <- all(res)
    res
})

# storePaths ####

#' @name storePaths
#' @title Store Paths
#' @description
#' Return the artifact-level path(s) for a store. For composite stores
#' (e.g. union stores), one path is returned per substore.
#' @param x `store` object
#' @param ... additional params to pass
#' @returns `character` vector of paths
#' @export
setMethod("storePaths", signature("fileStore"), function(x, ...) {
    x@path
})

#' @rdname storePaths
#' @export
setMethod("storePaths", signature("unionParquetStore"), function(x, ...) {
    vapply(x@stores, function(s) s@path, FUN.VALUE = character(1L))
})

# storeUID ####

#' @name storeUID
#' @title Store UIDs
#' @description
#' Return the UID(s) for a store. For union stores, one UID is returned per
#' substore. For `overlapPointDisk`, returns a named list with `poly` and
#' `feat` UID vectors capturing the provenance of the overlap computation.
#' @param x `store` or `overlapPointDisk` object
#' @param ... unused
#' @returns `character` vector of UIDs, or named `list` for `overlapPointDisk`
#' @export
setMethod("storeUID", signature("fileStore"), function(x, ...) {
    x@uid
})

#' @rdname storeUID
#' @export
setMethod("storeUID", signature("unionParquetStore"), function(x, ...) {
    vapply(x@stores, function(s) s@uid, FUN.VALUE = character(1L))
})

#' @rdname storeUID
#' @export
setMethod("storeUID", signature("overlapPointDisk"), function(x, ...) {
    list(poly = x@poly_uids, feat = x@feat_uids)
})

# internals ####
.test_empty_path <- function(path) {
    length(path) == 0L || is.na(path) || path == ""
}

.test_store_exists_remote <- function(x, uri) {
    if (is.null(x@read_fun)) return(TRUE)
    # when option FALSE: assume valid; let read-time errors surface
    test_remote <- getOption("giottodisk.uri_test_exists", FALSE)
    if (!isTRUE(test_remote)) return(TRUE)

    # Expensive: actually try to open it
    remote_status <- tryCatch({
        dataset <- x@read_fun(uri)
        !is.null(dataset)
    }, error = function(e) {
        FALSE  # Could be auth error OR missing file
    })
    return(remote_status)
}
