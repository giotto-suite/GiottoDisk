#' @include class-dataStore.R
#' @include class-parquetStore.R
#' @include class-parquetExprStore.R
#' @include class-parquetEdgeStore.R
NULL

#' @name fingerprint_data-stores
#' @title Fingerprint a store from its handle
#' @description
#' [GiottoClass::fingerprint_data()] methods for disk-backed stores, used by
#' [GiottoClass::objManifest()]. A store is fingerprinted from its handle and
#' never read: `@uid` names the artifact on disk, and the remaining slots are
#' the view of it (queued ops, selected ids and indices, crop and window,
#' tile filter, the schema recorded at write time). Union stores fingerprint
#' their substores.
#'
#' `@path` is left out, so moving an artifact into a project vault keeps its
#' fingerprint, as are `@read_fun`, the cached marginals in `@stats` and the
#' tile plan, which is fixed when the artifact is written. A handle holding
#' a value with no canonical form, such as a function, errors, and the
#' manifest reports the store as unfingerprinted.
#' @param x a store
#' @param ... ignored
#' @returns character scalar
#' @keywords internal
NULL

# Slots that do not describe what the handle reads. Everything else is
# fingerprinted, so a new state slot is covered by default; a cache or
# location slot has to be listed here.
.FP_SKIP_SLOTS <- c("path", "read_fun", "stats", "tiles")

.fp_store <- function(x, ...) {
    sn <- setdiff(methods::slotNames(x), .FP_SKIP_SLOTS)
    parts <- lapply(sn, function(s) .fp_canon(methods::slot(x, s)))
    fingerprint_data(c(
        list(class = class(x)[[1L]]),
        stats::setNames(parts, sn)
    ))
}

# Reduce a handle value to atomic vectors, lists and payloads that
# fingerprint_data() knows, so GiottoClass does the hashing. Queued ops hold
# calls with inlined values, arrow id tables and other stores.
.fp_canon <- function(x) {
    if (is.null(x) || is.atomic(x)) return(x)
    if (is.name(x)) return(list(symbol = as.character(x)))
    if (is.call(x)) return(list(call = lapply(as.list(x), .fp_canon)))
    if (is.data.frame(x)) return(x)
    if (inherits(x, c("Table", "RecordBatch"))) return(as.data.frame(x))
    if (inherits(x, c("Array", "ChunkedArray"))) return(as.vector(x))
    if (is.list(x)) return(lapply(x, .fp_canon))
    if (isS4(x)) return(fingerprint_data(x))
    stop("[fingerprint_data] no canonical form for ", class(x)[[1L]],
        call. = FALSE)
}

#' @rdname fingerprint_data-stores
#' @importFrom GiottoClass fingerprint_data
#' @export
setMethod("fingerprint_data", signature("dataStore"), .fp_store)

# union parquet stores extend parquetBase but not dataStore
#' @rdname fingerprint_data-stores
#' @export
setMethod("fingerprint_data", signature("parquetBase"), .fp_store)
