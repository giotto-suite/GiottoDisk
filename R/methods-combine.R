#' @name rbind
#' @title rbind stores
#' @description
#' Combine stores by rbinding. Methods provided for parquetStore
#' which implement this behavior using union partquet dataset
#' access. Note that rbind for parquetStore only works on
#' pristine stores with no lazy ops attached. Please write out
#' the lazy steps first via [storeWrite()] if rbinding stores
#' with lazy ops.
#' @param x,y store objects to rbind
#' @param src (optional) `gsource` object to write to if rbinding
#'   the same store to itself. In these cases, a second copy is written
#'   to disk to make sure that all data is repesented on disk.
#'   Providing `src` writes to the `gsource`. Otherwise, it is written
#'   to the default dump location (see [artifact_dump])
#' @param ... additional params to pass
NULL 

#' @rdname rbind
#' @export
setMethod("rbind2", signature("parquetStore", "parquetStore"), function(x, y, src = NULL, ...) {
    .guard_rbind_geom_classes(x, y)
    .guard_pstore_lazy_ops(x, y)
    y <- .rbind_dedup(y, src = src, existing_uids = .pstore_all_uids(x))
    create_class <- .pstore_get_union_class(x)
    z <- new(create_class, stores = list(x, y))
    
    # inherit tracked params
    if (!is.null(x@fields) || !is.null(y@fields)) {
        z@fields <- unique(c(x@fields, y@fields))
    }
  
    initialize(z)
})

#' @rdname rbind
#' @export
setMethod("rbind2", signature("parquetStore", "unionParquetStore"), function(x, y, src = NULL, ...) {
    .guard_rbind_geom_classes(x, y)
    .guard_pstore_lazy_ops(x, y)
    x <- .rbind_dedup(x, src = src, existing_uids = .pstore_all_uids(y))
    
    # inherit tracked params
    if (!is.null(x@fields) || !is.null(y@fields)) {
        y@fields <- unique(c(x@fields, y@fields))
    }
  
    y@stores <- c(list(x), y@stores)
  
    initialize(y)
})

#' @rdname rbind
#' @export
setMethod("rbind2", signature("unionParquetStore", "parquetStore"), function(x, y, src = NULL, ...) {
    .guard_rbind_geom_classes(x, y)
    .guard_pstore_lazy_ops(x, y)
    y <- .rbind_dedup(y, src = src, existing_uids = .pstore_all_uids(x))

    # inherit tracked params
    if (!is.null(x@fields) || !is.null(y@fields)) {
        x@fields <- unique(c(x@fields, y@fields))
    }
    
    x@stores <- c(x@stores, list(y))
    
    initialize(x)
})

#' @rdname rbind
#' @export
setMethod("rbind2", signature("unionParquetStore", "unionParquetStore"), function(x, y, src = NULL, ...) {
    .guard_rbind_geom_classes(x, y)
    .guard_pstore_lazy_ops(x, y)
    y <- .rbind_dedup(y, src = src, existing_uids = .pstore_all_uids(x))
    if (!inherits(y, "unionParquetStore")) return(rbind2(x, y, src = src, ...))
  
    # inherit tracked params
    if (!is.null(x@fields) || !is.null(y@fields)) {
        x@fields <- unique(c(x@fields, y@fields))
    }
  
    x@stores <- c(x@stores, y@stores)
  
    initialize(x)
})

#' @rdname rbind
#' @export
setMethod("rbind2", signature("parquetGeomStore", "parquetGeomStore"), function(x, y, ...) {
    .guard_rbind_geom_type(x, y)
    .guard_pstore_lazy_ops(x, y)
    z <- callNextMethod(x, y, ...)
    .pstore_disk_extent(z) <- .rbind_union_disk_ext(x, y)
    z@geomtype <- x@geomtype
    z
})

#' @rdname rbind
#' @export
setMethod("rbind2", signature("unionParquetGeomStore", "parquetGeomStore"), function(x, y, ...) {
    .guard_rbind_geom_type(x, y)
    .guard_pstore_lazy_ops(x, y)
    z <- callNextMethod(x, y, ...)
    .pstore_disk_extent(z) <- .rbind_union_disk_ext(x, y)
    z
})

#' @rdname rbind
#' @export
setMethod("rbind2", signature("parquetGeomStore", "unionParquetGeomStore"), function(x, y, ...) {
    .guard_rbind_geom_type(x, y)
    .guard_pstore_lazy_ops(x, y)
    z <- callNextMethod(x, y, ...)
    .pstore_disk_extent(z) <- .rbind_union_disk_ext(x, y)
    z
})

#' @rdname rbind
#' @export
setMethod("rbind2", signature("unionParquetGeomStore", "unionParquetGeomStore"), function(x, y, ...) {
    .guard_rbind_geom_type(x, y)
    .guard_pstore_lazy_ops(x, y)
    z <- callNextMethod(x, y, ...)
    .pstore_disk_extent(z) <- .rbind_union_disk_ext(x, y)
    z
})

# internals ####

.guard_rbind_geom_classes <- function(x, y) {
    if (xor(inherits(x, "parquetGeomBase"),
            inherits(y, "parquetGeomBase"))) {
        stop("[rbind] parquetGeomStore can only be rbinded to other parquetGeomStore-inheriting",
            call. = FALSE)
    }
}

.guard_rbind_geom_type <- function(x, y) {
    gx <- x@geomtype                                                
    gy <- y@geomtype
    if (!length(gx) || !length(gy)) {                               
        stop("[rbind] geomtype not set on one or both stores", call. = FALSE)                                                           
    }                                                               
    if (gx != gy) {
        stop("[rbind] geomtype of items must be the same.", call. = FALSE)                                                              
    }
}

.rbind_dedup <- function(store, src, existing_uids) {
    test_ids <- if (inherits(store, "unionParquetStore")) {
        vapply(store@stores, function(s) s@uid, FUN.VALUE = character(1L))
    } else {
        store@uid
    }

    if (!any(test_ids %in% existing_uids)) return(store)
  
    create_class <- if (inherits(store, "parquetGeomBase")) {
        "parquetGeomStore"
    } else {
        "parquetStore"
    }
  
    if (!is.null(src)) return(sourceWrite(src, store, type = create_class))
    storeWrite(storeCreate(type = create_class), data = store)
}

.rbind_union_disk_ext <- function(x, y) {
    .ext_to_num_vec(
        terra::union(
            ext(.pstore_disk_extent(x)), 
            ext(.pstore_disk_extent(y))
        )
    )
}

.pstore_all_uids <- function(x) {                                   
    if (inherits(x, "unionParquetStore")) {                         
        unlist(lapply(x@stores, function(s) s@uid))                 
    } else {                                                        
        x@uid                               
    }                                                               
}      

.pstore_get_union_class <- function(x) {
    if (inherits(x, "parquetGeomBase")) {
        return("unionParquetGeomStore")
    } else {
        return("unionParquetStore")
    }
}

.test_pstore_lazy_ops <- function(x) {
    if (inherits(x, "parquetGeomBase")) {
        active_ext <- .pstore_active_extent(x)
    } else {
        active_ext <- NULL
    }
    !is.null(active_ext) || length(x@ops) > 0L
}

.guard_pstore_lazy_ops <- function(x, y) {
    has_ops <- vapply(list(x, y), 
        .test_pstore_lazy_ops, 
        FUN.VALUE = logical(1L)
    )
    if (!any(has_ops)) return(invisible())
  
    stop("[rbind] stores with pending ops cannot be combined.\n",
        "  Call `storeWrite(storeRead(x), path)` first to materialize.",
        call. = FALSE)
}
