# subset ####

#' @name subset
#' @title Subset a parquet store
#' @description                                            
#' Filter rows of a parquet store by a logical expression. The 
#' operation is recorded lazily in the store's `@ops` slot and 
#' applied at read time.                                      
#'
#' Expressions are evaluated using non-standard evaluation. Local 
#' variables referenced in the expression are automatically inlined
#' at capture time, making the recorded operation self-contained 
#' and safe to pass across sessions or to parallel workers.    
#' @param x `parquetBase`-inheriting store object
#' @param subset logical `expression` to filter rows. Column names
#'  in the store are referenced directly. Local variables are 
#'  inlined automatically -- no `!!` injection needed.
#' @param select `expression`, indicating columns to keep. `-` can
#'  be used to drop columns.
#' @param negate `logical`. If `TRUE`, the filter expression is
#'  negated, keeping rows/cols that do NOT satisfy the condition.
#' @param quote `logical`. If `TRUE` (default), `subset` is captured
#'  via NSE. If `FALSE`, `subset` is treated as a pre-built R `call` 
#'  object, allowing programmatic construction of filter expressions.
#' @param ... additional arguments (ignored)
#' @returns the store with the filter step appended to `@ops`
#' @examples
#' # standard NSE
#' subset(store, gene == "EPCAM")
#'
#' # local variable -- inlined automatically
#' my_genes <- c("EPCAM", "CDH1")
#' subset(store, gene %in% my_genes)
#'
#' # negation
#' subset(store, gene == "EPCAM", negate = TRUE)
#'
#' # pre-built expression
#' expr <- quote(gene == "EPCAM")
#' subset(store, expr, quote = FALSE)
NULL

#' @rdname subset
#' @export
setMethod("subset", signature("parquetBase"), function(x, subset, select,
    negate = FALSE, quote = TRUE, ...) {
      
    if (!missing(subset)) {
        if (quote) {
            q <- rlang::enquo(subset)
            expr <- rlang::quo_get_expr(q)
            env <- rlang::quo_get_env(q)
            # Use effective_schema (NOT colnames(x) or disk_fields) as the
            # col universe -- a prior `[, j]` may have narrowed @fields,
            # and queued join ops bring in y-side cols. Filter exprs need
            # to reference any of these without being treated as local
            # vars. The lazy_fields path re-widens at storeRead time so
            # the upstream projection includes referenced cols, and
            # .apply_op's join handler materializes y-side cols post-join.
            expr <- .inline_local_vars(expr, .pstore_effective_schema(x), env)
        } else {
            expr <- subset  # pre-quoted expression passed directly
        }
        if (negate) expr <- call("!", expr)
        x@ops <- c(x@ops, list(list(type = "filter", expr = expr)))
    }
  
    if (!missing(select)) {
        nl <- as.list(seq_along(colnames(x)))
        names(nl) <- colnames(x)
        q <- rlang::enquo(select)
        vars <- eval(rlang::quo_get_expr(q), nl, rlang::quo_get_env(q))
        select_cols <- colnames(x)[vars]
        if (negate) select_cols <- setdiff(colnames(x), select_cols)
        x <- x[, select_cols]
    }

    x
})

## internals ####

# Walk an R expression AST and inline any symbols that are not column names
# by evaluating them in the caller's environment. This ensures the recorded
# operation is self-contained with no external variable references.
.inline_local_vars <- function(expr, schema_cols, env) {           
    if (is.symbol(expr)) { # parse referenced vars
        nm <- as.character(expr)
        if (!nm %in% schema_cols)
            # convert referenced var to inline values if not coldata
            return(eval(expr, envir = env))
        # return call with no changes if coldata requested
        return(expr)
    }
    if (is.call(expr)) { # parse call
        expr[-1] <- lapply(expr[-1], .inline_local_vars,
            schema_cols = schema_cols, env = env)
        return(expr)
    }
    expr # literals pass through
}

# unique ####

#' @export
setMethod("unique", signature("parquetBase"), function(x, incomparables = FALSE, ...) {
    x@ops <- c(x@ops, list(list(type = "distinct", cols = colnames(x))))
    x
})

# head / tail ####
# S3 methods, not S4. setMethod("head", "parquetBase", ...) against utils's
# implicit S3 generic auto-promotes utils::head to S4 inside this namespace,
# which produces "generic for X loaded with a different signature" warnings
# whenever a downstream package that also touches utils::head loads. Bare
# `head(pgs)` from user scope also doesn't dispatch through the S4 method
# because utils::head is found first on the search path; it falls into
# head.default → head.array → `[`, which errors on the S4 object.
#
# Following GiottoClass's pattern (head.giottoBinPoints + @exportS3Method
# utils::head): S3 dispatch walks S4 class() inheritance, so a method on
# parquetBase still fires for concrete subclasses (parquetGeomStore,
# parquetGeomTileStore, parquetStore, parquetExprStore, ...).

#' @title Head and tail
#' @name parquet-headtail
#' @description Queue a `head` or `tail` op on a lazy parquet store. The
#'   op is applied at `storeRead()` time.
#' @param x parquetBase-inheriting store
#' @param n integer. Number of rows to keep.
#' @param ... additional arguments (ignored)
#' @returns the input store with a head / tail op queued on `@ops`
#' @exportS3Method utils::head
head.parquetBase <- function(x, n = 6L, ...) {
    x@ops <- c(x@ops, list(list(type = "head", n = n)))
    x
}

#' @rdname parquet-headtail
#' @exportS3Method utils::tail
tail.parquetBase <- function(x, n = 6L, ...) {
    x@ops <- c(x@ops, list(list(type = "tail", n = n)))
    x
}

# rowSample ####
#' @name rowSample
#' @title Subsample rows of a parquetStore
#' @description
#' Record a row subsampling operation on a store. The operation is lazy and
#' applied when reading via [storeRead()]. If the store already has fewer
#' rows than `size`, the store is returned unchanged.
#'
#' Sampling is performed by systematic (evenly-spaced) selection rather than
#' random sampling. Given a target of `size` rows from `n` total rows, every
#' `k`th row is selected where `k = ceiling(n / size)`. This is
#' deterministic and reproducible across sessions without a seed.
#' 
#' The actual number of rows returned is approximate since the rows to be
#' kept are determined based the value of the internal `row_index` col.
#' Depending on any preceding [subset()] operations, the number of rows
#' surviving the sample filter may be slightly above or below `size`.
#' @param x `parquetBase`-inheriting object
#' @param size `numeric(1)`. Maximum number of rows to return
#' @param ... additional params (none implemented)
#' @returns `x` with sampling op appended to `@ops`
#' @seealso [storeRead()], [subset()]
#' @examples
#' ps <- parquetStore()
#' ps <- storeWrite(ps, mtcars)
#' ps |> rowSample(10)         # lazy, no disk read
#' ps |> rowSample(10) |> storeRead(output = "tibble")
#' @export 
setMethod("rowSample", signature("parquetBase"), function(x, size, ...) {
    x@ops <- c(x@ops, list(list(type = "sample", size = size)))
    x
})

# crop ####

#' @name crop
#' @title Crop a `parquetGeomStore`
#' @description
#' Apply a spatial crop on geometry data. The crop is more accurately
#' an xy selection on centroids. Polygon geometries are not modified.
#' 
#' Crop operations are composable operations that limit where future
#' crops or [window()] can be placed.
#' @param x object to crop
#' @param y spatial extent to crop to. Accepts any object that works
#'   with `ext()`
#' @param ... additional params to pass (none implemented)
#' @returns `parquetGeomBase`-inheriting object
#' @export
#' @seealso [window()]
setMethod("crop", signature("parquetGeomBase", "ANY"), function(x, y, ...) {
    r <- .pgeom_resolve_extent(x, y)
    x <- r$x; y <- r$e; e <- r$base_e
    if (length(x@crop) > 0L)   e <- terra::intersect(e, ext(x@crop))
    if (length(x@window) > 0L) e <- terra::intersect(e, ext(x@window))
    e <- terra::intersect(e, ext(y))
    if (is.null(e)) stop("[crop] no geometries within requested extent\n", call. = FALSE)
    x@crop <- .ext_to_num_vec(e)
    x
})

# * crop (tile store) ####

setMethod("crop", signature("parquetGeomTileStore", "ANY"), function(x, y, ...) {
    x <- callNextMethod()  # parquetGeomBase: resolves transforms, composes @crop
    # When a transform is pending, back-project the original query polygon through
    # the inverse affine and use the exact parallelogram for tile intersection.
    # This avoids the AABB over-inclusiveness that would pull in corner tiles outside
    # the actual rotated/sheared crop region.
    # When no transform is pending, @crop is already exact — use it directly.
    aff <- .pgeom_pending_transform(x)
    query_region <- if (!is.null(aff)) {
        affine(terra::as.polygons(.clamp_ext_infinite(terra::ext(y))), aff, inv = TRUE)
    } else {
        terra::ext(x@crop)
    }
    tile_sel <- tilework::intersect(x@tiles, query_region)
    x@tile_filter <- as.integer(tile_sel$tile)
    x
})

# * window ####

#' @name window
#' @aliases `window<-`
#' @title Window a `parquetGeomStore`
#' @description
#' Similar to [crop()], but does not apply a permanent spatial subset
#' on the data. 
#' @param x object to window
#' @param ... additional params to pass (none implemented)
#' @param value extent to apply as a window. Can be any object that
#'   works with `ext()` or `NULL` to remove the window
NULL

#' @rdname window
#' @export
setMethod("window", signature("parquetGeomBase"), function(x, ...) {
    !length(x@window) == 0L
})

#' @rdname window
#' @export
setMethod("window<-", signature("parquetGeomBase"), function(x, ..., value) {
    if (is.null(value)) { # remove window
        x@window <- numeric(0L)
        return(x)
    }

    r <- .pgeom_resolve_extent(x, value)
    x <- r$x; value <- r$e; e <- r$base_e
    if (length(x@crop) > 0L) e <- terra::intersect(e, ext(x@crop))
    e <- terra::intersect(e, ext(value))
    if (is.null(e)) stop("[window] no geometries within requested extent\n", call. = FALSE)
    x@window <- .ext_to_num_vec(e)
    x
})

# Back-project e through any pending affine into intrinsic space, injecting a
# half-plane filter into x@ops when rotation/shear is pending. Also resolves
# the tightest available intrinsic baseline (disk_extent or live scan).
# Returns list(x, e, base_e): modified store, intrinsic AABB, intrinsic baseline.
.pgeom_resolve_extent <- function(x, e) {
    ops_before <- length(x@ops)
    aff <- .pgeom_pending_transform(x)
    if (!is.null(aff)) {
        # Clamp +/-Inf: terra::as.polygons() cannot handle infinite coordinates.
        inv_poly <- affine(terra::as.polygons(.clamp_ext_infinite(ext(e))), aff, inv = TRUE)
        corners  <- terra::crds(inv_poly)[seq(4L), ]
        if (.affine_has_rotation(aff)) {
            x@ops <- c(x@ops, list(list(type = "filter",
                expr = .affine_halfplane_expr(corners))))
        }
        e <- ext(.affine_aabb(corners))
    }
    base_e <- if (!is.null(.pstore_disk_extent(x)) && ops_before == 0L) {
        ext(.pstore_disk_extent(x))
    } else {
        .pgeom_ext_intrinsic(x)
    }
    list(x = x, e = ext(e), base_e = base_e)
}

# apply arrow-compatible operations
#
# op chains in stores store operations as an ordered list of
# param lists. Acero-compatible operations can be performed by passing
# a pristine arrow table through .apply_op with the ops in sequence to
# apply the lazy steps.
#
# `storeRead()` methods control and enforce the application of ops.
# This helper should only be used in low-level `storeRead` methods
# that define how to apply ops to the object being read from.
#
# atab : arrow table
# op   : op item (usually from store@ops)

.ptabular_apply_op <- function(atab, op) {
    type <- op$type
    switch(type,
        "filter" = dplyr::filter(atab, !!op$expr),
        "head"   = head(atab, op$n),
        "tail"   = tail(atab, op$n),
        "sample" = .arrow_sample_max_rows(atab, op$size),
        "distinct" = .op_distinct(atab, op),
        "join"   = .op_join(atab, op),
        # `spat_relate` is handled at the `.pbase_storeread_processing`
        # level (not via `.apply_op`) so it can route through sedonadb.
        # Reaching here means something bypassed the processing loop.
        "spat_relate" = stop(
            "[.apply_op] 'spat_relate' must be handled by ",
            ".pbase_storeread_processing; reached .apply_op unexpectedly",
            call. = FALSE),
        # Internal op type created ephemerally by the spat_relate
        # evaluation path to inject cached surviving ids without re-running
        # the spatial predicate. Carries an arrow Table of id cols and the
        # join keys. Not part of the public API.
        "id_filter" = dplyr::semi_join(atab, op$ids_tab, by = op$by),
        # fallback error
        stop(sprintf("[.apply_op] unknown op type: '%s'", type), call. = FALSE)
    )
}

.op_distinct <- function(atab, op) {
    dplyr::distinct(atab, dplyr::across(dplyr::all_of(op$cols)))
}

.op_join <- function(atab, op) {
    y_q <- storeRead(op$y, output = "query")
    # drop y's special cols except join keys
    y_drop <- setdiff(specialCols(op$y), unname(op$by))
    y_drop <- intersect(y_drop, names(y_q))
    if (length(y_drop) > 0L) {
        y_q <- dplyr::select(y_q, -dplyr::all_of(y_drop))
    }
    join_fn <- if (op$nomatch == "inner") {
        dplyr::inner_join 
    } else {
        dplyr::left_join
    }
    join_fn(atab, y_q, by = op$by)
}
