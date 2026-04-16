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
#'  inlined automatically — no `!!` injection needed.
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
#' # local variable — inlined automatically
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
            expr <- .inline_local_vars(expr, c(colnames(x), specialCols(x)), env)
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

# head ####

setMethod("head", signature("parquetBase"), function(x, n = 6, ...) {
    x@ops <- c(x@ops, list(list(type = "head", n = n)))
    x
})

# tail ####

setMethod("tail", signature("parquetBase"), function(x, n = 6, ...) {
    x@ops <- c(x@ops, list(list(type = "tail", n = n)))
    x
})

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
    e <- ext(x) # data extent
    if (length(x@crop) > 0L) { # crop extent
        e <- terra::intersect(e, ext(x@crop))
    }
    if (length(x@window) > 0L) { # window extent
        e <- terra::intersect(e, ext(x@window))
    }
    e <- terra::intersect(e, ext(y)) # new crop extent
    if (is.null(e)) {
        stop("[crop] no geometries within requested extent\n", call. = FALSE)
    }
    x@crop <- .ext_to_num_vec(e)
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

    e <- ext(x) # data extent
    if (length(x@crop) > 0L) { # crop extent
        e <- terra::intersect(e, ext(x@crop))
    }
    e <- terra::intersect(e, ext(value)) # new window extent
    if (is.null(e)) {
        stop("[window] no geometries within requested extent\n", call. = FALSE)
    }
    x@window <- .ext_to_num_vec(e)
    x
})

.do_op <- function(atab, op) {
    type <- op$type
    switch(type,
        "filter" = dplyr::filter(atab, !!op$expr),
        "head"   = head(atab, op$n),
        "tail"   = tail(atab, op$n),
        "sample" = .arrow_sample_max_rows(atab, op$size),
        "distinct" = dplyr::distinct(atab, dplyr::across(dplyr::all_of(op$cols))),
        "join"   = {
            y_q <- storeRead(op$y, output = "query")
            # drop all of y's special cols except those used as join keys —
            # they belong to y's store context and would conflict or be
            # meaningless in the joined result
            y_drop <- setdiff(specialCols(op$y), unname(op$by))
            y_drop <- intersect(y_drop, names(y_q))
            if (length(y_drop) > 0L) {
                y_q <- dplyr::select(y_q, -dplyr::all_of(y_drop))
            }
            join_fn <- if (op$nomatch == "inner") dplyr::inner_join else dplyr::left_join
            join_fn(atab, y_q, by = op$by)
        },
        stop(sprintf("[.do_op] unknown op type: '%s'", type), call. = FALSE)
    )
}
