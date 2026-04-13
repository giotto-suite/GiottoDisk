setMethod("show", signature("fileStore"), function(object) {
    cat(sprintf("<%s>\n", class(object)))
    cat(sprintf("path: %s\n", 
        paste0(str_abbreviate(object@path), collapse = "\n      ")
    ))
    if (!storeExists(object)) {
        cat("<empty>\n")
        return(invisible())
    }
})

setMethod("show", signature("parquetStore"), function(object) {
    callNextMethod(object)
    if (!storeExists(object)) return(invisible())
    .pbase_show(object)
})

setMethod("show", signature("unionParquetStore"), function(object) {
    cat(sprintf("<%s>\n", class(object)))
    cat(sprintf("substores: %s\n", length(object@stores)))
    .pbase_show(object)
})

setMethod("show", signature("unionParquetGeomStore"), function(object) {
    callNextMethod(object) # unionParquetStore show
    if (length(object@stores) == 0L) return(invisible())
    if (!is.null(.pstore_active_extent(object))) {
        e <- .ext_to_num_vec(.pstore_active_extent(object))
    } else {
        e <- ext(object)
    }
    cat(.format_extent(e))
})

.pbase_show <- function(object) {
    fields <- colnames(object)
  
    # check for type compat
    atypes <- object@params$arrow_types
    if (!is.null(atypes)) {
        is_custom <- fields %in% names(object@datatype)
        is_int64 <- atypes[fields] == "int64"
        bool <- is_int64 & !is_custom
        fields[bool] <- color_red(fields[bool])
    }
  
    internal_fields <- setdiff(specialCols(object), fields)
    fields_info <- toString(fields)
    internal_fields_info <- if (length(internal_fields) > 0L) {
        paste(toString(internal_fields), "(internal)") |>
            color_purple() |>
            sprintf(fmt = "[%s]")
    } else {
        NULL
    }
    fields_info <- paste(fields_info, internal_fields_info)

    cat((sprintf("columns: %s\n", fields_info)))
    if (length(object@ops) == 0L) {
        nr <- format(nrow(object), big.mark = ",", scientific = FALSE)
    } else {
        nr <- "??"
    }
    cat(sprintf("nrows: %s\n", nr))
  
    cat("\nlazy ops:\n")
    if (length(object@ops) == 0L) {
        cat("<none>\n")
    } else {
        for (step in object@ops) {
            cat(.format_op(step))
        }
    }
}

setMethod("show", signature("parquetGeomStore"), function(object) {
    callNextMethod(object)
    if (!storeExists(object)) return(invisible())
    if (!is.null(.pstore_active_extent(object))) {
        e <- .ext_to_num_vec(.pstore_active_extent(object))
    } else if (!is.null(object@params$disk_extent) && length(object@ops) == 0L) {
        e <- object@params$disk_extent
    } else {
        e <- ext(object)
    }
  
    cat(.format_extent(e))
})

.format_extent <- function(x) {
    if (inherits(x, "SpatExtent")) x <- .ext_to_num_vec(x)
    sprintf("extent: %s (xmin, xmax, ymin, ymax)\n",
        toString(round(x, digits = 3))
    )
}

setMethod("show", signature("parquetGeomTileStore"), function(object) {
    callNextMethod(object)
    if (!storeExists(object)) return(invisible())
    cat(sprintf("tiles: %d\n", length(object@tiles)))
})

setMethod("show", signature("h5ArrayStore"), function(object) {
    callNextMethod(object)
    if (!storeExists(object)) return(invisible())
    cat(sprintf("name: \"%s\"\n", object@params$name))
})

setMethod("show", signature("tileDBArrayStore"), function(object) {
    callNextMethod(object)
    if (!storeExists(object)) return(invisible())
    cat(sprintf("name: \"%s\"\n", object@params$name))
})

# internals ####

.format_op <- function(step) {
    args <- switch(step$type,
        "filter" = {
            expr_str <- deparse(step$expr,
                width.cutoff = min(getOption("width", 50L), 50) - 12L)
            if (length(expr_str) > 1L) {
                expr_str <- paste0(expr_str[[1L]], " ...")
            }
            expr_str
        },
        "head"   = format(step$n),
        "tail"   = format(step$n),
        "sample" = sprintf("size = %g", step$size),
        "select" = toString(step$cols),
        "join"   = {
            type_str <- step$nomatch
            keys_str <- paste(names(step$by), unname(step$by), sep = " = ",
                collapse = ", ")
            sprintf("%s on [%s]", type_str, keys_str)
        },
        "..."
    )
    sprintf("  %-8s: %s\n", step$type, args)
}