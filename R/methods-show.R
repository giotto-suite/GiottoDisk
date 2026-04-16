# show ####
# Public show methods are thin entry points — all content building happens in
# .show_info(). With .print = FALSE, each method returns a named list of raw
# values; callNextMethod() is called with .print = FALSE so parent levels
# contribute without printing. The outermost (.print = TRUE) call collects the
# full list and hands it to .print_show(), which prints a <ClassName> header,
# routes the kv pairs through GiottoUtils::print_list() in canonical order,
# then appends the lazy-ops section at the bottom.

setMethod("show", signature("fileStore"), function(object) {
    .show_info(object, .print = TRUE)
    invisible(NULL)
})

setMethod("show", signature("unionParquetStore"), function(object) {
    .show_info(object, .print = TRUE)
    invisible(NULL)
})

setMethod("show", signature("h5ArrayStore"), function(object) {
    .show_info(object, .print = TRUE)
    invisible(NULL)
})

setMethod("show", signature("tileDBArrayStore"), function(object) {
    .show_info(object, .print = TRUE)
    invisible(NULL)
})


# .show_info ####

setGeneric(".show_info", function(object, .print = TRUE)
    standardGeneric(".show_info"))

setMethod(".show_info", signature("fileStore"), function(object, .print = TRUE) {
    info <- list()
    info[["class"]] <- class(object)
    info[["path"]]  <- paste0(str_abbreviate(object@path), collapse = "\n      ")
    if (!storeExists(object)) info[["status"]] <- "<empty>"
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("parquetStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    if (!storeExists(object)) {
        if (.print) return(.print_show(object, info))
        return(invisible(info))
    }
    fields <- colnames(object)
    atypes <- object@params$arrow_types
    if (!is.null(atypes)) {
        is_custom <- fields %in% names(object@datatype)
        is_int64  <- atypes[fields] == "int64"
        fields[is_int64 & !is_custom] <- color_red(fields[is_int64 & !is_custom])
    }
    internal_fields <- setdiff(specialCols(object), colnames(object))
    internal_info <- if (length(internal_fields) > 0L) {
        sprintf("[%s]", color_purple(paste(toString(internal_fields), "(internal)")))
    }
    info[["columns"]] <- paste(toString(fields), internal_info)
    info[["nrows"]]   <- if (length(object@ops) == 0L) {
        format(nrow(object), big.mark = ",", scientific = FALSE)
    } else "??"
    info[["ops"]] <- object@ops
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("unionParquetStore"), function(object, .print = TRUE) {
    info <- list()
    info[["class"]]     <- class(object)
    info[["substores"]] <- length(object@stores)
    fields <- colnames(object)
    internal_fields <- setdiff(specialCols(object), fields)
    internal_info <- if (length(internal_fields) > 0L) {
        sprintf("[%s]", color_purple(paste(toString(internal_fields), "(internal)")))
    }
    info[["columns"]] <- paste(toString(fields), internal_info)
    info[["nrows"]]   <- if (length(object@ops) == 0L) {
        format(nrow(object), big.mark = ",", scientific = FALSE)
    } else "??"
    info[["ops"]] <- object@ops
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("parquetGeomStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    if (!storeExists(object)) {
        if (.print) return(.print_show(object, info))
        return(invisible(info))
    }
    e <- if (!is.null(.pstore_active_extent(object))) {
        .ext_to_num_vec(.pstore_active_extent(object))
    } else if (!is.null(object@params$disk_extent) && length(object@ops) == 0L) {
        object@params$disk_extent
    } else {
        ext(object)
    }
    if (nzchar(object@geomtype)) info[["geomtype"]] <- object@geomtype
    info[["extent"]] <- .format_extent(e)
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("unionParquetGeomStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    if (length(object@stores) == 0L) {
        if (.print) return(.print_show(object, info))
        return(invisible(info))
    }
    e <- if (!is.null(.pstore_active_extent(object))) {
        .ext_to_num_vec(.pstore_active_extent(object))
    } else {
        ext(object)
    }
    if (nzchar(object@geomtype)) info[["geomtype"]] <- object@geomtype
    info[["extent"]] <- .format_extent(e)
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("parquetGeomTileStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    if (!storeExists(object)) {
        if (.print) return(.print_show(object, info))
        return(invisible(info))
    }
    info[["tiles"]] <- length(object@tiles)
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("h5ArrayStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    if (!storeExists(object)) {
        if (.print) return(.print_show(object, info))
        return(invisible(info))
    }
    info[["name"]] <- object@params$name
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("tileDBArrayStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    if (!storeExists(object)) {
        if (.print) return(.print_show(object, info))
        return(invisible(info))
    }
    info[["name"]] <- object@params$name
    if (.print) return(.print_show(object, info))
    invisible(info)
})


# internals ####

.show_key_order <- c(
    "path", "substores", "geomtype", "extent",
    "columns", "nrows", "tiles", "name", "status"
)

.print_show <- function(object, info) {
    cat(sprintf("<%s>\n", info[["class"]]))
    kv <- info[intersect(.show_key_order, names(info))]
    if (length(kv) > 0L) GiottoUtils::print_list(kv)
    if ("ops" %in% names(info)) {
        ops <- info[["ops"]]
        cat("\nlazy ops:\n")
        if (length(ops) == 0L) cat("<none>\n")
        else for (op in ops) cat(.format_op(op))
    }
    invisible(NULL)
}

.format_extent <- function(x) {
    if (inherits(x, "SpatExtent")) x <- .ext_to_num_vec(x)
    sprintf("%s (xmin, xmax, ymin, ymax)", toString(round(x, digits = 3)))
}

.format_op <- function(step) {
    args <- switch(step$type,
        "filter" = {
            expr_str <- deparse(step$expr,
                width.cutoff = min(getOption("width", 50L), 50) - 12L)
            if (length(expr_str) > 1L) expr_str <- paste0(expr_str[[1L]], " ...")
            expr_str
        },
        "head"   = format(step$n),
        "tail"   = format(step$n),
        "sample" = sprintf("size = %g", step$size),
        "select"   = toString(step$cols),
        "distinct" = toString(step$cols),
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
