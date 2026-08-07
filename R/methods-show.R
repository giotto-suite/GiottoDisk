# show ####
# Public show methods are thin entry points -- all content building happens in
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

setMethod("show", signature("parquetExprStore"), function(object) {
    .show_info(object, .print = TRUE)
    invisible(NULL)
})

setMethod("show", signature("unionParquetExprStore"), function(object) {
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
    if (length(object@uid) > 0L && nzchar(object@uid)) info[["uid"]] <- object@uid
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
    # exact = FALSE: estimated extent -- no scan, projects through pending
    # transform if any. Appropriate for display; row-filter effects not reflected.
    e <- .ext_to_num_vec(ext(object, exact = FALSE))
    if (length(object@geomtype) > 0L && nzchar(object@geomtype)) info[["geomtype"]] <- object@geomtype
    info[["extent"]] <- .format_extent(e)
    info[["ops"]] <- c(info[["ops"]], object@post_ops)
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("unionParquetGeomStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    if (length(object@stores) == 0L) {
        if (.print) return(.print_show(object, info))
        return(invisible(info))
    }
    e <- .ext_to_num_vec(ext(object, exact = FALSE))
    if (length(object@geomtype) > 0L && nzchar(object@geomtype)) info[["geomtype"]] <- object@geomtype
    info[["extent"]] <- .format_extent(e)
    info[["ops"]] <- c(info[["ops"]], object@post_ops)
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

.pe_id_preview <- function(ids, n = 3L) {
    n_show <- min(n, length(ids))
    paste0(
        paste(ids[seq_len(n_show)], collapse = ", "),
        if (length(ids) > n) ", ..." else ""
    )
}

setMethod(".show_info", signature("parquetExprStore"), function(object, .print = TRUE) {
    info <- callNextMethod(object, .print = FALSE)
    info[["dim"]] <- sprintf("%s genes x %s cells",
        format(object@n_genes, big.mark = ",", scientific = FALSE),
        format(object@n_cells, big.mark = ",", scientific = FALSE))
    if (length(object@feat_ids) > 0L) {
        info[["feat_ids"]] <- .pe_id_preview(object@feat_ids)
    }
    if (length(object@cell_ids) > 0L) {
        info[["cell_ids"]] <- .pe_id_preview(object@cell_ids)
    }
    if (length(object@cell_idx) > 0L || length(object@gene_idx) > 0L) {
        info[["subset"]] <- sprintf("cell_idx[%d] gene_idx[%d]",
            length(object@cell_idx), length(object@gene_idx))
    }
    if (length(object@ops) > 0L) {
        info[["jit_ops"]] <- paste(
            vapply(object@ops, function(op) op$type, character(1L)),
            collapse = " -> "
        )
    }
    # Marginals are a cached property of the file; the streaming window is
    # derived from them per read, so there is no stored chunk size to print.
    nnz <- .pestore_view_nnz(object)
    if (isTRUE(is.finite(nnz))) {
        info[["nonzeros"]] <- sprintf("%s (%.1f%% dense)",
            format(round(nnz), big.mark = ",", scientific = FALSE),
            100 * nnz / max(as.numeric(object@n_cells) *
                            as.numeric(object@n_genes), 1))
    }
    if (.print) return(.print_show(object, info))
    invisible(info)
})

setMethod(".show_info", signature("unionParquetExprStore"), function(object, .print = TRUE) {
    info <- list()
    info[["class"]] <- "unionParquetExprStore"
    info[["substores"]] <- length(object@stores)
    info[["dim"]] <- sprintf("%s genes x %s cells",
        format(object@n_genes, big.mark = ",", scientific = FALSE),
        format(object@n_cells, big.mark = ",", scientific = FALSE))
    if (length(object@feat_ids) > 0L) {
        info[["feat_ids"]] <- .pe_id_preview(object@feat_ids)
    }
    if (length(object@cell_ids) > 0L) {
        info[["cell_ids"]] <- .pe_id_preview(object@cell_ids)
    }
    if (.print) return(.print_show(object, info))
    invisible(info)
})


# internals ####

.show_key_order <- c(
    "path", "uid", "substores", "geomtype", "extent",
    "columns", "nrows", "tiles", "name",
    "dim", "feat_ids", "cell_ids", "subset", "jit_ops", "chunk",
    "status"
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
        "spat_relate" = {
            y_str <- if (!is.null(step$y_wkt)) {
                wkt <- step$y_wkt
                if (nchar(wkt) > 40L) paste0(substr(wkt, 1L, 37L), "...") else wkt
            } else {
                "<store>"
            }
            sprintf("%s [%s] %s", step$relation, step$form, y_str)
        },
        "multiply" = sprintf("x %s", .format_axis_payload(step$factors, step$axis)),
        "add"      = sprintf("+ %s", .format_axis_payload(step$terms, step$axis)),
        "log"           = sprintf("base = %g", step$base %null% 2),
        "transform" = {
            aff <- step$affine2d
            parts <- character(0L)
            if (!isTRUE(all.equal(unname(aff@rotate), 0)))
                parts <- c(parts, sprintf("rotate=%.4g", aff@rotate))
            if (!isTRUE(all.equal(unname(aff@shear), c(0, 0))))
                parts <- c(parts, sprintf("shear=(%s)", paste(round(aff@shear, 4), collapse = ",")))
            if (!isTRUE(all.equal(unname(aff@scale), c(1, 1))))
                parts <- c(parts, sprintf("scale=(%s)", paste(round(aff@scale, 4), collapse = ",")))
            if (!isTRUE(all.equal(unname(aff@translate), c(0, 0))))
                parts <- c(parts, sprintf("translate=(%s)", paste(round(aff@translate, 4), collapse = ",")))
            if (length(parts) == 0L) "identity" else paste(parts, collapse = " ")
        },
        "..."
    )
    sprintf("  %-10s: %s\n", step$type, args)
}


# Compact description of a multiply / add payload for show().
.format_axis_payload <- function(payload, axis) {
    if (is.null(payload)) return("?")
    if (!is.list(payload)) return(format(as.numeric(payload), digits = 4))
    n <- sum(vapply(payload, function(v) sum(!is.na(v)), numeric(1L)))
    sprintf("%s %s", format(n, big.mark = ","),
            if (identical(axis, "feat")) "feats" else "cells")
}
