#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-enrich ####
#
# analyzeData(parquetExprBase, pageEnrichParam | hyperEnrichParam) --
# streaming PAGE and hypergeometric enrichment.
#
# The statistics are Giotto's, not a re-implementation. Both score each cell
# from that cell's own values and one per-feature reference level shared by
# all cells, so:
#
#   pass 1  the reference level, on the feature axis;
#   pass 2  one window of cells at a time, read as a sparse matrix and scored
#           by Giotto's in-memory method against that reference.
#
# Windows hold every feature, because PAGE's per-cell sd and the
# hypergeometric per-cell cutoff are taken over all of them. Giotto densifies
# within a window in bounded blocks of cells, so memory follows the window.
#
# Scaling (`output_enrichment = "zscore"`) is across all cells, so windows
# return raw scores and the scaling is applied once at the end.
#
# Exactness: per-cell results equal the in-memory method's; the reference
# level is summed in window order rather than at once, which agrees to
# floating point.


#' @name analyzeData-pageEnrichParam
#' @title Streaming PAGE and hypergeometric enrichment
NULL

#' @rdname analyzeData-pageEnrichParam
#' @description
#' [GiottoClass::analyzeData()] methods for a [Giotto::pageEnrichParam-class]
#' or [Giotto::hyperEnrichParam-class] on a disk-backed expression store. The
#' store is read one window of cells at a time and each window is scored by
#' Giotto's in-memory method against the reference level of all cells, so
#' memory follows the window, not the number of cells. Rank enrichment is not
#' supported: it ranks each feature across all cells at once.
#' @param x a `parquetExprBase` store.
#' @param param a [Giotto::pageEnrichParam-class] or
#'   [Giotto::hyperEnrichParam-class].
#' @param sign_matrix sign matrix, features x cell types.
#' @param ... additional arguments (none used).
#' @returns a `data.table` of `cell_ID` and one column per cell type, the same
#'   contract the in-memory method returns.
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "pageEnrichParam"),
    function(x, param, ..., sign_matrix) {
        .stream_enrich(pe = x, sign_matrix = sign_matrix, param = param)
    }
)

#' @rdname analyzeData-pageEnrichParam
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "hyperEnrichParam"),
    function(x, param, ..., sign_matrix) {
        .stream_enrich(pe = x, sign_matrix = sign_matrix, param = param)
    }
)


.stream_enrich <- function(pe, sign_matrix, param) {
    ref <- .stream_enrich_reference(pe, param)

    # raw scores per window; the scaling across all cells comes last
    zscore <- identical(param$output_enrichment, "zscore") &&
        !isTRUE(param$p_value)
    inner <- param
    inner$output_enrichment <- "original"

    windows <- .pe_windows(pe, .pe_window_cells(pe, bytes_per_nz = 150, k = 0L))
    out <- data.table::rbindlist(lapply(seq_along(windows), function(i) {
        M <- storeRead(.pe_window_store(windows[[i]]), output = "dgcmatrix",
            max_rows = Inf, max_cols = Inf)
        # cell types are kept or dropped by feature counts alone, so report
        # that once, not per window
        if (i > 1L && !is.null(inner$verbose)) inner$verbose <- FALSE
        analyzeData(M, inner, sign_matrix = sign_matrix, reference = ref)
    }))

    if (zscore) {
        types <- setdiff(names(out), "cell_ID")
        out[, (types) := lapply(.SD, function(v) as.numeric(scale(v))),
            .SDcols = types]
    }
    out[]
}

# Each feature's reference level over all cells, as Giotto defines it: with
# reverse_log_scale the mean of the de-logged values logged back in the same
# base, else the plain mean.
.stream_enrich_reference <- function(pe, param) {
    n_cells <- length(pe@cell_ids)
    m <- if (isTRUE(param$reverse_log_scale)) {
        base <- param$logbase
        src <- .pe_push_op_last(pe, list(type = "expm1", base = base))
        log(.stream_expr_accum(src, axis = "feat", stats = "sum")$sum /
                n_cells + 1, base = base)
    } else {
        .stream_expr_accum(pe, axis = "feat", stats = "sum")$sum / n_cells
    }
    stats::setNames(m, pe@feat_ids)
}


# ---- rank: refused, with the reason -----------------------------------------

#' @rdname analyzeData-pageEnrichParam
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "rankEnrichParam"),
    function(x, param, ..., sign_matrix) {
        stop("[analyzeData(parquetExprBase, rankEnrichParam)] ",
            "rank enrichment is not supported on a disk-backed store. It ",
            "ranks each gene across every cell (sparseMatrixStats::rowRanks), ",
            "so no window of cells can be scored without the others -- ",
            "unlike PAGE and hypergeometric, which score each cell on its ",
            "own. Use one of those, or materialize the store first.",
            call. = FALSE)
    }
)
