#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-enrich ####
#
# analyzeData(parquetExprBase, pageEnrichParam | hyperEnrichParam |
# rankEnrichParam) -- streaming sign-matrix enrichment.
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
#' @title Streaming sign-matrix enrichment
NULL

#' @rdname analyzeData-pageEnrichParam
#' @description
#' [GiottoClass::analyzeData()] methods for a [Giotto::pageEnrichParam-class]
#' or [Giotto::hyperEnrichParam-class] on a disk-backed expression store. The
#' store is read one window of cells at a time and each window is scored by
#' Giotto's in-memory method against the reference level of all cells, so
#' memory follows the window, not the number of cells.
#'
#' A [Giotto::rankEnrichParam-class] runs the same way, except that rank
#' enrichment first ranks each feature across all cells. Those ranks are
#' approximated from a per-feature histogram of the values, so scores agree
#' closely with the in-memory method rather than exactly. Rank `p_value =
#' TRUE` is not supported on a store.
#' @param x a `parquetExprBase` store.
#' @param param a [Giotto::pageEnrichParam-class],
#'   [Giotto::hyperEnrichParam-class] or [Giotto::rankEnrichParam-class].
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
    .stream_enrich_windows(pe, sign_matrix, param,
        window_args = function(M) list(reference = ref))
}

# Score one window of cells at a time with Giotto's in-memory method, handing
# it what each window needs from all cells (`window_args(M)`). Windows return
# raw scores (and with p_value their p-values, the scores as an attribute);
# scaling (`output_enrichment = "zscore"`) is across all cells, so it comes
# last. `dense_bytes` is the dense working memory per value of a
# window (features x cells), on top of reading it.
.stream_enrich_windows <- function(pe, sign_matrix, param, window_args,
                                   dense_bytes = 0) {
    zscore <- identical(param$output_enrichment, "zscore")
    inner <- param
    inner$output_enrichment <- "original"

    dens <- tryCatch(.pe_view_density(pe), error = function(e) 1)
    per_nz <- 150 + dense_bytes / max(min(1, dens), 1e-6)
    windows <- .pe_windows(pe, .pe_window_cells(pe, bytes_per_nz = per_nz,
                                                k = 0L))
    res <- lapply(seq_along(windows), function(i) {
        M <- storeRead(.pe_window_store(windows[[i]]), output = "dgcmatrix",
            max_rows = Inf, max_cols = Inf)
        # cell types are kept or dropped by feature counts alone, so report
        # that once, not per window
        if (i > 1L && !is.null(inner$verbose)) inner$verbose <- FALSE
        do.call(analyzeData,
            c(list(M, inner, sign_matrix = sign_matrix), window_args(M)))
    })
    out <- data.table::rbindlist(res)
    scores <- if (isTRUE(param$p_value)) {
        data.table::rbindlist(lapply(res, attr, "scores"))
    } else {
        out
    }

    if (zscore) {
        types <- setdiff(names(scores), "cell_ID")
        scores[, (types) := lapply(.SD, function(v) as.numeric(scale(v))),
            .SDcols = types]
    }
    if (isTRUE(param$p_value)) data.table::setattr(out, "scores", scores[])
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


# ---- rank: feature ranks from a histogram -----------------------------------
#
# Rank enrichment ranks each feature across all cells, then ranks those ranks
# within each cell. The first step needs every cell at once, so it is
# approximated: one pass counts each feature's values into equal bins over
# (0, max] (`.pe_value_hist()`), and a value's rank is the count below its bin
# plus its share of the bin. The score uses these ranks only through their
# order within a cell, so it stays close to the in-memory result; ties and
# values sharing a bin are where it differs. Windows of cells are then scored
# by Giotto's method from these ranks.

#' @rdname analyzeData-pageEnrichParam
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "rankEnrichParam"),
    function(x, param, ..., sign_matrix) {
        .stream_rank(pe = x, sign_matrix = sign_matrix, param = param)
    }
)

.stream_rank <- function(pe, sign_matrix, param, n_bins = 1024L) {
    if (isTRUE(param$p_value)) {
        stop("[analyzeData(parquetExprBase, rankEnrichParam)] ",
            "`p_value = TRUE` is not supported on a disk-backed store: its ",
            "null is fitted to the scores of every cell for n_times random ",
            "signatures. Use the scores, or an in-memory object.",
            call. = FALSE)
    }
    h <- .rank_hist(pe, n_bins)
    ties <- if (identical(param$ties_method, "max")) "min" else "average"
    # each window holds its dense rank matrix, about two copies while built
    .stream_enrich_windows(pe, sign_matrix, param, dense_bytes = 16,
        window_args = function(M) list(ranks = .rank_from_hist(
            M, h, match(rownames(M), pe@feat_ids), ties)))
}

# The histogram as lookup tables: per feature its max, the cells at or below
# 0, and per bin the count and the count of all lower bins.
.rank_hist <- function(pe, n_bins) {
    h <- .pe_value_hist(pe, n_bins = n_bins)
    G <- length(h$max)
    n <- matrix(0, G, n_bins)
    n[cbind(h$counts$pos, h$counts$bin + 1L)] <- h$counts$n
    below <- cbind(0, t(apply(n, 1L, cumsum)))[, seq_len(n_bins), drop = FALSE]
    list(n = n, below = below, max = h$max, n_zero = h$n_zero,
         n_bins = n_bins)
}

# Each feature's approximate rank across all cells, for the cells of `M`
# (features x cells, dgCMatrix). `pos` maps the rows of `M` into the
# histogram. Ranks follow `ties`: "average", or "min" for the max-ties rule.
.rank_from_hist <- function(M, h, pos, ties) {
    if (length(M@x) && min(M@x) < 0) {
        stop("[rank enrichment] disk-backed rank enrichment needs ",
            "non-negative expression values.", call. = FALSE)
    }
    nz <- h$n_zero[pos]
    zero <- if (ties == "average") (nz + 1) / 2 else rep(1, length(pos))
    R <- matrix(zero, nrow(M), ncol(M), dimnames = dimnames(M))
    on <- M@x > 0   # stored zeros keep the zero rank
    if (!any(on)) return(R)
    row <- (M@i + 1L)[on]
    col <- rep(seq_len(ncol(M)), diff(M@p))[on]
    g <- pos[row]
    s <- M@x[on] * h$n_bins / h$max[g]
    b <- pmin(pmax(ceiling(s) - 1, 0), h$n_bins - 1)
    k <- cbind(g, b + 1L)
    R[cbind(row, col)] <- h$n_zero[g] +
        h$below[k] + (s - b) * h$n[k] + if (ties == "average") 0.5 else 1
    R
}
