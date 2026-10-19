#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-deconv ####
#
# analyzeData(parquetExprBase, dwlsDeconvParam) -- streaming spatialDWLS.
#
# The algorithm is the spatialDWLS package's, not a re-implementation:
# `spatialDWLS::dwls_stream()` drives both rounds and asks for one window of
# spots at a time through a callback. This file only supplies that callback
# over the store's cell windows. spatialDWLS is suggested, as it is by Giotto.
#
# Why it streams. Everything spatialDWLS reads is either a per-cluster or
# per-feature reduction (cluster means, PAGE reference levels, the best PAGE z
# per cluster, which types survive round 1) or a per-spot solve that reads its
# own column. Reductions fold across windows; solves never look outside their
# window.
#
# What it reads. Only the signature features, typically a few hundred out of
# the whole panel: the store is narrowed with `[` once, so every window is a
# features-subset read and the gene predicate prunes at the parquet level.
#
# Cost. Four passes over the narrowed store (means, z-scores, round 1,
# round 2). Each window is densified only per cluster inside the solver, at
# signature-feature width.
#
# Exactness. The per-spot solves are identical to the in-memory method. The
# reductions sum in window order rather than in one go, so cluster means agree
# to floating point; the dampening constant they feed is an argmin over a
# handful of candidates and is not moved by that.


#' @name analyzeData-dwlsDeconvParam
#' @title Streaming spatialDWLS deconvolution
NULL

#' @rdname analyzeData-dwlsDeconvParam
#' @description
#' [GiottoClass::analyzeData()] method for a [Giotto::dwlsDeconvParam-class]
#' on a disk-backed expression store. The store is narrowed to the signature
#' features and read one window of cells at a time, so memory scales with the
#' window, not with the number of spots. Needs the \pkg{spatialDWLS} package
#' (`remotes::install_github("drieslab/spatialDWLS")`).
#' @param x a `parquetExprBase` store of log-scaled expression.
#' @param param a [Giotto::dwlsDeconvParam-class].
#' @param sign_matrix mean expression per cell type, features x cell types,
#'   as built by [Giotto::makeSignMatrixDWLS()].
#' @param cell_clusters cluster of every cell, named by cell ID, or unnamed in
#'   the store's cell order.
#' @param ... additional arguments (none used).
#' @returns a `data.table` of `cell_ID` and one proportion column per cell
#'   type, the same contract the in-memory method returns.
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "dwlsDeconvParam"),
    function(x, param, ..., sign_matrix, cell_clusters) {
        .stream_dwls(pe = x, sign_matrix = sign_matrix,
            cell_clusters = cell_clusters, param = param)
    }
)


.stream_dwls <- function(pe, sign_matrix, cell_clusters, param) {
    if (!requireNamespace("spatialDWLS", quietly = TRUE)) {
        stop("[analyzeData(parquetExprBase, dwlsDeconvParam)] needs the ",
            "spatialDWLS package:\n",
            '  remotes::install_github("drieslab/spatialDWLS")',
            call. = FALSE)
    }
    cell_ids <- pe@cell_ids
    if (is.null(names(cell_clusters))) {
        if (length(cell_clusters) != length(cell_ids)) {
            stop("[analyzeData(parquetExprBase, dwlsDeconvParam)] unnamed ",
                "`cell_clusters` must have one value per cell of the store.",
                call. = FALSE)
        }
        names(cell_clusters) <- cell_ids
    }

    feat_ids <- pe@feat_ids
    genes <- intersect(rownames(sign_matrix), feat_ids)
    if (length(genes) == 0L) {
        stop("[analyzeData(parquetExprBase, dwlsDeconvParam)] no features ",
            "shared between `sign_matrix` and the store.", call. = FALSE)
    }

    # Narrow once; every window below is then a signature-features read.
    sub <- pe[genes, ]
    windows <- .pe_windows(sub,
        .pe_window_cells(sub, bytes_per_nz = 48, k = 0L))

    read_window <- function(i, features) {
        M <- storeRead(.pe_window_store(windows[[i]]), output = "dgcmatrix",
            max_rows = Inf, max_cols = Inf)
        M[features, , drop = FALSE]
    }

    res <- spatialDWLS::dwls_stream(
        read_window = read_window,
        n_windows = length(windows),
        feat_ids = feat_ids,
        clusters = cell_clusters,
        sign_matrix = sign_matrix,
        n_cell = param$n_cell,
        cutoff = param$cutoff,
        logbase = param$logbase,
        na_as_zero = param$na_as_zero,
        # the param's own, else the session's, as in Giotto's in-memory method
        BPPARAM = param$BPPARAM %||% getOption("giotto.bpparam"),
        chunk_size = param$chunk_size
    )

    out <- data.table::data.table(cell_ID = colnames(res))
    cbind(out, data.table::as.data.table(t(res)))[]
}
