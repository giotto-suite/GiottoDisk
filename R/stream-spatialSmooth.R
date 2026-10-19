#' @include pkg_imports.R class-parquetExprStore.R class-parquetEdgeStore.R
NULL

# Spatial smoothing of an expression store, S = X M (GiottoClass defines
# spatialSmoothParam, the in-memory method and the operator M):
#
#   processData(parquetExprBase, spatialSmoothParam)  S written to a new store
#   analyzeData(parquetExprBase, featCorParam(smooth =))
#                                       correlation of S, never stored
#
# Both walk the same windows of smoothed values (.pe_smooth_prepare +
# .pe_smooth_block):
#
#   pairs in view       .edge_to_cell_map + .edge_view_pairs
#   local cell order    .pairs_local, else .pairs_bfs_order +
#                       .pe_rewrite_cell_order (which also bakes the chain)
#   operator            GiottoClass::spatial_smooth_operator, built directly in
#                       that order and from the pairs as they are (one per
#                       pair of cells), so no reordered or deduplicated copy
#   window              per range of target cells, the row range their
#                       neighbourhoods span, read with .pe_read_chunk_sub
#
# A smoothed store keeps the input's cell and feature order. It is denser
# than the input -- a cell gains every feature its neighbours have -- so
# smooth a feature subset (`x[feats, ]`) rather than a whole panel.

#' @rdname processData
#' @param network for a `spatialSmoothParam`: the spatial network as a
#'   `parquetEdgeStore`.
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "spatialSmoothParam"),
    function(x, param, ..., network) {
        prep <- .pe_smooth_prepare(x, param, network)
        on.exit(unlink(prep$tmp, recursive = TRUE), add = TRUE)
        .pe_smooth_write(prep, x@cell_ids, x@feat_ids)
    }
)

# The materialized source in a local cell order, the operator in that order
# and the windows of target cells.
.pe_smooth_prepare <- function(x, param, network) {
    if (!inherits(network, "parquetEdgeStore")) {
        stop("[spatialSmoothParam] `network` must be a parquetEdgeStore.",
            call. = FALSE)
    }
    n_cells <- as.integer(x@n_cells)
    pairs <- .edge_view_pairs(network, .edge_to_cell_map(network, x),
        isTRUE(param$weighted))
    budget <- .pe_smooth_window_cells(x, nrow(pairs))
    local <- budget >= n_cells || .pairs_local(pairs, budget)
    ord <- if (local) seq_len(n_cells) else .pairs_bfs_order(pairs, n_cells)
    inv <- integer(n_cells)
    inv[ord] <- seq_len(n_cells)
    Mo <- spatial_smooth_operator(inv[pairs$a], inv[pairs$b], pairs$w,
        n = n_cells, smoothing = param$smoothing, dedupe = FALSE)
    rm(pairs, inv)
    M1 <- Mo
    for (h in seq_len(param$hops - 1L)) Mo <- Mo %*% M1
    rm(M1)
    invisible(gc(verbose = FALSE))

    tmp <- .dump_tempfile()
    src <- if (local) {
        storeWrite(parquetExprStore(path = tmp), x)
    } else {
        .pe_rewrite_cell_order(x, ord, tmp)
    }
    list(src = src, info = list(sub = .exprbase_substores(src)[[1L]]$store),
        Mo = Mo, ord = ord, tmp = tmp,
        windows = .pe_chunk_ranges(1L, n_cells, budget))
}

# Target cells per window. A window holds its source values (read as a
# sparse matrix, ~150 bytes a value with the reader's buffers), the smoothed
# block -- about (neighbours + 1) times denser,
# up to every feature -- in sparse, triplet and product-temporary form, so
# the window is sized in bytes from those densities rather than from the
# sparse-read window. A pinned `giottodisk.chunk_size` keeps its meaning of
# a quarter-window, as for the other passes.
.pe_smooth_window_cells <- function(x, n_pairs) {
    n <- as.numeric(x@n_cells)
    pinned <- getOption("giottodisk.chunk_size")
    if (!is.null(pinned)) return(max(1L, as.integer(pinned) %/% 4L))
    free <- tryCatch(.sc_free_ram(), error = function(e) NA_real_)
    bytes <- if (isTRUE(is.finite(free)) && free > 0) {
        .chunk_budget(free, .ram_frac(), n_cells = 0, k = 0, p = 0)
    } else 2e9
    P <- length(x@feat_ids)
    dens <- tryCatch(.pe_view_density(x), error = function(e) 1)
    nz <- P * min(1, dens)
    s_nz <- min(P, nz * (2 * n_pairs / max(1, n) + 1))
    # the read as for crossprod(); the block in sparse, triplet and product
    # form
    per_cell <- 150 * nz + 12 * s_nz * 3
    max(1L, as.integer(min(n, bytes / per_cell)))
}

# Window `i` of smoothed values: S (features present x the window's target
# cells, in source order) and `vp`, those targets' positions in the view.
# NULL when the window's neighbourhoods hold no values.
.pe_smooth_block <- function(prep, i) {
    w <- prep$windows[[i]]
    tgt <- w[[1L]]:w[[2L]]
    blk <- prep$Mo[, tgt, drop = FALSE]
    rng <- range(blk@i + 1L)
    X <- .pe_read_chunk_sub(prep$info, rng[1L], rng[2L], list(),
        length(prep$src@feat_ids))
    if (is.null(X)) return(NULL)
    # columns by cell ID onto the full range: a cell without values in the
    # range has no column in the read
    Xt <- methods::as(X, "TsparseMatrix")
    cj <- match(colnames(X), prep$src@cell_ids[rng[1L]:rng[2L]])
    X <- Matrix::sparseMatrix(i = Xt@i + 1L, j = cj[Xt@j + 1L], x = Xt@x,
        dims = c(nrow(X), rng[2L] - rng[1L] + 1L),
        dimnames = list(rownames(X), NULL))
    list(S = X %*% blk[rng[1L]:rng[2L], , drop = FALSE], vp = prep$ord[tgt])
}

# S written under the cells' view positions, so the new store has the view's
# order.
.pe_smooth_write <- function(prep, cell_ids, feat_ids) {
    row_id <- col_id <- value <- NULL  # NSE
    n_cells <- length(cell_ids)
    fpos <- stats::setNames(seq_along(feat_ids), feat_ids)

    max_rows <- .ingest_batch_rows(96)
    dens <- tryCatch(.pe_view_density(prep$src), error = function(e) 1)
    n_rec <- as.numeric(n_cells) * length(feat_ids) * min(1, 4 * dens)
    n_buckets <- .ingest_n_buckets(n_rec, max_rows)
    per_bucket <- max(1L, as.integer(ceiling(n_cells / n_buckets)))

    produce <- function(i) {
        b <- .pe_smooth_block(prep, i)
        if (is.null(b)) return(NULL)
        S <- methods::as(b$S, "TsparseMatrix")
        if (!length(S@x)) return(NULL)
        vp <- b$vp[S@j + 1L]
        data.table::data.table(row_id = vp,
            col_id = unname(fpos[rownames(S)[S@i + 1L]]),
            value = S@x, b = (vp - 1L) %/% per_bucket + 1L)
    }
    finish <- function(dt) {
        dt <- dt[, list(row_id = as.integer(row_id), col_id = as.integer(col_id),
            value = as.numeric(value))]
        data.table::setkeyv(dt, c("row_id", "col_id"))
        dt
    }
    itr <- .ingest_spill_iterator(produce, length(prep$windows), n_buckets,
        finish, max_rows)
    itr$cell_ids <- function() cell_ids
    itr$feat_ids <- function() feat_ids
    itr$n_cells <- function() n_cells
    itr$n_genes <- function() length(feat_ids)

    path <- .dump_tempfile()
    out <- parquetExprStore(path = path)
    partition_dir <- .idpath(path, out@uid)
    dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)
    .pestore_write_itr(out, itr, partition_dir)
}

# x x^T and the row sums of the smoothed values, window by window: what
# featCorParam(smooth =) needs, without writing S anywhere.
.pe_smooth_moments <- function(prep, feat_ids) {
    P <- length(feat_ids)
    G <- matrix(0, P, P)
    s <- numeric(P)
    for (i in seq_along(prep$windows)) {
        b <- .pe_smooth_block(prep, i)
        if (is.null(b)) next
        k <- match(rownames(b$S), feat_ids)
        # dense sub-blocks: a dense product is far faster than a sparse one
        # on smoothed (dense) values, and the blocks bound its copy
        step <- max(1L, as.integer(2e7 %/% max(1L, length(k))))
        for (a in seq(1L, ncol(b$S), by = step)) {
            S <- as.matrix(b$S[, a:min(ncol(b$S), a + step - 1L), drop = FALSE])
            G[k, k] <- G[k, k] + tcrossprod(S)
            s[k] <- s[k] + rowSums(S)
        }
        rm(b, S)
    }
    list(G = G, s = s)
}
