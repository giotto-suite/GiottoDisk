# Giotto's neighborhoodPCAParam (BANKSY embedding) on expression stores, in
# two passes over windows of target cells: the joint matrix (each feature, its
# neighbour mean and, optionally, its gradient block) is never built. Pass 1
# adds up its per-feature sums and Gram matrix; per-feature scaling and the
# lambda weights act on the Gram matrix exactly, and its eigenvectors are the
# loadings. Pass 2 projects each window.
#
# Windows are spatial patches of target cells rather than cell ranges: each
# reads its targets and their neighbours through the store view (`[`), so the
# stored cell order does not matter and the store is never rewritten.

#' @rdname reduceData
#' @param spatial_locs for a `neighborhoodPCAParam`: cells x 2 coordinates
#'   with the cell IDs as row names.
#' @param network,agf_network for a `neighborhoodPCAParam`: the spatial
#'   network(s) of the neighbours, as `parquetEdgeStore`s; an undirected one
#'   gives both ends of each edge.
#' @importClassesFrom Giotto neighborhoodPCAParam
#' @importFrom Giotto neighborhood_operators neighborhood_agf
#' @export
setMethod("reduceData",
    signature(x = "parquetExprBase", param = "neighborhoodPCAParam"),
    function(x, param, ..., spatial_locs, network, agf_network = network) {
        prep <- .pe_npca_prepare(x, param, spatial_locs, network, agf_network)
        on.exit(unlink(prep$tmp, recursive = TRUE), add = TRUE)
        .pe_npca_pca(prep, param$npcs)
    }
)

# The neighbour operators, the block weights and the windows (spatial patches
# of target cells, as view positions).
.pe_npca_prepare <- function(x, param, spatial_locs, network, agf_network) {
    locs <- as.matrix(spatial_locs)[x@cell_ids, , drop = FALSE]
    if (anyNA(locs)) {
        stop("[neighborhoodPCAParam] some cells have no spatial location",
            call. = FALSE)
    }
    agf <- isTRUE(param$compute_agf)
    operators <- function(net) {
        if (!inherits(net, "parquetEdgeStore")) {
            stop("[neighborhoodPCAParam] `network` must be a parquetEdgeStore",
                call. = FALSE)
        }
        e <- .edge_view_ends(net, .edge_to_cell_map(net, x))
        # an undirected network (Delaunay) gives each pair to both cells
        if (!isTRUE(net@directed)) {
            e <- list(a = c(e$a, e$b), b = c(e$b, e$a))
        }
        neighborhood_operators(e$a, e$b, locs)
    }
    ops <- operators(network)
    ops1 <- if (agf && !identical(agf_network, network)) operators(agf_network)
    ops <- list(w = ops$w, a = (ops1 %||% ops)$a, wc = (ops1 %||% ops)$wc,
        ws = (ops1 %||% ops)$ws)
    if (!agf) ops <- ops["w"]
    P <- length(x@feat_ids)
    budget <- .pe_npca_window_cells(x, P * (2L + agf))
    lam <- if (agf) {
        c(1 - param$lambda, param$lambda * c(2, 1) / 3)
    } else {
        c(1 - param$lambda, param$lambda)
    }
    list(src = x, ops = ops, tmp = character(), lam = lam,
        absolute = isTRUE(param$absolute),
        windows = .pe_spatial_windows(locs, budget))
}

# Cells grouped into compact spatial patches of at most `budget` cells: a grid
# of about four tiles per patch, walked row by row in alternating direction,
# with consecutive tiles gathered until the next would exceed the budget.
.pe_spatial_windows <- function(locs, budget) {
    n <- nrow(locs)
    if (budget >= n) return(list(seq_len(n)))
    r <- apply(locs[, 1:2, drop = FALSE], 2L, range)
    side <- sqrt(prod(pmax(r[2L, ] - r[1L, ], 1e-12)) * (budget / 4) / n)
    gx <- floor((locs[, 1L] - r[1L, 1L]) / side)
    gy <- floor((locs[, 2L] - r[1L, 2L]) / side)
    gx <- ifelse(gy %% 2 == 0, gx, max(gx) - gx)
    tiles <- split(seq_len(n), gy * (max(gx) + 1) + gx)
    out <- list()
    cur <- integer()
    for (t in tiles) {
        if (length(cur) && length(cur) + length(t) > budget) {
            out[[length(out) + 1L]] <- sort(cur)
            cur <- integer()
        }
        cur <- c(cur, t)
    }
    if (length(cur)) out[[length(out) + 1L]] <- sort(cur)
    out
}

# Target cells per window: the dense blocks of a window (features x 2 or 3,
# with the product temporaries) plus its sparse read. A pinned
# `giottodisk.chunk_size` keeps its meaning of a quarter-window.
.pe_npca_window_cells <- function(x, n_joint) {
    n <- as.numeric(x@n_cells)
    pinned <- getOption("giottodisk.chunk_size")
    if (!is.null(pinned)) return(max(1L, as.integer(pinned) %/% 4L))
    free <- tryCatch(.sc_free_ram(), error = function(e) NA_real_)
    bytes <- if (isTRUE(is.finite(free)) && free > 0) {
        .chunk_budget(free, .ram_frac(), n_cells = 0, k = 0, p = 0)
    } else 2e9
    dens <- tryCatch(.pe_view_density(x), error = function(e) 1)
    nz <- length(x@feat_ids) * min(1, dens)
    per_cell <- 150 * nz + 8 * n_joint * 4
    max(1L, as.integer(min(n, bytes / per_cell)))
}

# Window `i` of the joint values: features x 2 or 3 blocks (rows: x, x_H0,
# x_H1) x the window's target cells, unscaled, and their view positions.
.pe_npca_block <- function(prep, i) {
    tgt <- prep$windows[[i]]
    rows <- lapply(prep$ops, function(o) o[tgt, , drop = FALSE])
    nbr <- unlist(lapply(rows, function(r) which(diff(r@p) > 0L)))
    idx <- sort(unique(c(tgt, nbr)))
    feats <- prep$src@feat_ids
    P <- length(feats)
    out <- matrix(0, P * length(prep$lam), length(tgt))
    X <- .pe_read_cells(prep$src, idx)
    if (is.null(X)) return(list(J = out, vp = tgt))
    op <- function(name) Matrix::t(rows[[name]][, idx, drop = FALSE])
    out[seq_len(P), ] <- as.matrix(X[, match(tgt, idx), drop = FALSE])
    h0 <- as.matrix(X %*% op("w"))
    out[P + seq_len(P), ] <- if (prep$absolute) abs(h0) else h0
    if (length(prep$lam) == 3L) {
        sub <- lapply(rows[c("a", "wc", "ws")], function(r) r[, idx, drop = FALSE])
        out[2L * P + seq_len(P), ] <- neighborhood_agf(X, sub)
    }
    list(J = out, vp = tgt)
}

# Values of the view's cells `idx` (sorted view positions) through the view
# (`[`): features (all of the view's, in order) x cells (in `idx` order);
# NULL when they hold no values.
.pe_read_cells <- function(x, idx) {
    X <- storeRead(x[, idx], output = "dgcmatrix", max_rows = Inf,
        max_cols = Inf)
    if (length(X@x) == 0L) return(NULL)
    Xt <- methods::as(X, "TsparseMatrix")
    ri <- match(rownames(X), x@feat_ids)
    cj <- match(colnames(X), x@cell_ids[idx])
    Matrix::sparseMatrix(i = ri[Xt@i + 1L], j = cj[Xt@j + 1L], x = Xt@x,
        dims = c(length(x@feat_ids), length(idx)),
        dimnames = list(x@feat_ids, NULL))
}

# Both passes and the eigen solve; the output of pcaParam's reduceData()
# methods (u: cells x npcs coordinates in view order, d, v, sdev,
# eigenvalues).
.pe_npca_pca <- function(prep, npcs) {
    n <- as.numeric(prep$src@n_cells)
    feats <- prep$src@feat_ids
    P <- length(feats)
    nb <- length(prep$lam)
    G <- 0
    s <- 0
    for (i in seq_along(prep$windows)) {
        J <- .pe_npca_block(prep, i)$J
        G <- G + tcrossprod(J)
        s <- s + rowSums(J)
    }
    mu <- s / n
    C <- G - n * tcrossprod(mu)
    sd <- sqrt(pmax(diag(C), 0) / (n - 1))
    keep <- sd > sqrt(.Machine$double.eps) * pmax(1, abs(mu))
    cf <- rep(sqrt(prep$lam), each = P) / sd
    cf[!keep] <- 0
    e <- eigen(C[keep, keep] * tcrossprod(cf[keep]), symmetric = TRUE)
    npcs <- min(npcs, sum(keep))
    V <- e$vectors[, seq_len(npcs), drop = FALSE]
    V <- sweep(V, 2L, sign(V[cbind(max.col(abs(t(V))), seq_len(npcs))]), "*")
    rownames(V) <- paste0(rep(feats, nb),
        rep(c("", "_H0", "_H1")[seq_len(nb)], each = P))[keep]
    colnames(V) <- paste0("Dim.", seq_len(npcs))
    u <- matrix(0, n, npcs, dimnames = list(NULL, colnames(V)))
    for (i in seq_along(prep$windows)) {
        b <- .pe_npca_block(prep, i)
        z <- (b$J[keep, , drop = FALSE] - mu[keep]) * cf[keep]
        u[b$vp, ] <- crossprod(z, V)
    }
    rownames(u) <- prep$src@cell_ids
    ev <- pmax(e$values[seq_len(npcs)], 0)
    list(u = u, d = sqrt(ev), v = V, sdev = sqrt(ev / (n - 1)),
        eigenvalues = ev / (n - 1))
}
