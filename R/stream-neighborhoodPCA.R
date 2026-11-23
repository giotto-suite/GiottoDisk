# Giotto's neighborhoodPCAParam (BANKSY embedding) on expression stores, in
# two passes over windows of target cells: the joint matrix (each feature, its
# neighbour mean and, optionally, its gradient block) is never built. Pass 1
# adds up its per-feature sums and Gram matrix; per-feature scaling and the
# lambda weights act on the Gram matrix exactly, and its eigenvectors are the
# loadings. Pass 2 projects each window.

#' @rdname reduceData
#' @param spatial_locs for a `neighborhoodPCAParam`: cells x 2 coordinates
#'   with the cell IDs as row names.
#' @importClassesFrom Giotto neighborhoodPCAParam
#' @export
setMethod("reduceData",
    signature(x = "parquetExprBase", param = "neighborhoodPCAParam"),
    function(x, param, ..., spatial_locs) {
        prep <- .pe_npca_prepare(x, param, spatial_locs)
        on.exit(unlink(prep$tmp, recursive = TRUE), add = TRUE)
        .pe_npca_pca(prep, param$npcs)
    }
)

# The materialized source in a local cell order, the neighbour operators in
# that order, the block weights and the windows of target cells.
.pe_npca_prepare <- function(x, param, spatial_locs) {
    n_cells <- as.integer(x@n_cells)
    locs <- as.matrix(spatial_locs)[x@cell_ids, , drop = FALSE]
    if (anyNA(locs)) {
        stop("[neighborhoodPCAParam] some cells have no spatial location",
            call. = FALSE)
    }
    k0 <- param$k_geom[1L]
    k1 <- param$k_geom[length(param$k_geom)]
    agf <- isTRUE(param$compute_agf)
    ops <- Giotto:::.neighbor_operators(locs, k0)
    ops1 <- if (agf && k1 != k0) Giotto:::.neighbor_operators(locs, k1)
    ops <- list(w = ops$w, a = (ops1 %||% ops)$a, wc = (ops1 %||% ops)$wc,
        ws = (ops1 %||% ops)$ws)
    if (!agf) ops <- ops["w"]
    nb <- Matrix::summary(Reduce(`+`, lapply(ops, abs)))
    pairs <- list(a = nb$i, b = nb$j)
    P <- length(x@feat_ids)
    budget <- .pe_npca_window_cells(x, P * (2L + agf))
    local <- budget >= n_cells || .pairs_local(pairs, budget)
    ord <- if (local) seq_len(n_cells) else .pairs_bfs_order(pairs, n_cells)
    rm(nb, pairs)
    if (!local) ops <- lapply(ops, function(o) o[ord, ord])
    tmp <- .dump_tempfile()
    # a single store in a local order is read in place; otherwise copied
    subs <- .exprbase_substores(x)
    src <- if (local && length(subs) == 1L) {
        x
    } else if (local) {
        storeWrite(parquetExprStore(path = tmp), x)
    } else {
        .pe_rewrite_cell_order(x, ord, tmp)
    }
    lam <- if (agf) {
        c(1 - param$lambda, param$lambda * c(2, 1) / 3)
    } else {
        c(1 - param$lambda, param$lambda)
    }
    if (!identical(src, x)) subs <- .exprbase_substores(src)
    list(src = src, info = list(sub = subs[[1L]]$store),
        ops = ops, ord = ord, tmp = tmp, lam = lam,
        absolute = isTRUE(param$absolute),
        windows = .pe_chunk_ranges(1L, n_cells, budget))
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
    w <- prep$windows[[i]]
    tgt <- w[[1L]]:w[[2L]]
    rows <- lapply(prep$ops, function(o) o[tgt, , drop = FALSE])
    nbr <- unlist(lapply(rows, function(r) which(diff(r@p) > 0L)))
    rng <- range(c(tgt, nbr))
    X <- .pe_read_cell_range(prep, rng)
    P <- length(prep$src@feat_ids)
    out <- matrix(0, P * length(prep$lam), length(tgt))
    if (is.null(X)) return(list(J = out, vp = prep$ord[tgt]))
    X <- X[prep$src@feat_ids, , drop = FALSE]
    rr <- rng[1L]:rng[2L]
    op <- function(name) Matrix::t(rows[[name]][, rr, drop = FALSE])
    out[seq_len(P), ] <- as.matrix(X[, tgt - rng[1L] + 1L, drop = FALSE])
    h0 <- as.matrix(X %*% op("w"))
    out[P + seq_len(P), ] <- if (prep$absolute) abs(h0) else h0
    if (length(prep$lam) == 3L) {
        m <- as.matrix(X %*% op("a"))
        re <- as.matrix(X %*% op("wc")) -
            sweep(m, 2L, Matrix::rowSums(rows$wc), "*")
        im <- as.matrix(X %*% op("ws")) -
            sweep(m, 2L, Matrix::rowSums(rows$ws), "*")
        out[2L * P + seq_len(P), ] <- sqrt(re^2 + im^2)
    }
    list(J = out, vp = prep$ord[tgt])
}

# Values of the source cells rng[1]:rng[2] (features x cells, all of the
# range's columns, in order); NULL when the range holds no values.
.pe_read_cell_range <- function(prep, rng) {
    X <- .pe_read_chunk_sub(prep$info, rng[1L], rng[2L], list(),
        length(prep$src@feat_ids))
    if (is.null(X)) return(NULL)
    Xt <- methods::as(X, "TsparseMatrix")
    cj <- match(colnames(X), prep$src@cell_ids[rng[1L]:rng[2L]])
    full <- Matrix::sparseMatrix(i = Xt@i + 1L, j = cj[Xt@j + 1L], x = Xt@x,
        dims = c(nrow(X), rng[2L] - rng[1L] + 1L),
        dimnames = list(rownames(X), NULL))
    miss <- setdiff(prep$src@feat_ids, rownames(full))
    if (length(miss)) {
        full <- rbind(full, Matrix::sparseMatrix(i = integer(), j = integer(),
            dims = c(length(miss), ncol(full)),
            dimnames = list(miss, NULL)))
    }
    full
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
    rownames(u) <- prep$src@cell_ids[order(prep$ord)]
    ev <- pmax(e$values[seq_len(npcs)], 0)
    list(u = u, d = sqrt(ev), v = V, sdev = sqrt(ev / (n - 1)),
        eigenvalues = ev / (n - 1))
}
