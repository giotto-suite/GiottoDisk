#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-pca ####
# Streaming randomized SVD (Halko, Martinsson & Tropp 2011) with streaming
# Cholesky-QR for parquetExprBase-backed expression (single or union).
# Plugs into GiottoClass's reduceData(x, randomPcaParam) dispatch via:
#
#   reduceData(parquetExprBase, randomPcaParam)
#       -> list(u, d, v, sdev, eigenvalues)
#
# Algorithm (mirrors scstream::sc_pca):
#   omega ~ Gaussian(P_hvg x (k + p))
#   Y     <- (A_norm - 1*means^T) * omega          # forward pass
#   for q power iterations:
#       Z <- A_norm^T * Y                           # backward
#       Q <- chol-QR(Z)                             # streaming
#       Y <- (A_norm - 1*means^T) * Q               # forward
#   Z      <- A_norm^T * Y                          # final backward
#   B      <- chol-QR(Z) -> Q^T A
#   svd(B) -> recover U, d, V; sign-correct V
#
# Centering is implicit: column means of normalized data are computed by
# `.stream_norm_hvg_means()` in a single streaming pass over the @ops chain
# (summed across substores for union stores), then subtracted analytically
# inside the forward/backward passes — no densification.
#
# Single (`parquetExprStore`) and union (`unionParquetExprStore`) collapse
# to one implementation via `.exprbase_substores()`: forward/backward
# iterate substores, reading per-substore cell chunks, and place results
# into the union-axis Y matrix at the substore's cumulative cell offset.
# For a single store the iterator yields one entry (offset 0) and the
# loop runs once — same behavior as before.
#
# irlbaPcaParam / exactPcaParam are NOT supported on parquetExprBase;
# they require Lanczos-style iteration on the full sparse matrix and have
# no streaming advantage.

# ---- randomPcaParam: streaming Halko ---------------------------------------

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprBase", param = "randomPcaParam"),
    function(x, param, ...) {
        if (!.pe_has_norm_op(x@ops)) {
            stop("[reduceData(parquetExprBase, randomPcaParam)] ",
                 "expression backend has no normalization recipe. Run ",
                 "normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE) ",
                 "first.", call. = FALSE)
        }
        if (isTRUE(param$scale)) {
            stop("[reduceData(parquetExprBase, randomPcaParam)] ",
                 "scale = TRUE (per-gene z-score) is not supported for ",
                 "streaming because it densifies the matrix. Pass ",
                 "scale = FALSE.", call. = FALSE)
        }
        feats <- param$feats_to_use
        if (is.null(feats)) {
            stop("[reduceData(parquetExprBase, randomPcaParam)] ",
                 "feats_to_use is required for the streaming PCA path. ",
                 "Pass the HVG feature IDs (typically rownames where ",
                 "@featMetadata$hvf == \"yes\").", call. = FALSE)
        }

        .stream_random_svd(
            pe            = x,
            k             = param$ncp,
            n_oversamples = param$n_oversamples,
            n_power_iter  = param$n_power_iter,
            feats_to_use  = feats,
            center        = isTRUE(param$center),
            set_seed      = isTRUE(param$set_seed),
            seed_number   = param$seed_number
        )
    }
)


# ---- Other pcaParam variants on parquet: clear error ----------------------

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprBase", param = "irlbaPcaParam"),
    function(x, param, ...) {
        stop("[reduceData(parquetExprBase, irlbaPcaParam)] ",
             "method = \"irlba\" is not supported for streaming. ",
             "Use method = \"random\" (Halko randomized SVD) instead.",
             call. = FALSE)
    }
)

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprBase", param = "exactPcaParam"),
    function(x, param, ...) {
        stop("[reduceData(parquetExprBase, exactPcaParam)] ",
             "method = \"exact\" is not supported for streaming. ",
             "Use method = \"random\" (Halko randomized SVD) instead.",
             call. = FALSE)
    }
)


# gramEigenPcaParam: streaming exact PCA via AᵀA. See .stream_gram_svd
# for the algorithm; delegates to Halko when the condition-squared
# blowup exceeds param$fallback_relerr. Currently only single-store
# (parquetExprStore) — union dispatch is a follow-up.

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprStore", param = "gramEigenPcaParam"),
    function(x, param, ...) {
        args <- as.list(param@param)
        args$method <- NULL
        do.call(.stream_gram_svd, c(list(pe = x), args))
    }
)


# autoPcaParam on parquetExprStore: gram-eigen if P²·8 fits the budget
# (default 5 GB → P ≤ ~25k, covers all realistic scRNA / imaging /
# multi-omic RNA regimes), Halko otherwise. Fixed budget rather than
# %-of-RAM: cache/TLB effects make Halko faster before large machines'
# RAM runs out. Override via option("giottodisk.pca_auto_budget_gb").
# dry_run = TRUE returns the resolved param without running PCA.

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprStore", param = "autoPcaParam"),
    function(x, param, ...) {
        # P estimate: feats_to_use if set (HVG-restricted), else full store
        n_feat <- if (!is.null(param$feats_to_use)) {
            length(param$feats_to_use)
        } else {
            length(x@feat_ids)
        }
        budget_gb  <- getOption("giottodisk.pca_auto_budget_gb", 5)
        gram_bytes <- as.numeric(n_feat)^2 * 8
        knobs <- as.list(param@param)
        knobs$method  <- NULL   # strip auto sentinel
        knobs$dry_run <- NULL   # not a concrete-flavor arg
        resolved <- if (gram_bytes < budget_gb * 1e9) {
            do.call(gramEigenPcaParam, knobs)
        } else {
            do.call(Giotto::pcaParam,
                c(list(method = "random"), knobs))
        }
        if (isTRUE(param$dry_run)) return(resolved)
        reduceData(x, resolved, ...)
    }
)


# ---- Streaming Halko core --------------------------------------------------

.stream_random_svd <- function(pe, k, n_oversamples = 10L, n_power_iter = 2L,
                                feats_to_use, center = TRUE,
                                set_seed = TRUE, seed_number = 1234L) {
    if (set_seed) set.seed(seed_number)

    n_cells <- as.integer(pe@n_cells)
    # chunk_size lives on parquetExprStore; the union doesn't carry one,
    # so fall back to the first substore's value (or a sane default).
    chunk_size <- as.integer(.exprbase_chunk_size(pe))

    # Map HVG feature IDs to integer col_ids on the union/feat axis.
    # feat_ids align across substores (union invariant), so this lookup
    # is unambiguous.
    hvg_idx <- match(feats_to_use, pe@feat_ids)
    if (anyNA(hvg_idx)) {
        bad <- feats_to_use[is.na(hvg_idx)]
        stop("[stream PCA] feats_to_use has IDs not in pe@feat_ids: ",
             toString(head(bad, 5L)), call. = FALSE)
    }
    P_hvg <- length(hvg_idx)
    k     <- as.integer(k)
    if (k >= P_hvg) {
        warning("[stream PCA] ncp (", k, ") >= n_HVG (", P_hvg,
                "), setting ncp = ", P_hvg - 1L, call. = FALSE)
        k <- P_hvg - 1L
    }
    k_total <- k + as.integer(n_oversamples)

    # Build per-substore record list: each entry carries the substore
    # (with parent ops projected on union path), its cumulative cell
    # offset into the union axis, n_sub, and the HVG-as-original col_id
    # vector for that substore (handles per-substore @gene_idx via
    # `.pe_orig_col`). For single store the list has one entry.
    parent_ops <- if (inherits(pe, "unionParquetExprStore")) pe@ops else list()
    sub_infos <- lapply(.exprbase_substores(pe), function(se) {
        sub <- .exprbase_inject_parent_ops(se$store, parent_ops)
        list(sub      = sub,
             offset   = as.integer(se$cell_offset),
             n_sub    = as.integer(sub@n_cells),
             hvg_orig = .pe_orig_col(hvg_idx, sub))
    })

    # ---- Compute per-HVG-gene normalized means (one streaming pass) -------
    means <- if (center) {
        .stream_norm_hvg_means(pe, hvg_idx, sub_infos)
    } else {
        numeric(P_hvg)
    }

    # ---- Per-substore chunk reader (cell-major within substore) ----------
    .read_chunk_sub <- function(info, sub_cs, sub_ce) {
        row_id <- col_id <- NULL  # NSE
        sub <- info$sub
        orig_rows <- .pe_orig_row(sub_cs:sub_ce, sub)
        df <- storeRead(sub, output = "query") |>
            dplyr::filter(row_id %in% !!orig_rows,
                           col_id %in% !!info$hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) return(NULL)
        chunk_n <- sub_ce - sub_cs + 1L
        gene_map <- match(df$col_id, info$hvg_orig)
        i_within <- match(df$row_id, orig_rows)
        Matrix::sparseMatrix(
            i = i_within, j = gene_map, x = as.double(df$v_norm),
            dims = c(chunk_n, P_hvg), repr = "C"
        )
    }

    # ---- Forward: Y = (A_norm - 1·means^T) · M  --------------------------
    # Y is sized to the UNION cell axis; per-substore reads fill the
    # appropriate row band at `offset + sub_cs:sub_ce`.
    .forward <- function(M) {
        m <- ncol(M)
        correction <- if (center) as.numeric(means %*% M)
                       else numeric(m)
        Y <- matrix(0.0, nrow = n_cells, ncol = m)
        for (info in sub_infos) {
            offset <- info$offset
            n_sub  <- info$n_sub
            cs <- 1L
            while (cs <= n_sub) {
                ce <- min(cs + chunk_size - 1L, n_sub)
                A  <- .read_chunk_sub(info, cs, ce)
                chunk_n <- ce - cs + 1L
                rows <- (offset + cs):(offset + ce)
                if (!is.null(A)) {
                    Yc <- as.matrix(A %*% M)
                    if (center) {
                        Yc <- Yc - matrix(correction, nrow = chunk_n,
                                          ncol = m, byrow = TRUE)
                    }
                    Y[rows, ] <- Yc
                } else if (center) {
                    Y[rows, ] <- -matrix(correction, nrow = chunk_n,
                                          ncol = m, byrow = TRUE)
                }
                cs <- ce + 1L
            }
        }
        Y
    }

    # ---- Backward: returns Z = A_norm^T · Y  +  Gram G = Y^T Y -----------
    .backward <- function(Y_mat) {
        m <- ncol(Y_mat)
        Z <- matrix(0.0, nrow = P_hvg, ncol = m)
        G <- matrix(0.0, nrow = m,     ncol = m)
        cs_Y <- numeric(m)
        for (info in sub_infos) {
            offset <- info$offset
            n_sub  <- info$n_sub
            cs <- 1L
            while (cs <= n_sub) {
                ce <- min(cs + chunk_size - 1L, n_sub)
                A  <- .read_chunk_sub(info, cs, ce)
                rows <- (offset + cs):(offset + ce)
                Yc <- Y_mat[rows, , drop = FALSE]
                G  <- G + crossprod(Yc)
                cs_Y <- cs_Y + colSums(Yc)
                if (!is.null(A)) {
                    Z <- Z + as.matrix(Matrix::crossprod(A, Yc))
                }
                cs <- ce + 1L
            }
        }
        if (center) Z <- Z - tcrossprod(means, cs_Y)  # implicit centering
        list(Z = Z, G = G)
    }

    # ---- Halko algorithm -------------------------------------------------
    omega <- matrix(stats::rnorm(P_hvg * k_total), nrow = P_hvg, ncol = k_total)

    Y <- .forward(omega)
    for (i in seq_len(n_power_iter)) {
        zg <- .backward(Y)
        R_chol <- chol(zg$G)
        Z_orth <- t(backsolve(R_chol, t(zg$Z), transpose = TRUE))
        Q_z    <- qr.Q(qr(Z_orth))
        Y      <- .forward(Q_z)
    }

    zg_final <- .backward(Y)
    R_chol   <- chol(zg_final$G)
    B        <- backsolve(R_chol, t(zg_final$Z), transpose = TRUE)

    sv <- svd(B, nu = k_total, nv = k_total)

    M_recover <- backsolve(R_chol, sv$u[, seq_len(k), drop = FALSE])
    U <- Y %*% M_recover
    D_k <- sv$d[seq_len(k)]
    V   <- sv$v[, seq_len(k), drop = FALSE]

    # Sign convention: largest |V[, j]| entry positive
    signs <- vapply(seq_len(k), function(j) {
        sign(V[which.max(abs(V[, j])), j])
    }, numeric(1L))
    signs[signs == 0] <- 1
    V <- sweep(V, 2L, signs, "*")
    U <- sweep(U, 2L, signs, "*")

    rownames(U) <- pe@cell_ids
    rownames(V) <- pe@feat_ids[hvg_idx]

    eigenvalues <- D_k^2 / (n_cells - 1L)
    list(
        u           = sweep(U, 2L, D_k, "*"),    # cells × k coords (u · d)
        d           = D_k,
        v           = V,
        sdev        = sqrt(eigenvalues),
        eigenvalues = eigenvalues
    )
}


# Helper: per-HVG-gene mean of normalized data (one streaming pass per
# substore, summed across substores for union stores). Used for implicit
# centering inside .forward / .backward.
.stream_norm_hvg_means <- function(pe, hvg_idx, sub_infos = NULL) {
    n_cells  <- as.integer(pe@n_cells)
    P_hvg    <- length(hvg_idx)

    g_sum <- numeric(P_hvg)
    col_id <- v_norm <- s <- NULL  # NSE

    if (is.null(sub_infos)) {
        parent_ops <- if (inherits(pe, "unionParquetExprStore")) {
            pe@ops
        } else {
            list()
        }
        sub_infos <- lapply(.exprbase_substores(pe), function(se) {
            sub <- .exprbase_inject_parent_ops(se$store, parent_ops)
            list(sub = sub, hvg_orig = .pe_orig_col(hvg_idx, sub))
        })
    }

    for (info in sub_infos) {
        agg <- storeRead(info$sub, output = "query") |>
            dplyr::filter(col_id %in% !!info$hvg_orig) |>
            dplyr::group_by(col_id) |>
            dplyr::summarise(s = sum(v_norm, na.rm = TRUE)) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(agg) > 0L) {
            g_idx <- match(agg$col_id, info$hvg_orig)
            keep <- !is.na(g_idx)
            g_sum[g_idx[keep]] <- g_sum[g_idx[keep]] + as.numeric(agg$s[keep])
        }
    }
    g_sum / n_cells
}


# Helper: pick a chunk_size for the streaming chunk reader. parquetExprStore
# carries @chunk_size directly; for a union store, defer to the first
# substore's value (substores can have different chunk sizes in principle,
# but the union-axis chunking is uniform so we use one — the first is a
# safe default given the constructor enforces compatible substores).
.exprbase_chunk_size <- function(pe) {
    if (inherits(pe, "parquetExprStore")) {
        return(pe@chunk_size %null% 250000L)
    }
    if (inherits(pe, "unionParquetExprStore") && length(pe@stores) > 0L) {
        return(pe@stores[[1L]]@chunk_size %null% 250000L)
    }
    250000L
}
# Streaming gram-eigen core. Three passes: (1) per-gene μ + σ,
# (2) accumulate G_raw = Σ chunkᵀchunk with G = G_raw − n·μμᵀ (and
# optional / (σσᵀ) for scale=TRUE), (3) coords = A_c · V. AᵀA squares
# κ(A); when pred rel-err `ε·(d1/dk)²/2 > fallback_relerr`, delegate
# to Halko.

.stream_gram_svd <- function(pe, ncp, feats_to_use = NULL, center = TRUE,
                              scale = FALSE,
                              fallback_relerr = 0.01,
                              set_seed = TRUE, seed_number = 1234L,
                              n_oversamples = 10L, n_power_iter = 2L, ...) {
    n_cells    <- as.integer(pe@n_cells)
    chunk_size <- as.integer(pe@chunk_size %null% 250000L)

    # feats_to_use = NULL -> use every feature.
    hvg_idx <- if (is.null(feats_to_use)) {
        seq_along(pe@feat_ids)
    } else {
        idx <- match(feats_to_use, pe@feat_ids)
        if (anyNA(idx)) {
            bad <- feats_to_use[is.na(idx)]
            stop("[stream gram PCA] feats_to_use has IDs not in pe@feat_ids: ",
                 toString(head(bad, 5L)), call. = FALSE)
        }
        idx
    }
    P_hvg <- length(hvg_idx)
    ncp   <- as.integer(ncp)
    if (ncp >= P_hvg) {
        warning("[stream gram PCA] ncp (", ncp, ") >= n_HVG (", P_hvg,
                "), setting ncp = ", P_hvg - 1L, call. = FALSE)
        ncp <- P_hvg - 1L
    }
    hvg_orig <- .pe_orig_col(hvg_idx, pe)

    # Pass 1: means (center) + sds (scale) in one pass over @ops.
    stats <- .stream_norm_hvg_stats(pe, hvg_idx)
    means <- if (center) stats$means else numeric(P_hvg)
    sds   <- if (scale)  stats$sds   else rep(1, P_hvg)

    # Chunk reader. storeRead(pe) returns the arrow query with @ops
    # composed in (v_norm already projected by any norm_libsize_log op).
    # When no norm op is queued, uses raw `value`.
    val_col <- if (.pe_has_norm_op(pe@ops)) "v_norm" else "value"
    .read_chunk_norm_hvg <- function(cell_start, cell_end) {
        row_id <- col_id <- NULL  # NSE
        orig_rows <- .pe_orig_row(cell_start:cell_end, pe)
        df <- storeRead(pe, output = "query") |>
            dplyr::filter(row_id %in% !!orig_rows,
                           col_id %in% !!hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) return(NULL)
        chunk_n  <- cell_end - cell_start + 1L
        gene_map <- match(df$col_id, hvg_orig)
        cell_map <- match(df$row_id, orig_rows)
        A <- Matrix::sparseMatrix(
            i = cell_map, j = gene_map, x = as.double(df[[val_col]]),
            dims = c(chunk_n, P_hvg), repr = "C"
        )
        A
    }

    # Pass 2: G_raw = Σ chunk^T chunk (one streaming crossprod per band)
    G_raw <- matrix(0.0, nrow = P_hvg, ncol = P_hvg)
    cs <- 1L
    while (cs <= n_cells) {
        ce <- min(cs + chunk_size - 1L, n_cells)
        A  <- .read_chunk_norm_hvg(cs, ce)
        if (!is.null(A)) {
            G_raw <- G_raw + as.matrix(Matrix::crossprod(A))
        }
        cs <- ce + 1L
    }
    # Center: (A - 1μᵀ)ᵀ(A - 1μᵀ) = AᵀA - n·μμᵀ. Scale: divide by σσᵀ.
    G <- if (center) G_raw - as.numeric(n_cells) * tcrossprod(means)
         else       G_raw
    if (scale) G <- G / tcrossprod(sds)

    # Eigendecomposition (descending)
    ei <- eigen(G, symmetric = TRUE)
    ncp_used <- min(ncp, length(ei$values))
    eigenvalues_G <- ei$values[seq_len(ncp_used)]
    V <- ei$vectors[, seq_len(ncp_used), drop = FALSE]
    eigenvalues_G[eigenvalues_G < 0] <- 0    # tiny-negative clamp
    d <- sqrt(eigenvalues_G)

    # Fallback trigger: condition-squared blowup on the trailing PCs
    eps <- .Machine$double.eps
    d1 <- d[1L]
    dk <- d[ncp_used]
    pred_relerr <- if (dk > 0) eps * (d1 / dk)^2 / 2 else Inf
    if (pred_relerr > fallback_relerr) {
        warning("[stream gram PCA] predicted rel-err on d_k = ",
            format(pred_relerr, digits = 3),
            " exceeds fallback_relerr (", fallback_relerr,
            "). Delegating to streaming random SVD (Halko).",
            call. = FALSE)
        return(.stream_random_svd(
            pe = pe, ncp = ncp,
            n_oversamples = n_oversamples,
            n_power_iter  = n_power_iter,
            feats_to_use  = feats_to_use, center = center, scale = scale,
            set_seed = set_seed, seed_number = seed_number
        ))
    }

    # Sign convention: largest |V[, j]| entry positive (same as Halko path)
    signs <- vapply(seq_len(ncp_used), function(j) {
        sign(V[which.max(abs(V[, j])), j])
    }, numeric(1L))
    signs[signs == 0] <- 1
    V <- sweep(V, 2L, signs, "*")

    # Pass 3: coords = A_std · V. Absorb σ into V (V/σ, row-wise), so
    # per-chunk work is identical for scale on/off.
    V_use <- V / sds
    coords <- matrix(0.0, nrow = n_cells, ncol = ncp_used)
    correction <- if (center) {
        as.numeric(means %*% V_use)
    } else {
        numeric(ncp_used)
    }
    cs <- 1L
    while (cs <= n_cells) {
        ce <- min(cs + chunk_size - 1L, n_cells)
        A  <- .read_chunk_norm_hvg(cs, ce)
        chunk_n <- ce - cs + 1L
        if (!is.null(A)) {
            Cc <- as.matrix(A %*% V_use)
            if (center) {
                Cc <- Cc - matrix(correction, nrow = chunk_n,
                                   ncol = ncp_used, byrow = TRUE)
            }
            coords[cs:ce, ] <- Cc
        } else if (center) {
            coords[cs:ce, ] <- -matrix(correction, nrow = chunk_n,
                                        ncol = ncp_used, byrow = TRUE)
        }
        cs <- ce + 1L
    }

    rownames(coords) <- pe@cell_ids
    rownames(V)      <- pe@feat_ids[hvg_idx]

    eigenvalues <- d^2 / (n_cells - 1L)
    list(
        u           = coords,   # cells × k coords (= u · d in SVD notation)
        d           = d,
        v           = V,
        sdev        = sqrt(eigenvalues),
        eigenvalues = eigenvalues
    )
}


# Per-HVG-gene mean + sd of normalized data in one arrow-side pass.
# Variance via SS identity: σ² = (SS − n·μ²)/(n−1). `hvg_idx` NULL uses
# every column. Reads `v_norm` when an @ops norm op is queued; otherwise
# aggregates `value` directly.
.stream_norm_hvg_stats <- function(pe, hvg_idx, ...) {
    n_cells  <- as.integer(pe@n_cells)

    if (is.null(hvg_idx)) hvg_idx <- seq_along(pe@feat_ids)
    P_hvg <- length(hvg_idx)

    hvg_orig <- .pe_orig_col(hvg_idx, pe)

    # NSE bindings
    col_id <- value <- v_norm <- s <- ss <- NULL

    # storeRead returns the arrow query with @ops composed in — v_norm is
    # already projected by any norm_libsize_log op. Per-gene aggregation
    # runs at arrow layer; only the (P_hvg, 2)-sized result transfers to R.
    val_expr <- if (.pe_has_norm_op(pe@ops)) rlang::expr(v_norm)
                 else rlang::expr(value)
    agg <- storeRead(pe, output = "query") |>
        dplyr::filter(col_id %in% !!hvg_orig) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(
            s  = sum(!!val_expr, na.rm = TRUE),
            ss = sum((!!val_expr) * (!!val_expr), na.rm = TRUE)
        ) |>
        dplyr::collect() |>
        data.table::as.data.table()

    g_sum  <- numeric(P_hvg)
    g_ssum <- numeric(P_hvg)
    if (nrow(agg) > 0L) {
        g_idx <- match(agg$col_id, hvg_orig)
        keep  <- !is.na(g_idx)
        g_sum[g_idx[keep]]  <- as.numeric(agg$s[keep])
        g_ssum[g_idx[keep]] <- as.numeric(agg$ss[keep])
    }

    means <- g_sum / n_cells
    vars  <- (g_ssum - as.numeric(n_cells) * means^2) /
        max(as.numeric(n_cells) - 1, 1)
    vars[vars < 0] <- 0   # tiny-negative clamp (numerical)
    sds   <- sqrt(vars)
    # Zero-variance genes: keep σ = 1 for safe division downstream
    # (a truly constant feature contributes nothing to any PC anyway).
    sds[sds == 0] <- 1
    list(means = means, sds = sds)
}
