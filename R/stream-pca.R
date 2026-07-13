#' @include class-parquetExprStore.R
#' @include pca-param.R
NULL

# stream-pca ####
# Streaming randomized SVD (Halko, Martinsson & Tropp 2011) with streaming
# Cholesky-QR for parquetExprStore-backed expression. Plugs into
# GiottoClass's reduceData(x, randomPcaParam) dispatch via:
#
#   reduceData(parquetExprStore, randomPcaParam)
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
# Centering is implicit: column means of normalized data are stored on the
# store via the HVG step (pe@params$norm_means_hvg) and subtracted
# analytically inside the forward/backward passes — no densification.
#
# irlbaPcaParam / exactPcaParam are NOT supported on parquetExprStore;
# they require Lanczos-style iteration on the full sparse matrix and have
# no streaming advantage.

# ---- randomPcaParam: streaming Halko ---------------------------------------

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprStore", param = "randomPcaParam"),
    function(x, param, ...) {
        .stream_random_svd(
            pe            = x,
            k             = param$ncp,
            n_oversamples = param$n_oversamples,
            n_power_iter  = param$n_power_iter,
            feats_to_use  = param$feats_to_use,
            center        = isTRUE(param$center),
            scale         = isTRUE(param$scale),
            set_seed      = isTRUE(param$set_seed),
            seed_number   = param$seed_number
        )
    }
)


# ---- Other pcaParam variants on parquet: clear error ----------------------

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprStore", param = "irlbaPcaParam"),
    function(x, param, ...) {
        stop("[reduceData(parquetExprStore, irlbaPcaParam)] ",
             "method = \"irlba\" is not supported for streaming. ",
             "Use method = \"random\" (Halko randomized SVD) instead.",
             call. = FALSE)
    }
)

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprStore", param = "exactPcaParam"),
    function(x, param, ...) {
        stop("[reduceData(parquetExprStore, exactPcaParam)] ",
             "method = \"exact\" is not supported for streaming. ",
             "Use method = \"random\" (Halko randomized SVD) instead.",
             call. = FALSE)
    }
)


# gramEigenPcaParam: streaming exact PCA via AᵀA. See .stream_gram_svd
# for the algorithm; delegates to Halko when the condition-squared
# blowup exceeds param$fallback_relerr.

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprStore", param = "gramEigenPcaParam"),
    function(x, param, ...) {
        .stream_gram_svd(
            pe               = x,
            k                = param$ncp,
            feats_to_use     = param$feats_to_use,
            center           = isTRUE(param$center),
            scale            = isTRUE(param$scale),
            fallback_relerr  = param$fallback_relerr,
            set_seed         = isTRUE(param$set_seed),
            seed_number      = param$seed_number,
            n_oversamples    = param$n_oversamples,
            n_power_iter     = param$n_power_iter
        )
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
        resolved <- if (gram_bytes < budget_gb * 1e9) {
            gramEigenPcaParam(
                ncp             = param$ncp,
                center          = param$center,
                scale           = param$scale,
                feats_to_use    = param$feats_to_use,
                n_oversamples   = param$n_oversamples,
                n_power_iter    = param$n_power_iter,
                set_seed        = param$set_seed,
                seed_number     = param$seed_number
                # fallback_relerr uses factory default; users pin
                # explicitly via gramEigenPcaParam(fallback_relerr = ...)
            )
        } else {
            Giotto::pcaParam(
                method        = "random",
                ncp           = param$ncp,
                center        = param$center,
                scale         = param$scale,
                feats_to_use  = param$feats_to_use,
                n_oversamples = param$n_oversamples,
                n_power_iter  = param$n_power_iter,
                set_seed      = param$set_seed,
                seed_number   = param$seed_number
            )
        }
        if (isTRUE(param$dry_run)) return(resolved)
        reduceData(x, resolved, ...)
    }
)


# ---- Streaming Halko core --------------------------------------------------

.stream_random_svd <- function(pe, k, n_oversamples = 10L, n_power_iter = 2L,
                                feats_to_use = NULL, center = TRUE,
                                scale = FALSE,
                                set_seed = TRUE, seed_number = 1234L) {
    if (set_seed) set.seed(seed_number)

    norm     <- pe@params$norm
    has_norm <- !is.null(norm) && !is.null(norm$scale_factors)
    scalef   <- if (has_norm) norm$scale_factors else NULL
    log_norm <- isTRUE(norm$log)
    log_base <- norm$base %null% 2

    n_cells <- as.integer(pe@n_cells)
    chunk_size <- as.integer(pe@chunk_size %null% 250000L)

    # feats_to_use = NULL -> use every feature.
    hvg_idx <- if (is.null(feats_to_use)) {
        seq_along(pe@feat_ids)
    } else {
        idx <- match(feats_to_use, pe@feat_ids)
        if (anyNA(idx)) {
            bad <- feats_to_use[is.na(idx)]
            stop("[stream PCA] feats_to_use has IDs not in pe@feat_ids: ",
                 toString(head(bad, 5L)), call. = FALSE)
        }
        idx
    }
    P_hvg <- length(hvg_idx)
    k     <- as.integer(k)
    if (k >= P_hvg) {
        warning("[stream PCA] ncp (", k, ") >= n_HVG (", P_hvg,
                "), setting ncp = ", P_hvg - 1L, call. = FALSE)
        k <- P_hvg - 1L
    }
    k_total <- k + as.integer(n_oversamples)

    # Streaming stats: means (center) + sds (scale) in one pass.
    stats <- .stream_norm_hvg_stats(pe, hvg_idx)
    means <- if (center) stats$means else numeric(P_hvg)
    sds   <- if (scale)  stats$sds   else rep(1, P_hvg)

    hvg_orig <- .pe_orig_col(hvg_idx, pe)

    # Chunk reader. `has_norm` FALSE returns values as-written.
    .read_chunk_norm_hvg <- function(cell_start, cell_end) {
        ds <- pe@read_fun(pe@path)
        row_id <- col_id <- value <- NULL  # NSE
        orig_rows <- .pe_orig_row(cell_start:cell_end, pe)
        df <- ds |>
            dplyr::filter(row_id %in% !!orig_rows,
                           col_id %in% !!hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) return(NULL)
        chunk_n  <- cell_end - cell_start + 1L
        gene_map <- match(df$col_id, hvg_orig)
        cell_map <- match(df$row_id, orig_rows)
        A <- Matrix::sparseMatrix(
            i = cell_map, j = gene_map, x = as.double(df$value),
            dims = c(chunk_n, P_hvg), repr = "C"
        )
        if (has_norm) {
            scalef_chunk <- scalef[cell_start:cell_end]
            A@x <- A@x * scalef_chunk[A@i + 1L]
            if (log_norm) A@x <- log1p(A@x) / log(log_base)
        }
        A
    }

    # Forward: Y = ((A - 1μᵀ) · diag(1/σ)) · M. σ absorbed into M via
    # row-wise M/σ (sds = 1 when scale=FALSE, so no-op).
    .forward <- function(M) {
        m <- ncol(M)
        M_use <- M / sds
        correction <- if (center) as.numeric(means %*% M_use)
                       else numeric(m)
        Y <- matrix(0.0, nrow = n_cells, ncol = m)
        cs <- 1L
        while (cs <= n_cells) {
            ce <- min(cs + chunk_size - 1L, n_cells)
            A  <- .read_chunk_norm_hvg(cs, ce)
            chunk_n <- ce - cs + 1L
            if (!is.null(A)) {
                Yc <- as.matrix(A %*% M_use)
                if (center) {
                    Yc <- Yc - matrix(correction, nrow = chunk_n,
                                      ncol = m, byrow = TRUE)
                }
                Y[cs:ce, ] <- Yc
            } else if (center) {
                Y[cs:ce, ] <- -matrix(correction, nrow = chunk_n,
                                       ncol = m, byrow = TRUE)
            }
            cs <- ce + 1L
        }
        Y
    }

    # Backward: Z = ((A - 1μᵀ) · diag(1/σ))ᵀ · Y. Accumulate unscaled Z,
    # divide rows by σ once at the end.
    .backward <- function(Y_mat) {
        m <- ncol(Y_mat)
        Z <- matrix(0.0, nrow = P_hvg, ncol = m)
        G <- matrix(0.0, nrow = m,     ncol = m)
        cs_Y <- numeric(m)
        cs <- 1L
        while (cs <= n_cells) {
            ce <- min(cs + chunk_size - 1L, n_cells)
            A  <- .read_chunk_norm_hvg(cs, ce)
            chunk_n <- ce - cs + 1L
            Yc <- Y_mat[cs:ce, , drop = FALSE]
            G  <- G + crossprod(Yc)
            cs_Y <- cs_Y + colSums(Yc)
            if (!is.null(A)) {
                Z <- Z + as.matrix(Matrix::crossprod(A, Yc))
            }
            cs <- ce + 1L
        }
        if (center) Z <- Z - tcrossprod(means, cs_Y)
        Z <- Z / sds
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


# Streaming gram-eigen core. Three passes: (1) per-gene μ + σ,
# (2) accumulate G_raw = Σ chunkᵀchunk with G = G_raw − n·μμᵀ (and
# optional / (σσᵀ) for scale=TRUE), (3) coords = A_c · V. AᵀA squares
# κ(A); when pred rel-err `ε·(d1/dk)²/2 > fallback_relerr`, delegate
# to Halko.

.stream_gram_svd <- function(pe, k, feats_to_use = NULL, center = TRUE,
                              scale = FALSE,
                              fallback_relerr = 0.01,
                              set_seed = TRUE, seed_number = 1234L,
                              n_oversamples = 10L, n_power_iter = 2L) {
    norm     <- pe@params$norm
    has_norm <- !is.null(norm) && !is.null(norm$scale_factors)
    scalef   <- if (has_norm) norm$scale_factors else NULL
    log_norm <- isTRUE(norm$log)
    log_base <- norm$base %null% 2

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
    k     <- as.integer(k)
    if (k >= P_hvg) {
        warning("[stream gram PCA] ncp (", k, ") >= n_HVG (", P_hvg,
                "), setting ncp = ", P_hvg - 1L, call. = FALSE)
        k <- P_hvg - 1L
    }
    hvg_orig <- .pe_orig_col(hvg_idx, pe)

    # Pass 1: means (center) + sds (scale) in one pass.
    stats <- .stream_norm_hvg_stats(pe, hvg_idx)
    means <- if (center) stats$means else numeric(P_hvg)
    sds   <- if (scale)  stats$sds   else rep(1, P_hvg)

    # Chunk reader. `has_norm` FALSE returns values as-written.
    .read_chunk_norm_hvg <- function(cell_start, cell_end) {
        ds <- pe@read_fun(pe@path)
        row_id <- col_id <- value <- NULL  # NSE
        orig_rows <- .pe_orig_row(cell_start:cell_end, pe)
        df <- ds |>
            dplyr::filter(row_id %in% !!orig_rows,
                           col_id %in% !!hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) return(NULL)
        chunk_n  <- cell_end - cell_start + 1L
        gene_map <- match(df$col_id, hvg_orig)
        cell_map <- match(df$row_id, orig_rows)
        A <- Matrix::sparseMatrix(
            i = cell_map, j = gene_map, x = as.double(df$value),
            dims = c(chunk_n, P_hvg), repr = "C"
        )
        if (has_norm) {
            scalef_chunk <- scalef[cell_start:cell_end]
            A@x <- A@x * scalef_chunk[A@i + 1L]
            if (log_norm) A@x <- log1p(A@x) / log(log_base)
        }
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
    ncp_used <- min(k, length(ei$values))
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
            pe = pe, k = k,
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


# Per-HVG-gene mean + sd of normalized data in one streaming pass.
# Variance via SS identity: σ² = (SS − n·μ²)/(n−1). `norm` NULL reads
# values as-written; `hvg_idx` NULL uses every column.
.stream_norm_hvg_stats <- function(pe, hvg_idx) {
    norm     <- pe@params$norm
    has_norm <- !is.null(norm) && !is.null(norm$scale_factors)
    scalef   <- if (has_norm) norm$scale_factors else NULL
    log_norm <- isTRUE(norm$log)
    log_base <- norm$base %null% 2
    n_cells  <- as.integer(pe@n_cells)

    if (is.null(hvg_idx)) hvg_idx <- seq_along(pe@feat_ids)
    P_hvg    <- length(hvg_idx)

    g_sum  <- numeric(P_hvg)
    g_ssum <- numeric(P_hvg)
    row_id <- col_id <- value <- v_norm <- s <- ss <- NULL  # NSE

    hvg_orig <- .pe_orig_col(hvg_idx, pe)

    ds_base <- pe@read_fun(pe@path)
    chunk_size <- as.integer(pe@chunk_size %null% 250000L)
    for (cs in seq.int(1L, n_cells, by = chunk_size)) {
        ce <- min(cs + chunk_size - 1L, n_cells)
        orig_rows <- .pe_orig_row(cs:ce, pe)
        df <- ds_base |>
            dplyr::filter(row_id %in% !!orig_rows,
                           col_id %in% !!hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) next
        if (has_norm) {
            df[, row_id := .pe_remap_row(row_id, pe)]
            df[, v_norm := value * scalef[row_id]]
            if (log_norm) df[, v_norm := log1p(v_norm) / log(log_base)]
        } else {
            df[, v_norm := value]
        }
        agg <- df[, .(s = sum(v_norm), ss = sum(v_norm * v_norm)),
                  by = col_id]
        g_idx <- match(agg$col_id, hvg_orig)
        keep <- !is.na(g_idx)
        g_sum[g_idx[keep]]  <- g_sum[g_idx[keep]]  + agg$s[keep]
        g_ssum[g_idx[keep]] <- g_ssum[g_idx[keep]] + agg$ss[keep]
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
