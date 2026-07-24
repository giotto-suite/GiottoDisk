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
# `.stream_norm_hvg_means()` — one collect + R-side @post_ops apply per
# substore (summed across substores for union stores), then subtracted
# analytically inside the forward/backward passes — no densification.
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
        if (!.pe_has_norm_op(x)) {
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
# blowup exceeds param$fallback_relerr. Works on parquetExprBase
# (single store or union) via per-substore chunk iteration.

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprBase", param = "gramEigenPcaParam"),
    function(x, param, ...) {
        args <- as.list(param@param)
        args$method <- NULL
        do.call(.stream_gram_svd, c(list(pe = x), args))
    }
)


# autoPcaParam on parquetExprBase: gram-eigen if P²·8 fits the budget
# (default 5 GB → P ≤ ~25k, covers all realistic scRNA / imaging /
# multi-omic RNA regimes), Halko otherwise. Fixed budget rather than
# %-of-RAM: cache/TLB effects make Halko faster before large machines'
# RAM runs out. Override via option("giottodisk.pca_auto_budget_gb").
# dry_run = TRUE returns the resolved param without running PCA.

#' @rdname reduceData
#' @export
setMethod("reduceData",
    signature(x = "parquetExprBase", param = "autoPcaParam"),
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
    post_ops <- pe@post_ops
    sub_infos <- lapply(.exprbase_substores(pe), function(se) {
        sub <- .exprbase_inject_parent_ops(se$store, parent_ops)
        # Pre-extract per-substore positional scalef vectors for each
        # post op that has per-cell state on this substore.
        scalef_vecs <- .pe_scalef_vecs_for_sub(post_ops, sub@uid)
        list(sub         = sub,
             offset      = as.integer(se$cell_offset),
             n_sub       = as.integer(sub@n_cells),
             hvg_orig    = .pe_orig_col(hvg_idx, sub),
             scalef_vecs = scalef_vecs)
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
        A <- Matrix::sparseMatrix(
            i = i_within, j = gene_map, x = as.double(df$value),
            dims = c(chunk_n, P_hvg), repr = "C"
        )
        # Apply @post_ops R-side (mutates A@x). scalef_vecs are pre-sliced
        # to this substore; chunk-local slice happens inside the executor.
        .pe_apply_post_ops_mat(A, post_ops, info$scalef_vecs, sub_cs, sub_ce)
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


# Helper: per-HVG-gene mean of normalized data. Collects per-substore
# triplets (HVG-column-filtered), applies @post_ops R-side, aggregates
# per-gene sums via data.table by-group, sums across substores.
.stream_norm_hvg_means <- function(pe, hvg_idx, sub_infos = NULL) {
    n_cells  <- as.integer(pe@n_cells)
    P_hvg    <- length(hvg_idx)

    g_sum <- numeric(P_hvg)
    col_id <- value <- s <- NULL  # NSE

    post_ops <- pe@post_ops
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
        df <- storeRead(info$sub, output = "query") |>
            dplyr::filter(col_id %in% !!info$hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) next
        # Apply @post_ops R-side; mutates df$value.
        df <- .pe_apply_post_ops_df(df, post_ops)
        agg <- df[, .(s = sum(value, na.rm = TRUE)), by = col_id]
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

    # ---- Materialize: HVG-only, @post_ops baked, union collapsed ----------
    # Writes a fresh temp pestore whose value column is the normalized
    # (JIT-applied) HVG-subset triplet stream.  Downstream 3-pass PCA reads
    # this temp store directly — no per-band @post_ops apply, no per-band
    # `%in%` HVG filter, no cross-substore composition.  The temp store is
    # cleaned up on function exit.
    cell_ids_final <- pe@cell_ids
    pe_narrow <- pe[hvg_idx, ]
    pe_mat_path <- tempfile("gd_pca_mat_", fileext = ".parquet")
    pe_mat <- storeWrite(
        parquetExprStore(path = pe_mat_path),
        pe_narrow
    )
    on.exit(unlink(pe_mat_path, recursive = TRUE), add = TRUE)

    # Rebuild sub_infos + related state against the materialized store.
    # pe_mat is a fresh single-store parquetExprStore with empty @ops /
    # @post_ops, so parent_ops is empty, hvg_orig is identity, and
    # scalef_vecs is empty.
    pe         <- pe_mat
    n_cells    <- as.integer(pe@n_cells)
    chunk_size <- as.integer(.exprbase_chunk_size(pe))
    hvg_idx    <- seq_len(P_hvg)          # identity on materialized axis
    post_ops   <- list()                   # baked in
    sub_infos  <- lapply(.exprbase_substores(pe), function(se) {
        sub <- se$store
        list(sub         = sub,
             offset      = as.integer(se$cell_offset),
             n_sub       = as.integer(sub@n_cells),
             hvg_orig    = seq_len(P_hvg),
             scalef_vecs = list())
    })

    # Pass 1: means (center) + sds (scale) in one pass over @post_ops.
    stats <- .stream_norm_hvg_stats(pe, hvg_idx, sub_infos = sub_infos)
    means <- if (center) stats$means else numeric(P_hvg)
    sds   <- if (scale)  stats$sds   else rep(1, P_hvg)

    # Per-substore chunk reader (cell-major within substore).
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
        chunk_n  <- sub_ce - sub_cs + 1L
        gene_map <- match(df$col_id, info$hvg_orig)
        cell_map <- match(df$row_id, orig_rows)
        A <- Matrix::sparseMatrix(
            i = cell_map, j = gene_map, x = as.double(df$value),
            dims = c(chunk_n, P_hvg), repr = "C"
        )
        .pe_apply_post_ops_mat(A, post_ops, info$scalef_vecs,
            sub_cs, sub_ce)
    }

    # Pass 2: G_raw = Σ chunk^T chunk across all substores × chunks.
    # Associative accumulation — substore boundaries don't matter.
    # Kept serial: an mclapply-over-bands experiment (2026-07-24) added
    # 3-5 s vs serial due to per-band arrow filter + sparseMatrix rebuild
    # overhead × 8 workers, plus arrow C++ thread-pool oversubscription
    # (n_workers × cpu_count threads competing for cpu_count cores).
    # Pass 2 band-parallel: since input is now the small materialized
    # store, per-band arrow reads are cheap (no `%in%` HVG filter, tiny
    # row-group stats).  Split n_sub into n_workers bands, each returns
    # a partial G_raw = ΣchunkᵀAchunk; reduce by summing.
    n_workers <- .par_workers()
    G_raw <- matrix(0.0, nrow = P_hvg, ncol = P_hvg)
    info1 <- sub_infos[[1L]]
    n_sub1 <- info1$n_sub
    band_size <- max(1L, as.integer(ceiling(n_sub1 / n_workers)))
    band_starts <- seq.int(1L, n_sub1, by = band_size)
    bands <- lapply(band_starts, function(bs)
        c(bs, min(bs + band_size - 1L, n_sub1)))

    gram_band <- function(rng) {
        G_local <- matrix(0.0, nrow = P_hvg, ncol = P_hvg)
        cs <- rng[1L]
        ce_stop <- rng[2L]
        while (cs <= ce_stop) {
            ce <- min(cs + chunk_size - 1L, ce_stop)
            A  <- .read_chunk_sub(info1, cs, ce)
            if (!is.null(A)) {
                G_local <- G_local + as.matrix(Matrix::crossprod(A))
            }
            cs <- ce + 1L
        }
        G_local
    }
    partials <- if (.Platform$OS.type == "unix" && n_workers > 1L) {
        parallel::mclapply(bands, gram_band,
            mc.cores = n_workers, mc.preschedule = TRUE)
    } else {
        lapply(bands, gram_band)
    }
    G_raw <- Reduce("+", partials, init = G_raw)
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
    # per-chunk work is identical for scale on/off. coords is sized to
    # the union cell axis; each band writes its own disjoint row range,
    # main thread copies partials into the global coords matrix.
    V_use <- V / sds
    coords <- matrix(0.0, nrow = n_cells, ncol = ncp_used)
    correction <- if (center) {
        as.numeric(means %*% V_use)
    } else {
        numeric(ncp_used)
    }
    # Pass 3 band-parallel: same shape as Pass 2.  Each worker returns
    # its band's coord slice; main thread copies into the global matrix.
    coords_band <- function(rng) {
        band_n <- rng[2L] - rng[1L] + 1L
        band_coords <- matrix(0.0, nrow = band_n, ncol = ncp_used)
        cs <- rng[1L]
        ce_stop <- rng[2L]
        while (cs <= ce_stop) {
            ce <- min(cs + chunk_size - 1L, ce_stop)
            A  <- .read_chunk_sub(info1, cs, ce)
            chunk_n <- ce - cs + 1L
            in_band <- (cs - rng[1L] + 1L):(cs - rng[1L] + chunk_n)
            if (!is.null(A)) {
                Cc <- as.matrix(A %*% V_use)
                if (center) {
                    Cc <- Cc - matrix(correction, nrow = chunk_n,
                                       ncol = ncp_used, byrow = TRUE)
                }
                band_coords[in_band, ] <- Cc
            } else if (center) {
                band_coords[in_band, ] <- -matrix(correction,
                    nrow = chunk_n, ncol = ncp_used, byrow = TRUE)
            }
            cs <- ce + 1L
        }
        list(rows = (info1$offset + rng[1L]):(info1$offset + rng[2L]),
             band = band_coords)
    }
    coord_partials <- if (.Platform$OS.type == "unix" && n_workers > 1L) {
        parallel::mclapply(bands, coords_band,
            mc.cores = n_workers, mc.preschedule = TRUE)
    } else {
        lapply(bands, coords_band)
    }
    for (p in coord_partials) coords[p$rows, ] <- p$band

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


# Per-HVG-gene mean + sd of normalized data.  Iterates substores,
# collects HVG-column-filtered triplets per substore, applies @post_ops
# R-side, aggregates via data.table by-group, sums across substores.
# Variance via SS identity: σ² = (SS − n·μ²)/(n−1).  `hvg_idx` NULL uses
# every column.  `sub_infos` may be reused from a caller that already
# built them (must include $sub and $hvg_orig).
#
# NOTE: kept serial intentionally. Pass 1's data volume is small (HVG-
# column filter reduces to ~11% of full triplets), so mclapply's fork/
# join overhead + sparseMatrix sort per band cost more than they save.
# Verified by bench: parallel Pass 1 pushed PCA from 20 s → 25 s.
.stream_norm_hvg_stats <- function(pe, hvg_idx, sub_infos = NULL, ...) {
    n_cells  <- as.integer(pe@n_cells)

    if (is.null(hvg_idx)) hvg_idx <- seq_along(pe@feat_ids)
    P_hvg <- length(hvg_idx)

    post_ops <- pe@post_ops
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

    # NSE bindings
    col_id <- value <- s <- ss <- NULL

    g_sum  <- numeric(P_hvg)
    g_ssum <- numeric(P_hvg)

    for (info in sub_infos) {
        df <- storeRead(info$sub, output = "query") |>
            dplyr::filter(col_id %in% !!info$hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) next
        df <- .pe_apply_post_ops_df(df, post_ops)
        agg <- df[, .(
            s  = sum(value, na.rm = TRUE),
            ss = sum(value * value, na.rm = TRUE)
        ), by = col_id]
        if (nrow(agg) > 0L) {
            g_idx <- match(agg$col_id, info$hvg_orig)
            keep  <- !is.na(g_idx)
            g_sum[g_idx[keep]]  <- g_sum[g_idx[keep]] + as.numeric(agg$s[keep])
            g_ssum[g_idx[keep]] <- g_ssum[g_idx[keep]] + as.numeric(agg$ss[keep])
        }
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
