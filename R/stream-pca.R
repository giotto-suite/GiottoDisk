#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
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
# Centering is implicit: column means of normalized data are computed by
# .stream_norm_hvg_means() in a single streaming pass over the @ops chain,
# then subtracted analytically inside the forward/backward passes — no
# densification.
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
        if (!.pe_has_norm_op(x@ops)) {
            stop("[reduceData(parquetExprStore, randomPcaParam)] ",
                 "expression backend has no normalization recipe. Run ",
                 "normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE) ",
                 "first.", call. = FALSE)
        }
        if (isTRUE(param$scale)) {
            stop("[reduceData(parquetExprStore, randomPcaParam)] ",
                 "scale = TRUE (per-gene z-score) is not supported for ",
                 "streaming because it densifies the matrix. Pass ",
                 "scale = FALSE.", call. = FALSE)
        }
        feats <- param$feats_to_use
        if (is.null(feats)) {
            stop("[reduceData(parquetExprStore, randomPcaParam)] ",
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


# ---- Streaming Halko core --------------------------------------------------

.stream_random_svd <- function(pe, k, n_oversamples = 10L, n_power_iter = 2L,
                                feats_to_use, center = TRUE,
                                set_seed = TRUE, seed_number = 1234L) {
    if (set_seed) set.seed(seed_number)

    n_cells <- as.integer(pe@n_cells)
    chunk_size <- as.integer(pe@chunk_size %null% 250000L)

    # Map HVG feature IDs to integer col_ids
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

    # ---- Compute per-HVG-gene normalized means (one streaming pass) -------
    means <- if (center) {
        .stream_norm_hvg_means(pe, hvg_idx)
    } else {
        numeric(P_hvg)
    }

    # Translate subset positions to ORIGINAL parquet col_ids that
    # correspond to the requested HVG features. (When pe is not
    # subsetted, hvg_idx values ARE original col_ids.)
    hvg_orig <- .pe_orig_col(hvg_idx, pe)

    # ---- Streaming chunk reader (cell-major) -----------------------------
    # storeRead(pe) returns the arrow query with @ops composed in (v_norm
    # already projected by the norm_libsize_log op). Tighten the query
    # with the cell band + HVG col filter, then collect just this chunk.
    .read_chunk_norm_hvg <- function(cell_start, cell_end) {
        row_id <- col_id <- NULL  # NSE
        orig_rows <- .pe_orig_row(cell_start:cell_end, pe)
        df <- storeRead(pe, output = "query") |>
            dplyr::filter(row_id %in% !!orig_rows,
                           col_id %in% !!hvg_orig) |>
            dplyr::collect() |>
            data.table::as.data.table()
        if (nrow(df) == 0L) return(NULL)
        chunk_n <- cell_end - cell_start + 1L
        # Map original col_id -> HVG position; original row_id -> within-band
        # position. row_id values from the arrow query are still in
        # original-parquet coords (sf join keyed off orig_row_id).
        gene_map <- match(df$col_id, hvg_orig)
        i_within <- match(df$row_id, orig_rows)
        Matrix::sparseMatrix(
            i = i_within, j = gene_map, x = as.double(df$v_norm),
            dims = c(chunk_n, P_hvg), repr = "C"
        )
    }

    # ---- Forward: Y = (A_norm - 1·means^T) · M  --------------------------
    .forward <- function(M) {
        m <- ncol(M)
        correction <- if (center) as.numeric(means %*% M)
                       else numeric(m)
        Y <- matrix(0.0, nrow = n_cells, ncol = m)
        cs <- 1L
        while (cs <= n_cells) {
            ce <- min(cs + chunk_size - 1L, n_cells)
            A  <- .read_chunk_norm_hvg(cs, ce)
            chunk_n <- ce - cs + 1L
            if (!is.null(A)) {
                Yc <- as.matrix(A %*% M)
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

    # ---- Backward: returns Z = A_norm^T · Y  +  Gram G = Y^T Y -----------
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


# Helper: per-HVG-gene mean of normalized data (one streaming pass).
# Used for implicit centering inside .forward / .backward.
.stream_norm_hvg_means <- function(pe, hvg_idx) {
    n_cells  <- as.integer(pe@n_cells)
    P_hvg    <- length(hvg_idx)

    g_sum <- numeric(P_hvg)
    col_id <- v_norm <- s <- NULL  # NSE

    # Translate hvg_idx (subset positions) to ORIGINAL parquet col_ids
    hvg_orig <- .pe_orig_col(hvg_idx, pe)

    # storeRead returns the arrow query with @ops composed in. Restrict to
    # HVG cols, aggregate sum(v_norm) per gene at arrow layer, then map
    # back to HVG positions.
    agg <- storeRead(pe, output = "query") |>
        dplyr::filter(col_id %in% !!hvg_orig) |>
        dplyr::group_by(col_id) |>
        dplyr::summarise(s = sum(v_norm, na.rm = TRUE)) |>
        dplyr::collect() |>
        data.table::as.data.table()

    if (nrow(agg) > 0L) {
        g_idx <- match(agg$col_id, hvg_orig)
        keep <- !is.na(g_idx)
        g_sum[g_idx[keep]] <- as.numeric(agg$s[keep])
    }
    g_sum / n_cells
}
