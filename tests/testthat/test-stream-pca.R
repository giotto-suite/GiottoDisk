# Tests for streaming Halko PCA dispatch:
#   reduceData(parquetExprStore, randomPcaParam)
#       -> list(u, d, v, sdev, eigenvalues)

.tiny_mat <- function(n_genes = 80L, n_cells = 400L,
                       density = 0.4, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", sprintf("%03d", seq_len(n_genes)))
    colnames(m) <- paste0("c", sprintf("%04d", seq_len(n_cells)))
    m
}


.setup_normalized_pe <- function(mat, hvg_idx) {
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # Push a norm_libsize_log op directly (matches in-memory math the
    # tests compare against). Slice constructor lives on the internal
    # namespace; the full op record fuses libsize + log = TRUE.
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    slice <- GiottoDisk:::.pe_norm_libsize_scalef_slice(
        pe, scalef = 1e4 / libsz)
    pe@post_ops <- list(list(
        type   = "norm_libsize_log",
        scalef = slice,
        log    = TRUE,
        base   = 2
    ))
    pe
}


test_that("randomPcaParam errors on a store with no normalize recipe", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    # Raw store, empty @ops -- reduceData requires a norm op.
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        mat)
    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("random", ncp = 5,
                feats_to_use = rownames(mat)[1:20], scale = FALSE)),
        "no normalization recipe"
    )
})


test_that("randomPcaParam errors when feats_to_use is missing", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 2)
    pe  <- .setup_normalized_pe(mat)
    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("random", ncp = 5, scale = FALSE,
                set_seed = TRUE, seed_number = 42L)),
        "feats_to_use is required"
    )
})


test_that("streaming Halko top-k singular values match irlba (Pearson r > 0.99)", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("irlba")

    set.seed(99)
    mat <- .tiny_mat(n_genes = 100, n_cells = 500, density = 0.4, seed = 99)
    pe  <- .setup_normalized_pe(mat)

    # Pick the top 50 most-variable genes as "HVG"
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    sf <- 1e4 / libsz
    mat_norm <- log1p(t(t(mat) * sf)) / log(2)
    vars <- as.numeric(apply(as.matrix(mat_norm), 1, var))
    hvg <- rownames(mat)[order(vars, decreasing = TRUE)][1:50]
    mat_hvg <- mat_norm[hvg, , drop = FALSE]
    hvg_means <- as.numeric(Matrix::rowMeans(mat_hvg))

    NCP <- 10
    ir <- irlba::irlba(t(mat_hvg), nv = NCP, nu = NCP, center = hvg_means)

    pca_pq <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = NCP,
                          feats_to_use = hvg,
                          center = TRUE, scale = FALSE,
                          set_seed = TRUE, seed_number = 42L,
                          n_oversamples = 10L, n_power_iter = 2L))

    expect_equal(length(pca_pq$d), NCP)
    expect_equal(nrow(pca_pq$u), ncol(mat))
    expect_equal(nrow(pca_pq$v), length(hvg))

    pearson <- cor(ir$d, pca_pq$d)
    expect_gt(pearson, 0.99)

    # Top-3 PCs should be very close. On a synthetic 100×500 matrix
    # with q=2 power iterations, ~2 % is the expected Halko range.
    rel_err_top3 <- max(abs(ir$d[1:3] - pca_pq$d[1:3]) / ir$d[1:3])
    expect_lt(rel_err_top3, 0.05)
})


test_that("irlba and exact pcaParam variants error on parquet backend", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 7)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:30]

    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("irlba", ncp = 5, feats_to_use = hvg)),
        "not supported for streaming"
    )
    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("exact", ncp = 5, feats_to_use = hvg)),
        "not supported for streaming"
    )
})


# autoPcaParam: substrate routes to gram-eigen or random depending on
# Gram-fits-in-budget; dry_run = TRUE returns the resolved param.

test_that("reduceData(parquetExprStore, auto + dry_run) resolves gram-eigen when P is small", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:20]

    # Default budget (4 GB); 20^2 * 8 = 3.2 KB -- easily fits, picks gram
    resolved <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("auto", ncp = 5, feats_to_use = hvg,
            dry_run = TRUE))
    expect_s4_class(resolved, "gramEigenPcaParam")
    expect_equal(resolved$ncp, 5L)
    expect_equal(resolved$feats_to_use, hvg)
})

test_that("reduceData(parquetExprStore, auto + dry_run) resolves random when budget too tight", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:20]

    # Force a tiny budget so 20^2 * 8 = 3.2 KB exceeds it -> Halko
    withr::with_options(list(giottodisk.pca_auto_budget_gb = 1e-9), {
        resolved <- GiottoClass::reduceData(pe,
            Giotto::pcaParam("auto", ncp = 5, feats_to_use = hvg,
                dry_run = TRUE))
        expect_s4_class(resolved, "randomPcaParam")
    })
})


# gramEigenPcaParam: constructor lives in GiottoDisk, class inherits
# from Giotto::pcaParam.

test_that("gramEigenPcaParam() constructor populates all knobs", {
    p <- gramEigenPcaParam(ncp = 30, feats_to_use = c("g1", "g2"),
        center = TRUE, scale = FALSE, fallback_relerr = 0.005,
        n_oversamples = 15L, n_power_iter = 3L, seed_number = 7L)
    expect_s4_class(p, "gramEigenPcaParam")
    expect_true(is(p, "pcaParam"))
    expect_equal(p$method, "gram-eigen")
    expect_equal(p$ncp, 30L)
    expect_equal(p$feats_to_use, c("g1", "g2"))
    expect_equal(p$fallback_relerr, 0.005)
    expect_equal(p$n_oversamples, 15L)
    expect_equal(p$seed_number, 7L)
})


# Streaming gram-eigen: parity against dense svd() on the same
# normalized + centered data.

test_that("streaming gram-eigen singular values match svd() to 1e-8", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 11)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:30]

    # Reference: apply the same normalization the streaming path sees,
    # subset HVG rows, center columns, and run dense svd().
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    sf <- 1e4 / libsz
    mat_norm <- log1p(t(t(mat) * sf)) / log(2)
    A <- as.matrix(t(mat_norm[hvg, , drop = FALSE]))   # cells x HVG
    A_c <- scale(A, center = TRUE, scale = FALSE)
    ref <- svd(A_c, nu = 5, nv = 5)

    stream_res <- GiottoClass::reduceData(pe,
        gramEigenPcaParam(ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = FALSE))

    # Gram-eigen on well-conditioned data matches svd() to machine
    # precision (well under the CONVENTIONS 1% bar).
    rel <- max(abs(stream_res$d - ref$d[seq_len(5)]) / ref$d[seq_len(5)])
    expect_lt(rel, 1e-8)
    expect_equal(nrow(stream_res$v), 30L)     # loadings restricted to HVG
    expect_equal(rownames(stream_res$v), hvg)

    # Per-PC score/loading correlation (allow sign flip)
    ref_coords <- ref$u %*% diag(ref$d[seq_len(5)])
    for (k in seq_len(5)) {
        rho_v <- abs(cor(stream_res$v[, k], ref$v[, k]))
        rho_u <- abs(cor(stream_res$u[, k], ref_coords[, k]))
        expect_gt(rho_v, 0.999)
        expect_gt(rho_u, 0.999)
    }
})

test_that("streaming gram-eigen top-k singular values match irlba", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("irlba")
    set.seed(99)
    mat <- .tiny_mat(n_genes = 100, n_cells = 500, density = 0.4, seed = 99)
    pe  <- .setup_normalized_pe(mat)
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    sf <- 1e4 / libsz
    mat_norm <- log1p(t(t(mat) * sf)) / log(2)
    vars <- as.numeric(apply(as.matrix(mat_norm), 1, var))
    hvg <- rownames(mat)[order(vars, decreasing = TRUE)][1:50]
    mat_hvg <- mat_norm[hvg, , drop = FALSE]
    hvg_means <- as.numeric(Matrix::rowMeans(mat_hvg))

    NCP <- 10
    ir <- irlba::irlba(t(mat_hvg), nv = NCP, nu = NCP, center = hvg_means)
    pca_pq <- GiottoClass::reduceData(pe,
        gramEigenPcaParam(ncp = NCP, feats_to_use = hvg,
            center = TRUE, scale = FALSE))

    expect_equal(length(pca_pq$d), NCP)
    # Gram-eigen matches irlba to machine precision on well-conditioned
    # top-k (well below the 1% CONVENTIONS bar).
    rel <- max(abs(ir$d - pca_pq$d) / ir$d)
    expect_lt(rel, 1e-8)
})


# scale = TRUE parity vs svd(scale(A, center = TRUE, scale = TRUE)).

test_that("streaming gram-eigen with scale=TRUE matches svd(scale(A))", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 21)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:30]

    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    sf <- 1e4 / libsz
    mat_norm <- log1p(t(t(mat) * sf)) / log(2)
    A <- as.matrix(t(mat_norm[hvg, , drop = FALSE]))     # cells x HVG
    A_cs <- scale(A, center = TRUE, scale = TRUE)         # standardized
    ref <- svd(A_cs, nu = 5, nv = 5)

    stream_res <- GiottoClass::reduceData(pe,
        gramEigenPcaParam(ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = TRUE))

    rel <- max(abs(stream_res$d - ref$d[seq_len(5)]) / ref$d[seq_len(5)])
    expect_lt(rel, 1e-8)
    ref_coords <- ref$u %*% diag(ref$d[seq_len(5)])
    for (k in seq_len(5)) {
        rho_v <- abs(cor(stream_res$v[, k], ref$v[, k]))
        rho_u <- abs(cor(stream_res$u[, k], ref_coords[, k]))
        expect_gt(rho_v, 0.999)
        expect_gt(rho_u, 0.999)
    }
})

# Gram-eigen on unionParquetExprStore: substore-iterated Pass 2 + Pass 3.
# Correctness bar mirrors the single-store case (machine-precision d;
# per-PC cor > 0.999 vs dense svd on concatenated normalized data).

test_that("streaming gram-eigen on union matches svd() on concat data", {
    skip_if_not_installed("Giotto")
    set.seed(42)
    m1 <- Matrix::rsparsematrix(50, 30, density = 0.4,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m1) <- paste0("g", sprintf("%03d", seq_len(50)))
    colnames(m1) <- paste0("a_c", seq_len(30))
    m2 <- Matrix::rsparsematrix(50, 20, density = 0.4,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m2) <- rownames(m1)
    colnames(m2) <- paste0("b_c", seq_len(20))

    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m1)
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m2)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library", scalefactor = 1e4)) |>
        GiottoClass::processData(Giotto::normParam("log", base = 2, offset = 1))

    hvg <- rownames(m1)[seq_len(20)]
    res <- GiottoClass::reduceData(u,
        gramEigenPcaParam(ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = FALSE))

    # Dense reference: concat + normalize with the same libsize+log2 recipe.
    m_all <- cbind(m1, m2)
    libsz <- as.numeric(Matrix::colSums(m_all))
    libsz[libsz == 0] <- 1
    sf <- 1e4 / libsz
    mat_norm <- log1p(t(t(m_all) * sf)) / log(2)
    A <- as.matrix(t(mat_norm[hvg, , drop = FALSE]))
    A_c <- scale(A, center = TRUE, scale = FALSE)
    ref <- svd(A_c, nu = 5, nv = 5)

    rel <- max(abs(res$d - ref$d[seq_len(5)]) / ref$d[seq_len(5)])
    expect_lt(rel, 1e-8)
    expect_equal(nrow(res$v), 20L)
    expect_equal(nrow(res$u), 50L)  # union cell axis (30 + 20)

    ref_coords <- ref$u %*% diag(ref$d[seq_len(5)])
    for (k in seq_len(5)) {
        rho_v <- abs(cor(res$v[, k], ref$v[, k]))
        rho_u <- abs(cor(res$u[, k], ref_coords[, k]))
        expect_gt(rho_v, 0.999)
        expect_gt(rho_u, 0.999)
    }
})


test_that("streaming Halko errors on scale=TRUE (densification guard)", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:20]
    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("random", ncp = 5, feats_to_use = hvg,
                center = TRUE, scale = TRUE)),
        "scale = TRUE.*not supported"
    )
})
