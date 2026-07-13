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

    # Apply normalize recipe to the store (matches in-memory math)
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    pe@params$norm <- list(
        method        = "library_size",
        scalefactor   = 1e4,
        scale_factors = 1e4 / libsz,
        log           = TRUE,
        base          = 2,
        offset        = 1
    )
    pe
}


test_that("randomPcaParam works on a store with no normalize recipe", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    # Raw store, no @params$norm set — chunk reader treats values as-is.
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        mat)
    res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 5,
            feats_to_use = rownames(mat)[1:20], scale = FALSE))
    expect_named(res, c("u", "d", "v", "sdev", "eigenvalues"),
        ignore.order = TRUE)
    expect_equal(length(res$d), 5L)
})


test_that("randomPcaParam works with feats_to_use = NULL (all features)", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 2)
    pe  <- .setup_normalized_pe(mat)
    res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 5, scale = FALSE,
            set_seed = TRUE, seed_number = 42L))
    # Loadings restricted to full feat set
    expect_equal(nrow(res$v), nrow(mat))
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

test_that("streaming Halko with scale=TRUE matches svd(scale(A)) on top-k", {
    skip_if_not_installed("Giotto")
    set.seed(31)
    mat <- .tiny_mat(n_genes = 100, n_cells = 500, density = 0.4, seed = 31)
    pe  <- .setup_normalized_pe(mat)
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    sf <- 1e4 / libsz
    mat_norm <- log1p(t(t(mat) * sf)) / log(2)
    vars <- as.numeric(apply(as.matrix(mat_norm), 1, var))
    hvg <- rownames(mat)[order(vars, decreasing = TRUE)][1:50]
    A <- as.matrix(t(mat_norm[hvg, , drop = FALSE]))
    A_cs <- scale(A, center = TRUE, scale = TRUE)
    ref <- svd(A_cs, nu = 3, nv = 3)

    pca_pq <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = TRUE,
            set_seed = TRUE, seed_number = 42L,
            n_oversamples = 10L, n_power_iter = 3L))

    # Halko is approximate; singular value magnitudes are the correctness
    # signal for scale=TRUE wiring. On random-Poisson data the top-k
    # spectrum is near-degenerate (d1/d3 ~ 1.05) and individual loadings
    # rotate freely among near-tied eigenvalues -- correlate only top-1,
    # where the gap is largest.
    rel_top3 <- max(abs(pca_pq$d[1:3] - ref$d[1:3]) / ref$d[1:3])
    expect_lt(rel_top3, 0.05)
    rho_v1 <- abs(cor(pca_pq$v[, 1], ref$v[, 1]))
    expect_gt(rho_v1, 0.90)
})
