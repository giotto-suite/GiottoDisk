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


test_that("randomPcaParam errors without normalize recipe", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("random", ncp = 5,
                              feats_to_use = rownames(mat)[1:20])),
        "no normalization recipe"
    )
})


test_that("randomPcaParam requires feats_to_use", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 2)
    pe  <- .setup_normalized_pe(mat)
    # scale = FALSE so the feats_to_use check fires (scale check is earlier)
    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("random", ncp = 5, scale = FALSE)),
        "feats_to_use is required"
    )
})


test_that("randomPcaParam errors with scale = TRUE", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- .setup_normalized_pe(mat)
    expect_error(
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("random", ncp = 5,
                              feats_to_use = rownames(mat)[1:20],
                              scale = TRUE)),
        "scale = TRUE"
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


# autoPcaParam ####
# "auto" defers method selection to the substrate's reduceData
# (parquetExprStore, autoPcaParam) method. For parquetExprStore the
# choice is currently "random" (Halko) -- the only streaming-safe path.
# When gramEigenPcaParam lands, that method's body grows a branch.
# `dry_run = TRUE` returns the resolved concrete pcaParam without
# running PCA.

test_that("reduceData(parquetExprStore, auto + dry_run) returns randomPcaParam", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- .setup_normalized_pe(mat)
    resolved <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("auto", ncp = 5,
            feats_to_use = rownames(mat)[1:20], dry_run = TRUE))
    expect_s4_class(resolved, "randomPcaParam")
    expect_equal(resolved$ncp, 5L)
    expect_equal(resolved$feats_to_use, rownames(mat)[1:20])
    # dry_run stripped from the concrete param
    expect_null(resolved$dry_run)
})

test_that("reduceData(parquetExprStore, autoPcaParam) matches randomPcaParam byte-for-byte", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 3)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:30]

    auto_res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("auto", ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = FALSE,
            set_seed = TRUE, seed_number = 42L))
    ref_res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = FALSE,
            set_seed = TRUE, seed_number = 42L))

    expect_equal(auto_res$d, ref_res$d, tolerance = 1e-10)
    expect_equal(dim(auto_res$u), dim(ref_res$u))
    expect_equal(dim(auto_res$v), dim(ref_res$v))
    for (k in seq_along(ref_res$d)) {
        rho <- abs(cor(auto_res$u[, k], ref_res$u[, k]))
        expect_gt(rho, 0.999)
    }
})
