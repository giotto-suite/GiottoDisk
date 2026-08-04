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


# A matrix with real PC structure: `nfac` gene modules co-varying with `nfac`
# interleaved cell groups on a Poisson background. `.tiny_mat` is structureless
# noise whose singular values past the first sit within ~1-3 % of each other,
# so its singular VECTORS are not identifiable and no loading comparison is
# well posed on it. Here the top three PCs are separated by wide margins
# (measured gaps ~6.2, ~1.31, ~2.44) and their loadings ARE identifiable;
# PC4 onward is degenerate, hence the top-3 bound in tests that use this.
.structured_mat <- function(n_genes = 60L, n_cells = 300L, seed = 1L,
                             nfac = 3L, strength = 6L) {
    set.seed(seed)
    grp  <- rep(seq_len(nfac), length.out = n_cells)
    base <- matrix(rpois(n_genes * n_cells, 3L), n_genes, n_cells)
    mods <- split(seq_len(n_genes), rep(seq_len(nfac), length.out = n_genes))
    for (f in seq_len(nfac)) {
        rows <- mods[[f]]
        cols <- grp == f
        base[rows, cols] <- base[rows, cols] +
            rpois(length(rows) * sum(cols), strength * f)
    }
    m <- methods::as(methods::as(base, "dgCMatrix"), "generalMatrix")
    dimnames(m) <- list(paste0("g", sprintf("%03d", seq_len(n_genes))),
                        paste0("c", sprintf("%04d", seq_len(n_cells))))
    m
}

# Dense normalized reference matching what the streaming path sees, on an
# HVG-RANKED (deliberately non-ascending) selection so the gene axis exercises
# the unsorted-index path through `[` / .pe_axis_pred.
.dense_hvg_ref <- function(mat, n_hvg) {
    libsz <- as.numeric(Matrix::colSums(mat))
    libsz[libsz == 0] <- 1
    mat_norm <- log1p(t(t(mat) * (1e4 / libsz))) / log(2)
    vars <- as.numeric(apply(as.matrix(mat_norm), 1, var))
    hvg  <- rownames(mat)[order(vars, decreasing = TRUE)][seq_len(n_hvg)]
    list(hvg = hvg, A = as.matrix(mat_norm[hvg, , drop = FALSE]))
}


# Halko imposes neither a normalization requirement nor an HVF requirement --
# PCA of whatever the store holds is the caller's business, and Halko is the
# large-feature-space fallback, so demanding a feature selection would be
# backwards. Both used to error; these pin the freedom.
test_that("randomPcaParam runs on a raw (un-normalized) store", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        mat)   # raw store, empty @ops / @post_ops
    hvg <- rownames(mat)[1:20]

    res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = FALSE, set_seed = TRUE, seed_number = 42L))

    # Reference: dense SVD of the same raw counts, centered.
    A <- as.matrix(t(mat[hvg, , drop = FALSE]))
    ref <- svd(scale(A, center = TRUE, scale = FALSE))$d[seq_len(5)]
    expect_lt(abs(res$d[1] - ref[1]) / ref[1], 0.05)
    expect_identical(rownames(res$v), hvg)
    expect_equal(nrow(res$u), ncol(mat))
})


test_that("randomPcaParam runs with feats_to_use = NULL (all features)", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 2)
    pe  <- .setup_normalized_pe(mat)

    res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 5, scale = FALSE,
            center = TRUE, set_seed = TRUE, seed_number = 42L))

    # Every feature is retained, in store order.
    expect_equal(nrow(res$v), nrow(mat))
    expect_identical(rownames(res$v), rownames(mat))
    expect_equal(nrow(res$u), ncol(mat))
    expect_equal(length(res$d), 5L)
})


# The transient bake is a performance choice only -- both sides of the gate
# must produce the same answer. `giottodisk.pca_bake_max_ratio` is forced to
# each extreme to take the two branches on identical input.
test_that("transient bake and un-baked paths agree", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 5)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:20]

    run <- function(ratio) withr::with_options(
        list(giottodisk.pca_bake_max_ratio = ratio),
        GiottoClass::reduceData(pe,
            Giotto::pcaParam("random", ncp = 5, feats_to_use = hvg,
                center = TRUE, scale = FALSE, set_seed = TRUE,
                seed_number = 42L)))

    baked   <- run(1)   # ratio 20/80 = 0.25 <= 1  -> bakes
    unbaked <- run(0)   #                  > 0     -> does not

    expect_equal(baked$d, unbaked$d, tolerance = 1e-10)
    expect_equal(baked$u, unbaked$u, tolerance = 1e-10)
    expect_equal(baked$v, unbaked$v, tolerance = 1e-10)
    expect_identical(rownames(baked$v), rownames(unbaked$v))
    expect_identical(rownames(baked$u), rownames(unbaked$u))
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


# scale = TRUE used to be rejected on Halko as a "densification guard". That
# reasoning was wrong: per-gene scaling is a column scale of a sparse matrix
# and preserves sparsity. Centering is what densifies, and it was already
# handled analytically. Gram supported scaling all along.
test_that("streaming Halko accepts scale=TRUE", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:20]
    res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 5, feats_to_use = hvg,
            center = TRUE, scale = TRUE))
    expect_equal(length(res$d), 5L)
    expect_equal(nrow(res$u), ncol(mat))
    expect_identical(rownames(res$v), hvg)
    expect_true(all(is.finite(res$d)))
})


# Halko with center = FALSE. Separate from the centered parity test above
# because the two take different branches: with center = FALSE no per-gene
# means are accumulated and the rank-1 centering correction is never applied,
# so a fault in the deferred-centering logic need not show up on the centered
# path. Reference is an uncentered dense SVD.
#
# Uses `.structured_mat` so the top PCs are separated and their loadings are
# identifiable. That is what lets the per-PC `v` comparison below carry weight:
# `d` alone is invariant to a gene permutation, so a fault that reordered the
# gene axis while leaving rownames HVG-ranked would satisfy every singular
# value assertion and still return mislabeled loadings.
test_that("streaming Halko with center=FALSE matches uncentered svd()", {
    skip_if_not_installed("Giotto")

    mat <- .structured_mat(seed = 7L)
    pe  <- .setup_normalized_pe(mat)
    ref <- .dense_hvg_ref(mat, n_hvg = 30L)

    NCP <- 6L
    sv  <- svd(t(ref$A))

    res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = NCP, feats_to_use = ref$hvg,
            center = FALSE, scale = FALSE, set_seed = TRUE,
            seed_number = 42L, n_oversamples = 10L, n_power_iter = 2L))

    expect_lt(max(abs(res$d[1:3] - sv$d[1:3]) / sv$d[1:3]), 0.02)
    expect_identical(rownames(res$v), ref$hvg)
    for (k in seq_len(3)) {
        expect_gt(abs(cor(res$v[, k], sv$v[, k])), 0.99)
    }
    expect_equal(nrow(res$u), ncol(mat))
})


# Halko on a union store. The suite's other union PCA test is gram-eigen,
# which takes a different route (it bakes a transient store first). Halko
# reads the substores directly, so this covers parent @ops / @post_ops
# injection and the per-substore gene narrowing.
test_that("streaming Halko on union matches svd() on concat data", {
    skip_if_not_installed("Giotto")

    m1 <- .structured_mat(n_genes = 60L, n_cells = 180L, seed = 11L)
    m2 <- .structured_mat(n_genes = 60L, n_cells = 120L, seed = 12L)
    rownames(m2) <- rownames(m1)
    colnames(m1) <- paste0("a_c", seq_len(ncol(m1)))
    colnames(m2) <- paste0("b_c", seq_len(ncol(m2)))
    both <- cbind(m1, m2)

    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m1)
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m2)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library", scalefactor = 1e4)) |>
        GiottoClass::processData(Giotto::normParam("log", base = 2, offset = 1))

    ref <- .dense_hvg_ref(both, n_hvg = 30L)
    NCP <- 6L
    sv  <- svd(scale(t(ref$A), center = TRUE, scale = FALSE))

    res <- GiottoClass::reduceData(u,
        Giotto::pcaParam("random", ncp = NCP, feats_to_use = ref$hvg,
            center = TRUE, scale = FALSE, set_seed = TRUE,
            seed_number = 42L, n_oversamples = 10L, n_power_iter = 2L))

    expect_lt(max(abs(res$d[1:3] - sv$d[1:3]) / sv$d[1:3]), 0.02)

    # Cells must land on the union axis in substore order, and gene loadings
    # must match their labels row-for-row.
    #
    # Loadings are checked on the top TWO PCs only: concatenating two
    # differently-seeded blocks gives d = ~105, 70, 31, 29, ... so the gap
    # behind PC3 is only ~1.08 and Halko at q = 2 resolves it to |cor| ~0.988.
    # That is randomized-method accuracy on a near-degenerate direction, not a
    # fault -- gram-eigen on the same input returns 1.0000 for all four.
    expect_equal(nrow(res$u), ncol(both))
    expect_identical(rownames(res$u), colnames(both))
    expect_identical(rownames(res$v), ref$hvg)
    for (k in seq_len(2)) {
        expect_gt(abs(cor(res$v[, k], sv$v[, k])), 0.99)
    }
})


# scale = TRUE on Halko. Per-gene scaling does NOT densify -- only centering
# does, and that is already analytic -- so this is supported the same way gram
# supports it: sigma is absorbed into M (A_std*M == A*(D^-1 M)) and divides Z's
# rows afterward. Sigma cannot be deferred like the centering term because it
# sits inside the product, so a stats pass supplies it first.
test_that("streaming Halko with scale=TRUE matches svd(scale(A))", {
    skip_if_not_installed("Giotto")

    mat <- .structured_mat(seed = 21L)
    pe  <- .setup_normalized_pe(mat)
    ref <- .dense_hvg_ref(mat, n_hvg = 30L)

    sv <- svd(scale(t(ref$A), center = TRUE, scale = TRUE))
    res <- GiottoClass::reduceData(pe,
        Giotto::pcaParam("random", ncp = 6L, feats_to_use = ref$hvg,
            center = TRUE, scale = TRUE, set_seed = TRUE,
            seed_number = 42L, n_oversamples = 10L, n_power_iter = 2L))

    expect_lt(max(abs(res$d[1:3] - sv$d[1:3]) / sv$d[1:3]), 0.05)
    expect_identical(rownames(res$v), ref$hvg)
    expect_equal(nrow(res$u), ncol(mat))
    for (k in seq_len(2)) {
        expect_gt(abs(cor(res$v[, k], sv$v[, k])), 0.99)
    }
})


# gram delegates to Halko when the condition-squared blowup exceeds
# fallback_relerr. That call passed `ncp` and `scale`, neither of which
# .stream_random_svd accepted, so the path errored for every input. No test
# reached it because the default threshold is rarely tripped.
test_that("gram-eigen falls back to Halko without erroring", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat(seed = 9)
    pe  <- .setup_normalized_pe(mat)
    hvg <- rownames(mat)[1:25]

    for (sc in c(FALSE, TRUE)) {
        res <- suppressWarnings(GiottoClass::reduceData(pe,
            gramEigenPcaParam(ncp = 5, feats_to_use = hvg,
                center = TRUE, scale = sc,
                fallback_relerr = 0)))   # 0 -> always delegate
        expect_equal(length(res$d), 5L)
        expect_equal(nrow(res$u), ncol(mat))
        expect_identical(rownames(res$v), hvg)
    }

    # and it warns about the delegation rather than doing it silently
    expect_warning(
        GiottoClass::reduceData(pe, gramEigenPcaParam(ncp = 5,
            feats_to_use = hvg, center = TRUE, scale = FALSE,
            fallback_relerr = 0)),
        "Delegating to streaming random SVD"
    )
})
