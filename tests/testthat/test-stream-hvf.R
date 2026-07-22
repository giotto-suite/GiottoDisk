# Tests for streaming HVF dispatch:
#   analyzeData(parquetExprStore, covLoessParam) -> per-feature stats
#     (cov_diff = residual COV above LOESS fit). Selection is a separate
#     downstream step under processData.

.tiny_mat <- function(n_genes = 50L, n_cells = 200L,
                       density = 0.4, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", sprintf("%03d", seq_len(n_genes)))
    colnames(m) <- paste0("c", sprintf("%04d", seq_len(n_cells)))
    m
}


test_that("analyzeData(parquetExprStore, covLoessParam) requires JIT recipe", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(
        GiottoClass::analyzeData(pe, Giotto::analyzeParam("cov_loess")),
        "no normalization recipe"
    )
})


test_that("streaming covLoessParam matches in-memory selection on a tiny matrix", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")

    mat <- .tiny_mat(seed = 11)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # In-memory baseline
    locs <- data.frame(cell_ID = colnames(mat),
                        sdimx = runif(ncol(mat)),
                        sdimy = runif(ncol(mat)))
    g_mem <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                         verbose = FALSE)
    g_mem <- suppressWarnings(Giotto::normalizeGiotto(g_mem,
        scale_feats = FALSE, scale_cells = FALSE, verbose = FALSE))
    g_mem <- suppressWarnings(Giotto::calculateHVF(g_mem,
        method = "cov_loess",
        show_plot = FALSE, return_plot = FALSE, save_plot = FALSE,
        verbose = FALSE))
    fm_mem <- GiottoClass::getFeatureMetadata(g_mem, output = "data.table")
    sel_mem <- fm_mem$feat_ID[!is.na(fm_mem$hvf) & fm_mem$hvf == "yes"]

    # Parquet path
    g_pq <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                        verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g_pq <- GiottoClass::setExpression(g_pq, x = eo, name = "raw",
                                         verbose = FALSE)
    g_pq <- suppressWarnings(Giotto::normalizeGiotto(g_pq,
        scale_feats = FALSE, scale_cells = FALSE, verbose = FALSE))
    g_pq <- suppressWarnings(Giotto::calculateHVF(g_pq,
        method = "cov_loess",
        show_plot = FALSE, return_plot = FALSE, save_plot = FALSE,
        verbose = FALSE))
    fm_pq <- GiottoClass::getFeatureMetadata(g_pq, output = "data.table")
    sel_pq <- fm_pq$feat_ID[!is.na(fm_pq$hvf) & fm_pq$hvf == "yes"]

    # Jaccard >= 0.9 (we expect ~1.0 in practice — the math is identical)
    inter <- length(intersect(sel_mem, sel_pq))
    union <- length(union(sel_mem, sel_pq))
    jaccard <- if (union == 0L) 1 else inter / union
    expect_gte(jaccard, 0.9)
})


test_that("varParam errors clearly on parquet backend (cov_groups is supported)", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 19)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe@post_ops <- list(list(
        type   = "norm_libsize_log",
        scalef = GiottoDisk:::.pe_norm_libsize_scalef_slice(
            pe, scalef = rep(1, ncol(mat))),
        log    = FALSE,
        base   = 2
    ))

    # cov_groups is supported by the streaming backend — returns stats
    # data.table without erroring.
    res <- GiottoClass::analyzeData(pe, Giotto::analyzeParam("cov_groups"))
    expect_s3_class(res, "data.table")
    expect_true("cov_group_zscore" %in% colnames(res))

    # var (per-feature variance on a scaled matrix) requires densifying
    # the streaming reads and is intentionally unsupported.
    expect_error(
        GiottoClass::analyzeData(pe, Giotto::analyzeParam("var")),
        "not supported for streaming"
    )
})


test_that("covLoessParam output schema matches Giotto convention", {
    mat <- .tiny_mat(seed = 23)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe@post_ops <- list(list(
        type   = "norm_libsize_log",
        scalef = GiottoDisk:::.pe_norm_libsize_scalef_slice(
            pe, scalef = rep(1, ncol(mat))),
        log    = FALSE,
        base   = 2
    ))

    dt <- GiottoClass::analyzeData(pe, Giotto::analyzeParam("cov_loess"))

    expect_s3_class(dt, "data.table")
    # intermediate prediction column is dropped (matches Giotto)
    expect_false("pred_cov" %in% colnames(dt))
    expect_false("pred_cov_feats" %in% colnames(dt))
    expect_true("cov_diff" %in% colnames(dt))
    # sorted descending by cov_diff
    expect_equal(dt$cov_diff, sort(dt$cov_diff, decreasing = TRUE))
    # zero-detection features filtered out before fit
    expect_true(all(dt$nr_cells > 0))
})


test_that("covGroupsParam output schema matches Giotto convention", {
    mat <- .tiny_mat(seed = 29)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe@post_ops <- list(list(
        type   = "norm_libsize_log",
        scalef = GiottoDisk:::.pe_norm_libsize_scalef_slice(
            pe, scalef = rep(1, ncol(mat))),
        log    = FALSE,
        base   = 2
    ))

    dt <- GiottoClass::analyzeData(pe, Giotto::analyzeParam("cov_groups"))

    expect_s3_class(dt, "data.table")
    # internal binning column is dropped before return
    expect_false("expr_groups" %in% colnames(dt))
    expect_true("cov_group_zscore" %in% colnames(dt))
    # zero-detection features filtered out before binning
    expect_true(all(dt$nr_cells > 0))
})


test_that("streaming + in-memory analyzeData share column schema", {
    mat <- .tiny_mat(seed = 31)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe@post_ops <- list(list(
        type   = "norm_libsize_log",
        scalef = GiottoDisk:::.pe_norm_libsize_scalef_slice(
            pe, scalef = rep(1, ncol(mat))),
        log    = FALSE,
        base   = 2
    ))

    for (method in c("cov_loess", "cov_groups")) {
        p <- Giotto::analyzeParam(method)
        disk_dt <- GiottoClass::analyzeData(pe, p)
        mem_dt  <- GiottoClass::analyzeData(mat, p)
        expect_setequal(colnames(disk_dt), colnames(mem_dt))
    }
})
