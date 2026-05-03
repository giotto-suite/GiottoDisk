# Tests for streaming HVG dispatch:
#   processData(parquetExprStore, covLoessHvgParam) -> data.table with selected

.tiny_mat <- function(n_genes = 50L, n_cells = 200L,
                       density = 0.4, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", sprintf("%03d", seq_len(n_genes)))
    colnames(m) <- paste0("c", sprintf("%04d", seq_len(n_cells)))
    m
}


test_that("processData(parquetExprStore, covLoessHvgParam) requires JIT recipe", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_error(
        GiottoClass::processData(pe, Giotto::hvgParam("cov_loess")),
        "no normalization recipe"
    )
})


test_that("streaming covLoessHvgParam matches in-memory selection on a tiny matrix", {
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


test_that("cov_groups and var_p_resid error clearly on parquet backend", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 19)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    locs <- data.frame(cell_ID = colnames(mat),
                        sdimx = runif(ncol(mat)),
                        sdimy = runif(ncol(mat)))
    g <- Giotto::createGiottoObject(expression = mat, spatial_locs = locs,
                                     verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)
    g  <- suppressWarnings(Giotto::normalizeGiotto(g,
        scale_feats = FALSE, scale_cells = FALSE, verbose = FALSE))

    expect_error(
        suppressWarnings(Giotto::calculateHVF(g, method = "cov_groups",
            show_plot = FALSE, return_plot = FALSE, save_plot = FALSE,
            verbose = FALSE)),
        "not implemented for the streaming"
    )
    expect_error(
        suppressWarnings(Giotto::calculateHVF(g, method = "var_p_resid",
            show_plot = FALSE, return_plot = FALSE, save_plot = FALSE,
            verbose = FALSE)),
        "not supported for streaming"
    )
})
