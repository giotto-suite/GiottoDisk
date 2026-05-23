# Tests for streaming filter dispatch:
#   processData(parquetExprStore, filterParam) -> list(feats_keep, cells_keep)
# Result must match Giotto's in-memory two-stage filter exactly.

.tiny_mat <- function(n_genes = 20L, n_cells = 80L,
                       density = 0.5, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}


test_that("processData(parquetExprStore, filterParam) matches in-memory masks", {
    skip_if_not_installed("Giotto")

    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    fp <- Giotto::filterParam(
        expression_threshold   = 1,
        feat_det_in_min_cells  = 5,
        min_det_feats_per_cell = 3)

    masks_pq  <- GiottoClass::processData(pe,  fp)
    masks_mem <- GiottoClass::processData(mat, fp)

    expect_setequal(masks_pq$feats_keep, masks_mem$feats_keep)
    expect_setequal(masks_pq$cells_keep, masks_mem$cells_keep)
})


test_that("two-stage logic — streaming recounts after gene mask is applied", {
    # Build a matrix where this matters: cells whose nfeats *over all genes*
    # would pass min_det_feats_per_cell but fail once low-detected genes
    # are removed by feat_det_in_min_cells.
    skip_if_not_installed("Giotto")
    set.seed(99)
    mat <- Matrix::rsparsematrix(40, 100, density = 0.3,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(mat) <- paste0("g", 1:40)
    colnames(mat) <- paste0("c", 1:100)

    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    fp <- Giotto::filterParam(
        expression_threshold   = 1,
        feat_det_in_min_cells  = 30,
        min_det_feats_per_cell = 8)

    masks_pq  <- GiottoClass::processData(pe,  fp)
    masks_mem <- GiottoClass::processData(mat, fp)

    expect_setequal(masks_pq$feats_keep, masks_mem$feats_keep)
    expect_setequal(masks_pq$cells_keep, masks_mem$cells_keep)
})


test_that("empty-result branches return empty character vectors", {
    skip_if_not_installed("Giotto")
    mat <- .tiny_mat(seed = 7)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # Impossibly tight threshold — nothing should be kept
    fp <- Giotto::filterParam(
        expression_threshold   = 1000,    # no value reaches this
        feat_det_in_min_cells  = 1,
        min_det_feats_per_cell = 1)
    masks <- GiottoClass::processData(pe, fp)
    expect_length(masks$feats_keep, 0L)
    expect_length(masks$cells_keep, 0L)
})
