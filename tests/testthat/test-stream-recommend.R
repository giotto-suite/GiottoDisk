# Tests for sc_recommend_chunk — chunk-size suggestion for streaming
# expression-matrix conversions.

test_that("sc_recommend_chunk returns a sensible integer", {
    rec <- sc_recommend_chunk(n_cells = 10000, n_genes = 1000, density = 0.1,
                                verbose = FALSE)
    expect_true(is.integer(rec))
    expect_true(rec >= 10000)   # floor
})

test_that("storeWrite finalizer auto-tunes pe@chunk_size", {
    # tiny dense matrix → density 1.0; sc_recommend_chunk should give
    # a number well above the prototype default of 250000 on a 4x6
    # matrix (RAM-bounded recommendation will saturate at n_cells).
    mat <- Matrix::Matrix(matrix(seq_len(24), nrow = 4, ncol = 6),
                          sparse = TRUE)
    pe  <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    expect_true(is.numeric(pe@chunk_size))
    expect_gte(pe@chunk_size, pe@n_cells)
})
