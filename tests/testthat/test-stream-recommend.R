# Tests for sc_recommend_chunk — chunk-size suggestion for streaming
# expression-matrix conversions.

test_that("sc_recommend_chunk returns a sensible integer", {
    rec <- sc_recommend_chunk(n_cells = 10000, n_genes = 1000, density = 0.1,
                                verbose = FALSE)
    expect_true(is.integer(rec))
    expect_true(rec >= 10000)   # floor
})
