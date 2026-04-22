test_that(".is_sparse_matrix: dense base matrix is FALSE", {
    m <- matrix(1:9, nrow = 3)
    expect_false(.is_sparse_matrix(m))
})

test_that(".is_sparse_matrix: mostly-zero base matrix is TRUE", {
    m <- matrix(0, nrow = 10, ncol = 10)
    m[1, 1] <- 1  # 1% non-zero
    expect_true(.is_sparse_matrix(m))
})

test_that(".is_sparse_matrix: Matrix sparse class is TRUE", {
    skip_if_not_installed("Matrix")
    m <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = 1:3, dims = c(10, 10))
    expect_true(.is_sparse_matrix(m))
})

test_that(".is_sparse_matrix: non-matrix returns FALSE", {
    expect_false(.is_sparse_matrix(data.frame(x = 1:3)))
})
