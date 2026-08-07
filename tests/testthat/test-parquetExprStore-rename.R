# Tests for renaming a parquetExprStore's axes, with and without a pending
# subset.
#
# The invariant under test: `@feat_ids` / `@cell_ids` name the CURRENT VIEW,
# while `@gene_idx` / `@cell_idx` hold on-disk positions and `@stats` holds
# file-coordinate marginals. A rename touches only the first pair, so values
# keep reading correctly and the file-side state is untouched.

.rn_mat <- function(n_genes = 6L, n_cells = 5L, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.8,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    m
}

.rn_store <- function(m) {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}


test_that("rename on an unsliced store updates both axes", {
    m  <- .rn_mat()
    pe <- .rn_store(m)

    rownames(pe) <- paste0("F", seq_len(nrow(pe)))
    colnames(pe) <- paste0("C", seq_len(ncol(pe)))

    expect_equal(dimnames(pe), list(paste0("F", 1:6), paste0("C", 1:5)))
    got <- as.matrix(storeRead(pe, output = "dgcmatrix"))
    expect_equal(rownames(got), paste0("F", 1:6))
    expect_equal(colnames(got), paste0("C", 1:5))
    expect_equal(unname(got), unname(as.matrix(m)))
})


test_that("rename on a sliced store renames the view, not the file mapping", {
    m  <- .rn_mat()
    pe <- .rn_store(m)
    s  <- pe[c(2L, 5L), c(1L, 3L)]

    gi <- s@gene_idx
    ci <- s@cell_idx
    st <- s@stats

    rownames(s) <- c("A", "B")
    colnames(s) <- c("X", "Y")

    # names are the view's; on-disk positions and file marginals untouched
    expect_equal(s@feat_ids, c("A", "B"))
    expect_equal(s@cell_ids, c("X", "Y"))
    expect_identical(s@gene_idx, gi)
    expect_identical(s@cell_idx, ci)
    expect_identical(s@stats, st)

    # values still resolve through the unchanged on-disk mapping
    got <- as.matrix(storeRead(s, output = "dgcmatrix"))
    expect_equal(unname(got), unname(as.matrix(m)[c(2L, 5L), c(1L, 3L)]))
    expect_equal(dimnames(got), list(c("A", "B"), c("X", "Y")))
})


test_that("a renamed sliced store subsets by its new names", {
    m  <- .rn_mat()
    s  <- .rn_store(m)[c(2L, 5L), c(1L, 3L)]
    dimnames(s) <- list(c("A", "B"), c("X", "Y"))

    sa <- s["A", ]
    expect_equal(sa@feat_ids, "A")
    expect_equal(sa@gene_idx, 2L)          # still points at on-disk gene 2
    expect_equal(unname(as.matrix(storeRead(sa, output = "dgcmatrix"))),
                 unname(as.matrix(m)[2L, c(1L, 3L), drop = FALSE]))

    sy <- s[, "Y"]
    expect_equal(sy@cell_ids, "Y")
    expect_equal(sy@cell_idx, 3L)          # still points at on-disk cell 3
})


test_that("renames survive storeWrite and reset the on-disk mapping", {
    m  <- .rn_mat()
    s  <- .rn_store(m)[c(2L, 5L), c(1L, 3L)]
    dimnames(s) <- list(c("A", "B"), c("X", "Y"))

    w <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), s)

    expect_equal(w@feat_ids, c("A", "B"))
    expect_equal(w@cell_ids, c("X", "Y"))
    # written fresh, so the view IS the file: no subset state left over
    expect_length(w@gene_idx, 0L)
    expect_length(w@cell_idx, 0L)
    # marginals are recomputed against the new file, not inherited
    expect_length(w@stats$col_nnz, 2L)
    expect_length(w@stats$row_nnz, 2L)
    expect_equal(unname(as.matrix(storeRead(w, output = "dgcmatrix"))),
                 unname(as.matrix(m)[c(2L, 5L), c(1L, 3L)]))
})


test_that("assigning NULL for one axis leaves the other alone", {
    pe <- .rn_store(.rn_mat())
    dimnames(pe) <- list(paste0("F", 1:6), NULL)

    expect_equal(pe@feat_ids, paste0("F", 1:6))
    expect_equal(pe@cell_ids, paste0("c", 1:5))
})


test_that("a wrong-length rename is rejected at assignment", {
    m  <- .rn_mat()
    pe <- .rn_store(m)

    expect_error(rownames(pe) <- letters[1:5], "feature names length")
    expect_error(colnames(pe) <- letters[1:9], "cell names length")

    # and on a sliced store, where the view is what the length must match
    s <- pe[c(2L, 5L), c(1L, 3L)]
    expect_error(rownames(s) <- letters[1:6], "feature names length")
    expect_error(colnames(s) <- letters[1:5], "cell names length")
    expect_silent(rownames(s) <- c("A", "B"))

    # the failed assignments left the store intact
    expect_equal(pe@feat_ids, paste0("g", 1:6))
    expect_equal(pe@cell_ids, paste0("c", 1:5))
})
