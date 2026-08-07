# Tests for chunk sizing: the internal window computation and the
# storeChunkInfo report built on it.

test_that(".recommend_chunk_size returns a sensible integer", {
    rec <- GiottoDisk:::.recommend_chunk_size(
        n_cells = 10000, n_genes = 1000, density = 0.1)
    expect_true(is.integer(rec))
    expect_true(rec >= 10000)   # floor
})

test_that("storeWrite finalizer caches on-disk marginals on @stats", {
    # 4 genes x 6 cells, fully dense
    mat <- Matrix::Matrix(matrix(seq_len(24), nrow = 4, ncol = 6),
                          sparse = TRUE)
    pe  <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    m <- as.matrix(mat)
    expect_equal(pe@stats$col_nnz, unname(rowSums(m != 0)))
    expect_equal(pe@stats$row_nnz, unname(colSums(m != 0)))
    # marginals are keyed by on-disk id, so `[` never invalidates them
    expect_identical(pe[c(1L, 3L), ]@stats, pe@stats)
})


test_that("storeChunkInfo reports both read shapes from cached marginals", {
    mat <- Matrix::rsparsematrix(20, 12, density = 0.5,
        rand.x = function(n) as.double(seq_len(n)))
    rownames(mat) <- paste0("g", seq_len(20))
    colnames(mat) <- paste0("c", seq_len(12))
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    info <- storeChunkInfo(pe, verbose = FALSE)
    expect_s3_class(info, "data.table")
    expect_setequal(unique(info$pass), c("sparse matrix", "triplet frame"))
    # the configured ram_frac is always tabulated and marked
    expect_equal(sum(info$current == "<--"), 2L)
    # windows never exceed the store and shrink as the budget shrinks
    expect_true(all(info$chunk_rows <= pe@n_cells))
    sparse <- info[info$pass == "sparse matrix", ]
    expect_false(is.unsorted(sparse$chunk_rows[order(sparse$ram_frac)]))

    # the triplet pass budgets ~4x the memory per stored value, so its window
    # is never larger than the sparse pass at the same budget
    m <- merge(info[info$pass == "sparse matrix", c("ram_frac", "chunk_rows")],
               info[info$pass == "triplet frame", c("ram_frac", "chunk_rows")],
               by = "ram_frac")
    expect_true(all(m$chunk_rows.y <= m$chunk_rows.x))
})

test_that("storeChunkInfo errors on a store with no shape", {
    expect_error(
        storeChunkInfo(parquetExprStore(path = tempfile())),
        "no cells or no features"
    )
})
