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

test_that("plan workers are capped by the memory budget, explicit ones are not", {
    skip_if_not_installed("future")
    withr::defer(future::plan(future::sequential))
    future::plan(future::multisession, workers = 4L)
    withr::local_options(list(giottodisk.par_workers = NULL,
                              giottodisk.chunk_ram_frac = NULL))

    # a budget below one worker's share runs serially
    withr::local_options(list(giottodisk.chunk_ram_frac = 1e-6))
    expect_identical(.par_workers(), 1L)
    # never more than the plan offers
    expect_lte(.workers_in_budget(4L), 4L)
    # an explicit count is used as given
    withr::local_options(list(giottodisk.par_workers = 3L))
    expect_identical(.par_workers(), 3L)
})

test_that("the window floors are options", {
    withr::local_options(giottodisk.chunk_min_rows = NULL,
        giottodisk.chunk_min_frac = NULL)
    expect_equal(GiottoDisk:::.chunk_min_rows(), 10000)
    expect_equal(GiottoDisk:::.chunk_min_frac(), 0.05)

    # a budget too small for even one row lands on the floor
    local_mocked_bindings(.sc_free_ram = function() 1e6, .package = "GiottoDisk")
    rec <- function() GiottoDisk:::.recommend_chunk_size(
        n_cells = 1e6, n_genes = 2e4, density = 0.1, k = 0L,
        bytes_per_nz = 64)
    expect_identical(rec(), 10000L)
    withr::local_options(giottodisk.chunk_min_rows = 500)
    expect_identical(rec(), 500L)

    # the budget floor is a share of free RAM
    withr::local_options(giottodisk.chunk_min_frac = 0.5)
    expect_equal(GiottoDisk:::.chunk_budget(1e9, 0.01, n_cells = 1e7), 5e8)
})

test_that("storeChunkProbe measures one window and sizes the budget", {
    mat <- Matrix::rsparsematrix(30, 400, density = 0.3,
        rand.x = function(n) as.double(seq_len(n)))
    rownames(mat) <- paste0("g", seq_len(30))
    colnames(mat) <- paste0("c", seq_len(400))
    pe <- storeWrite(
        parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    local_mocked_bindings(.sc_free_ram = function() 4e9, .package = "GiottoDisk")
    p <- storeChunkProbe(pe, cells = 100, verbose = FALSE)
    expect_s3_class(p, "data.table")
    expect_identical(nrow(p), 1L)
    expect_identical(p$cells, 100L)
    expect_identical(p$stored_values, length(mat[, 1:100]@x))
    expect_true(p$seconds >= 0 && p$peak_mb >= 0)
    expect_equal(p$free_mb, 4000L)
    expect_true(p$suggested_cells >= 1L && p$suggested_cells <= 400L)
    expect_identical(p$n_windows, as.integer(ceiling(400 / p$suggested_cells)))

    # more cells than the store has probes the whole store
    expect_identical(storeChunkProbe(pe, cells = 1e6, verbose = FALSE)$cells,
        400L)
    expect_output(storeChunkProbe(pe, cells = 50), "storeChunkProbe")
})
