# The parquet write helpers must not persist data.table's in-memory state:
# arrow restores R attributes on collect, and a stale `sorted` / `index` makes
# data.table answer keyed and indexed subsets wrongly.

.keyed_dt <- function() {
    dt <- data.table::data.table(row_id = c(2L, 1L, 3L), col_id = c(5L, 4L, 6L),
                                 value = as.double(1:3))
    data.table::setkey(dt, row_id, col_id)
    invisible(dt[col_id == 4L])                 # builds an auto-index on col_id
    dt
}

test_that(".write_parquet_file writes no data.table sort or index state", {
    dt <- .keyed_dt()
    expect_false(is.null(attr(dt, "sorted")))
    f <- tempfile(fileext = ".parquet")
    on.exit(unlink(f), add = TRUE)
    .write_parquet_file(dt, f)
    back <- dplyr::collect(arrow::open_dataset(f))
    expect_null(attr(back, "sorted"))
    expect_null(attr(back, "index"))
    expect_equal(as.data.frame(back)[order(back$row_id), ], as.data.frame(dt), ignore_attr = TRUE)
    expect_false(is.null(attr(dt, "sorted")))   # the caller's table is untouched
})

test_that(".write_dataset writes no data.table sort or index state", {
    dt <- .keyed_dt()
    d <- tempfile("ds_")
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    .write_dataset(dt, d)
    back <- dplyr::collect(arrow::open_dataset(d))
    expect_null(attr(back, "sorted"))
    expect_null(attr(back, "index"))
    expect_equal(sort(back$value), sort(dt$value))
})

test_that("a store read back filters correctly after a keyed write", {
    dt <- .keyed_dt()
    f <- tempfile(fileext = ".parquet")
    on.exit(unlink(f), add = TRUE)
    .write_parquet_file(dt, f)
    back <- data.table::setDT(dplyr::collect(arrow::open_dataset(f)))
    expect_equal(back[col_id == 4L, value], 2)
})
