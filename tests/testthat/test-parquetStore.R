test_that("parquetStore: unwritten store does not exist", {
    ps <- parquetStore()
    expect_false(storeExists(ps))
})

test_that("parquetStore: write and existence", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_true(storeExists(ps))
})

test_that("parquetStore: dimensions after write", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_equal(nrow(ps), 32)
    expect_equal(ncol(ps), ncol(mtcars))
    expect_equal(dim(ps), c(32, ncol(mtcars)))
})

test_that("parquetStore: special cols hidden from colnames", {
    ps <- parquetStore() |> storeWrite(mtcars)
    cn <- colnames(ps)
    expect_false(any(c("row_index", "source_id") %in% cn))
    expect_true(all(colnames(mtcars) %in% cn))
})

test_that("parquetStore: storeRead query output is lazy arrow object", {
    ps <- parquetStore() |> storeWrite(mtcars)
    q <- storeRead(ps)
    expect_true(inherits(q, "ArrowObject"))
})

test_that("parquetStore: storeRead tibble roundtrip", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- storeRead(ps, output = "tibble")
    expect_s3_class(tbl, "data.frame")
    expect_equal(nrow(tbl), 32)
    expect_equal(sort(colnames(tbl)), sort(colnames(mtcars)))
})

test_that("parquetStore: subset op is recorded lazily", {
    ps <- parquetStore() |> storeWrite(mtcars)
    ps2 <- subset(ps, cyl == 4)
    expect_length(ps2@ops, 1L)
    # original store unaffected
    expect_length(ps@ops, 0L)
})

test_that("parquetStore: subset filters correctly at read time", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- subset(ps, cyl == 4) |> storeRead(output = "tibble")
    expect_true(all(tbl$cyl == 4))
    expect_equal(nrow(tbl), sum(mtcars$cyl == 4))
})

test_that("parquetStore: subset with local variable is inlined", {
    ps <- parquetStore() |> storeWrite(mtcars)
    target_cyl <- 6
    tbl <- subset(ps, cyl == target_cyl) |> storeRead(output = "tibble")
    expect_true(all(tbl$cyl == 6))
})

test_that("parquetStore: head limits rows", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- head(ps, 5) |> storeRead(output = "tibble")
    expect_lte(nrow(tbl), 5)
})

test_that("parquetStore: [,j] column selection", {
    ps <- parquetStore() |> storeWrite(mtcars)
    ps2 <- ps[, c("mpg", "cyl")]
    tbl <- storeRead(ps2, output = "tibble")
    expect_equal(sort(colnames(tbl)), sort(c("mpg", "cyl")))
})

test_that("parquetStore: nrow returns numeric not integer", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_type(nrow(ps), "double")
})
