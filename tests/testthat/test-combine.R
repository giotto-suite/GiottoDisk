test_that("rbind2: two parquetStores produce unionParquetStore", {
    ps1 <- parquetStore() |> storeWrite(mtcars[1:10, ])
    ps2 <- parquetStore() |> storeWrite(mtcars[11:20, ])
    u <- rbind2(ps1, ps2)
    expect_s4_class(u, "unionParquetStore")
    expect_equal(nrow(u), 20)
})

test_that("rbind2: tibble output combines all rows", {
    ps1 <- parquetStore() |> storeWrite(mtcars[1:10, ])
    ps2 <- parquetStore() |> storeWrite(mtcars[11:20, ])
    tbl <- rbind2(ps1, ps2) |> storeRead(output = "tibble")
    expect_equal(nrow(tbl), 20)
    expect_equal(sort(colnames(tbl)), sort(colnames(mtcars)))
})

test_that("rbind2: errors on store with pending ops", {
    ps1 <- parquetStore() |> storeWrite(mtcars[1:10, ])
    ps2 <- parquetStore() |> storeWrite(mtcars[11:20, ])
    ps1_sub <- subset(ps1, cyl == 4)
    expect_error(rbind2(ps1_sub, ps2), regexp = "pending ops")
})

test_that("rbind2: self-rbind doubles rows", {
    ps <- parquetStore() |> storeWrite(mtcars)
    u <- rbind2(ps, ps)
    expect_s4_class(u, "unionParquetStore")
    expect_equal(nrow(u), 64)
})

test_that("rbind2: parquetStore + unionParquetStore extends union", {
    ps1 <- parquetStore() |> storeWrite(mtcars[1:10, ])
    ps2 <- parquetStore() |> storeWrite(mtcars[11:20, ])
    ps3 <- parquetStore() |> storeWrite(mtcars[21:32, ])
    u <- rbind2(ps1, ps2)
    u2 <- rbind2(u, ps3)
    expect_s4_class(u2, "unionParquetStore")
    expect_equal(nrow(u2), 32)
})

test_that("rbind2: cannot combine parquetGeomStore with parquetStore", {
    ps <- parquetStore() |> storeWrite(mtcars)
    pts <- terra::vect(
        data.frame(x = 1:3, y = 1:3),
        geom = c("x", "y"), crs = ""
    )
    pgs <- parquetGeomStore() |> storeWrite(pts)
    expect_error(rbind2(ps, pgs))
})
