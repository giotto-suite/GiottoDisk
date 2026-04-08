make_pts <- function(n = 5) {
    terra::vect(
        data.frame(
            x = seq_len(n),
            y = seq_len(n),
            cell_ID = paste0("cell_", seq_len(n))
        ),
        geom = c("x", "y"),
        crs = ""
    )
}

test_that("parquetGeomStore: unwritten store does not exist", {
    pgs <- parquetGeomStore()
    expect_false(storeExists(pgs))
})

test_that("parquetGeomStore: write SpatVector and existence", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts())
    expect_true(storeExists(pgs))
})

test_that("parquetGeomStore: nrow after write", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    expect_equal(nrow(pgs), 5)
})

test_that("parquetGeomStore: tibble output with fields drops unselected special cols", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    tbl <- storeRead(pgs, fields = "cell_ID", output = "tibble")
    special <- c("x_index", "y_index", "geom", "row_index", "source_id")
    expect_false(any(special %in% colnames(tbl)))
    expect_true("cell_ID" %in% colnames(tbl))
})

test_that("parquetGeomStore: terra output is SpatVector", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    sv <- storeRead(pgs, output = "terra")
    expect_s4_class(sv, "SpatVector")
    expect_equal(nrow(sv), 5)
    expect_true("cell_ID" %in% names(sv))
})

test_that("parquetGeomStore: sf output is sf object", {
    skip_if_not_installed("sf")
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    sfobj <- storeRead(pgs, output = "sf")
    expect_true(inherits(sfobj, "sf"))
    expect_equal(nrow(sfobj), 5)
})

test_that("parquetGeomStore: extent filtering reduces rows", {
    pts <- terra::vect(
        data.frame(x = 1:10, y = 1:10),
        geom = c("x", "y"),
        crs = ""
    )
    pgs <- parquetGeomStore() |> storeWrite(pts)
    e <- terra::ext(c(3, 7, 3, 7))
    tbl <- storeRead(pgs, extent = e, output = "tibble")
    expect_gt(nrow(tbl), 0)
    expect_lt(nrow(tbl), 10)
})

test_that("parquetGeomStore: special cols not in colnames", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts(3))
    cn <- colnames(pgs)
    expect_false(any(c("row_index", "source_id", "x_index", "y_index", "geom") %in% cn))
})
