# Starting grid of quadtree tile plans ####
#
# The seed grid has four tiles per unit of aspect ratio; a source whose
# extent has zero width or height (one point, or every point on one line)
# made that Inf or NaN and the write failed.

options("tilework.warn_sequential" = FALSE)

.seed_qs <- function(df) {
    f <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, f)
    as(fileStore(path = f, read_fun = function(x, ...) arrow::open_dataset(x)),
        "queryableStore")
}

test_that(".seed_tile_count is four per unit of aspect ratio", {
    expect_equal(.seed_tile_count(c(100, 50)), 8)
    expect_equal(.seed_tile_count(c(10, 10)), 4)
})

test_that(".seed_tile_count starts from one tile for a zero side", {
    expect_identical(.seed_tile_count(c(0, 10)), 1L)
    expect_identical(.seed_tile_count(c(10, 0)), 1L)
    expect_identical(.seed_tile_count(c(0, 0)), 1L)
})

describe("parquetGeomTileStore <- queryableStore with a degenerate extent", {

    test_that("points on one vertical line write one tile", {
        qs <- .seed_qs(data.frame(x = 5, y = seq(0, 40, by = 10),
            feat = letters[1:5]))
        pgts <- parquetGeomTileStore() |> storeWrite(qs, type = "points",
            id_col = "feat", sdimx = "x", sdimy = "y", threshold = 10L)
        expect_equal(nrow(pgts), 5)
        expect_equal(nrow(pgts@tiles@bounds), 1L)
    })

    test_that("a single point writes one tile", {
        qs <- .seed_qs(data.frame(x = 5, y = 3, feat = "a"))
        pgts <- parquetGeomTileStore() |> storeWrite(qs, type = "points",
            id_col = "feat", sdimx = "x", sdimy = "y", threshold = 10L)
        expect_equal(nrow(pgts), 1)
    })
})
