# Helpers ####

options("tilework.warn_sequential" = FALSE)


# 5x5 grid of points covering [0, 50] x [0, 50]
make_tile_pts <- function() {
    coords <- data.frame(
        x = rep(seq(5, 45, by = 10), 5),
        y = rep(seq(5, 45, by = 10), each = 5)
    )
    terra::vect(
        data.frame(coords, id = seq_len(nrow(coords))),
        geom = c("x", "y"),
        crs = ""
    )
}

# Write a parquetGeomStore then tile it.
# threshold = 4 forces multiple tiles over the 25-point grid.
make_tile_store <- function() {
    pgs  <- parquetGeomStore() |> storeWrite(make_tile_pts())
    pgts <- parquetGeomTileStore()
    storeWrite(pgts, pgs, threshold = 4L)
}

# tile pruning via crop ####

describe("crop tile pruning", {

    test_that("@tile_filter is empty on a fresh tile store", {
        s <- make_tile_store()
        expect_equal(s@tile_filter, integer(0L))
    })

    test_that("crop populates @tile_filter", {
        s <- make_tile_store()
        s2 <- crop(s, terra::ext(0, 25, 0, 25))
        expect_gt(length(s2@tile_filter), 0L)
    })

    test_that("crop @tile_filter is a subset of all tile indices", {
        s <- make_tile_store()
        all_tiles <- s@tiles$tile  # all tile indices in the plan
        s2 <- crop(s, terra::ext(0, 25, 0, 25))
        expect_true(all(s2@tile_filter %in% all_tiles))
        expect_lt(length(s2@tile_filter), length(all_tiles))
    })

    test_that("cropped store returns only rows within the crop extent", {
        s  <- make_tile_store()
        s2 <- crop(s, terra::ext(0, 25, 0, 25))
        tbl <- storeRead(s2, output = "tibble")
        expect_gt(nrow(tbl), 0L)
        expect_true(all(tbl$x_index <= 25))
        expect_true(all(tbl$y_index <= 25))
        # verify rows outside crop are excluded
        expect_equal(nrow(storeRead(s, output = "tibble")), 25L)
        expect_lt(nrow(tbl), 25L)
    })

    test_that("second crop further restricts @tile_filter", {
        s  <- make_tile_store()
        s2 <- crop(s,  terra::ext(0, 40, 0, 40))
        s3 <- crop(s2, terra::ext(0, 20, 0, 20))
        expect_lte(length(s3@tile_filter), length(s2@tile_filter))
    })

    test_that("crop outside data extent errors", {
        s <- make_tile_store()
        expect_error(crop(s, terra::ext(200, 300, 200, 300)))
    })

    test_that("explicit tile_idx in storeRead overrides @tile_filter", {
        s  <- make_tile_store()
        s2 <- crop(s, terra::ext(0, 25, 0, 25))
        # supplying tile_idx = integer(0) explicitly returns empty result
        # (overrides the @tile_filter which would return rows)
        tbl_explicit <- storeRead(s2, tile_idx = integer(0L), output = "tibble")
        expect_equal(nrow(tbl_explicit), 0L)
    })

    test_that("uncropped storeRead with tile_idx still works", {
        s <- make_tile_store()
        all_tbl  <- storeRead(s,  output = "tibble")
        tile1_tbl <- storeRead(s, tile_idx = s@tiles$tile[[1L]], output = "tibble")
        expect_gt(nrow(all_tbl), nrow(tile1_tbl))
    })

})
