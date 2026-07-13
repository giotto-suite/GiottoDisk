# Tests for the native duckdb compile path (`.pstore_to_duckdb`).
# Mirrors test-sedonadb.R in scope so the two SQL backends stay parallel
# in behavior. Skipped if {duckdb} isn't available.

make_pts_dd <- function(n = 5) {
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


# output type ####

test_that("duckdb: storeRead returns tbl_dbi", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd())
    tbl <- storeRead(pgs, output = "duckdb")
    expect_s3_class(tbl, "tbl_dbi")
})

test_that("duckdb: collect returns correct row count", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(5))
    df <- storeRead(pgs, output = "duckdb") |> dplyr::collect()
    expect_equal(nrow(df), 5L)
})


# connection lifecycle ####

test_that("duckdb: user-supplied connection is honored", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(3))
    conn <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    tbl <- storeRead(pgs, output = "duckdb", duckdb_params = list(conn = conn))
    expect_identical(tbl$src$con, conn)
})

test_that("duckdb: ephemeral connection is created when none supplied", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(3))
    tbl <- storeRead(pgs, output = "duckdb")
    # ephemeral conn lives on tbl$src$con and is valid for the duration
    df <- dplyr::collect(tbl)
    expect_equal(nrow(df), 3L)
})


# filtering ####

test_that("duckdb: crop reduces row count", {
    skip_if_not_installed("duckdb")
    pts <- terra::vect(
        data.frame(x = 1:10, y = 1:10),
        geom = c("x", "y"), crs = ""
    )
    pgs <- parquetGeomStore() |> storeWrite(pts)
    df <- storeRead(crop(pgs, terra::ext(3, 7, 3, 7)), output = "duckdb") |>
        dplyr::collect()
    expect_gt(nrow(df), 0L)
    expect_lt(nrow(df), 10L)
})

test_that("duckdb: @ops filter translates to SQL", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(5))
    pgs_f <- subset(pgs, cell_ID %in% c("cell_1", "cell_2"))
    df <- storeRead(pgs_f, output = "duckdb") |> dplyr::collect()
    expect_equal(nrow(df), 2L)
    expect_setequal(df$cell_ID, c("cell_1", "cell_2"))
})


# affine transform ####

test_that("duckdb: pending affine applies via ST_Affine", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(5))
    aff <- GiottoClass::affine(diag(c(2, 1)))  # scale x by 2, y unchanged
    pgs2 <- affine(pgs, aff)
    tbl <- storeRead(pgs2, output = "duckdb")
    # Pull coords via SQL rather than decoding the WKB list-column R-side.
    coords <- DBI::dbGetQuery(tbl$src$con, sprintf(
        'SELECT ST_X(geom) AS x_post, ST_Y(geom) AS y_post FROM "%s"',
        dbplyr::remote_name(tbl)
    ))
    expect_equal(sort(coords$x_post), seq_len(5) * 2, tolerance = 1e-6)
    expect_equal(sort(coords$y_post), seq_len(5),     tolerance = 1e-6)
})

# Regression guard for the duckdb ST_Affine branch (untouched by the
# sedona transpose fix). Shear catches any future coefficient swap.
test_that("duckdb: pending affine applies via ST_Affine (shear, off-diagonal)", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(5))
    # Row-vector post-mult M = [[1, 0.5], [0, 1]]: (n, n) -> (n, 1.5*n).
    aff <- GiottoClass::affine(matrix(c(1, 0, 0.5, 1), nrow = 2L))
    pgs2 <- affine(pgs, aff)
    tbl <- storeRead(pgs2, output = "duckdb")
    coords <- DBI::dbGetQuery(tbl$src$con, sprintf(
        'SELECT ST_X(geom) AS x_post, ST_Y(geom) AS y_post FROM "%s"',
        dbplyr::remote_name(tbl)
    ))
    expect_equal(sort(coords$x_post), as.numeric(seq_len(5)),       tolerance = 1e-6)
    expect_equal(sort(coords$y_post), as.numeric(seq_len(5)) * 1.5, tolerance = 1e-6)
})


# spatial predicates ####

test_that("duckdb: spat_relate intersects matches expected ids", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(5))
    df <- pgs |>
        spatRelate("POLYGON ((0 0, 3 0, 3 3, 0 3, 0 0))", "intersects") |>
        storeRead(output = "duckdb") |>
        dplyr::collect()
    expect_setequal(df$cell_ID, c("cell_1", "cell_2", "cell_3"))
})


# hive partition cols ####

test_that("duckdb: flat parquetStore injects source_id only (no tile_index)", {
    skip_if_not_installed("duckdb")
    ps <- parquetStore() |> storeWrite(mtcars)
    df <- storeRead(ps, output = "duckdb") |> dplyr::collect()
    expect_true("source_id" %in% names(df))
    expect_false("tile_index" %in% names(df))
    expect_equal(length(unique(df$source_id)), 1L)
})

test_that("duckdb: parquetGeomStore injects source_id and tile_index=0", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(5))
    df <- storeRead(pgs, output = "duckdb") |> dplyr::collect()
    expect_true(all(c("source_id", "tile_index") %in% names(df)))
    expect_equal(length(unique(df$source_id)), 1L)
    expect_equal(unique(df$tile_index), 0L)
})


.make_tile_pts_dd <- function() {
    coords <- data.frame(
        x = rep(seq(5, 45, by = 10), 5),
        y = rep(seq(5, 45, by = 10), each = 5)
    )
    terra::vect(
        data.frame(coords, id = seq_len(nrow(coords))),
        geom = c("x", "y"), crs = ""
    )
}

test_that("duckdb: parquetGeomTileStore unions all tiles, composite key is unique", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(.make_tile_pts_dd())
    pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
    df <- storeRead(pgts, output = "duckdb") |> dplyr::collect()
    expect_equal(nrow(df), 25L)
    expect_lt(length(unique(df$row_index)), 25L)
    composite <- paste(df$tile_index, df$row_index)
    expect_equal(length(unique(composite)), 25L)
})

test_that("duckdb: @tile_filter prunes at file level", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(.make_tile_pts_dd())
    pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
    pgts_cropped <- crop(pgts, terra::ext(0, 25, 0, 25))
    expect_gt(length(pgts_cropped@tile_filter), 0L)
    df <- storeRead(pgts_cropped, output = "duckdb") |> dplyr::collect()
    expect_lt(nrow(df), 25L)
    expect_true(all(unique(df$tile_index) %in% pgts_cropped@tile_filter))
})

test_that("duckdb: tile_idx arg overrides filter, selects only that tile", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(.make_tile_pts_dd())
    pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
    df <- storeRead(pgts, output = "duckdb", tile_idx = 1L) |> dplyr::collect()
    expect_equal(unique(df$tile_index), 1L)
})

test_that("duckdb: tile_idx with no matching tile errors clearly", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_dd(5))
    expect_error(
        storeRead(pgs, output = "duckdb", tile_idx = 999L),
        "no tile directories match"
    )
})


# fields narrowing ####

test_that("duckdb: [, j] column selection narrows projection", {
    skip_if_not_installed("duckdb")
    ps <- parquetStore() |> storeWrite(mtcars)
    df <- ps[, c("mpg", "cyl")] |>
        storeRead(output = "duckdb") |>
        dplyr::collect()
    expect_setequal(names(df), c("mpg", "cyl"))
})

test_that("duckdb: fields= arg narrows projection", {
    skip_if_not_installed("duckdb")
    ps <- parquetStore() |> storeWrite(mtcars)
    df <- storeRead(ps, output = "duckdb", fields = c("mpg", "hp")) |>
        dplyr::collect()
    expect_setequal(names(df), c("mpg", "hp"))
})

test_that("duckdb: col-select + filter on different col composes", {
    skip_if_not_installed("duckdb")
    ps <- parquetStore() |> storeWrite(mtcars)
    df <- ps[, "mpg"] |>
        subset(cyl == 4) |>
        storeRead(output = "duckdb") |>
        dplyr::collect()
    expect_setequal(names(df), "mpg")
    expect_equal(nrow(df), sum(mtcars$cyl == 4))
})
