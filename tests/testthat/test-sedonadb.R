make_pts_sdb <- function(n = 5) {
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

test_that("sedonadb: storeRead returns sedonadb_dataframe", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb())
    sdf <- storeRead(pgs, output = "sedona")
    expect_s3_class(sdf, "sedonadb_dataframe")
})

test_that("sedonadb: collect returns correct row count", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    sdf <- storeRead(pgs, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_equal(nrow(df), 5L)
})

test_that("sedonadb: geom column is geoarrow_vctr after collect", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    sdf <- storeRead(pgs, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_true(inherits(df$geom, "geoarrow_vctr"))
})

# sd_view_ref ####

test_that("sedonadb: sd_view_ref returns double-quoted lowercase name", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb())
    sdf <- storeRead(pgs, output = "sedona")
    ref <- sd_view_ref(sdf)
    expect_type(ref, "character")
    expect_match(ref, '^"gd_[a-z0-9_]+"$')
})

test_that("sedonadb: sd_view_ref errors on object without view_name attribute", {
    expect_error(sd_view_ref(list()), "view_name")
})

test_that("sedonadb: view is queryable via sd_sql", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    sdf <- storeRead(pgs, output = "sedona")
    result <- sedonadb::sd_collect(
        sedonadb::sd_sql(sprintf("SELECT COUNT(*) AS n FROM %s", sd_view_ref(sdf)))
    )
    expect_equal(result$n, 5L)
})

# filtering ####

test_that("sedonadb: crop reduces row count", {
    skip_if_not_installed("sedonadb")
    pts <- terra::vect(
        data.frame(x = 1:10, y = 1:10),
        geom = c("x", "y"), crs = ""
    )
    pgs <- parquetGeomStore() |> storeWrite(pts)
    sdf <- storeRead(crop(pgs, terra::ext(3, 7, 3, 7)), output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_gt(nrow(df), 0L)
    expect_lt(nrow(df), 10L)
})

test_that("sedonadb: @ops filter translates to SQL", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    pgs_f <- subset(pgs, cell_ID %in% c("cell_1", "cell_2"))
    sdf <- storeRead(pgs_f, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_equal(nrow(df), 2L)
    expect_setequal(df$cell_ID, c("cell_1", "cell_2"))
})

# affine transform ####

test_that("sedonadb: pending affine applies via ST_Affine", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    aff <- GiottoClass::affine(diag(c(2, 1)))  # scale x by 2, y unchanged
    pgs2 <- affine(pgs, aff)
    sdf <- storeRead(pgs2, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    sv <- terra::vect(unclass(wk::as_wkb(df$geom)))
    coords <- terra::crds(sv)
    expect_equal(sort(coords[, "x"]), seq_len(5) * 2, tolerance = 1e-6)
    expect_equal(sort(coords[, "y"]), seq_len(5),     tolerance = 1e-6)
})

# Shear catches the ST_Affine argument-order transpose: sedona uses
# (a, b, d, e) with x' = a*x + d*y, PostGIS/duckdb uses the transposed
# convention. Diagonal scale matrices are transpose-invariant so the
# scale test above passes with either wiring; only off-diagonal entries
# expose the bug.
test_that("sedonadb: pending affine applies via ST_Affine (shear, off-diagonal)", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    # Row-vector post-mult M = [[1, 0.5], [0, 1]]:
    #   x' = x*M[1,1] + y*M[2,1] = x
    #   y' = x*M[1,2] + y*M[2,2] = 0.5*x + y
    # So (n, n) -> (n, 1.5*n). Buggy transpose would give (1.5*n, n).
    aff <- GiottoClass::affine(matrix(c(1, 0, 0.5, 1), nrow = 2L))
    pgs2 <- affine(pgs, aff)
    sdf <- storeRead(pgs2, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    sv <- terra::vect(unclass(wk::as_wkb(df$geom)))
    coords <- terra::crds(sv)
    expect_equal(sort(coords[, "x"]), as.numeric(seq_len(5)),       tolerance = 1e-6)
    expect_equal(sort(coords[, "y"]), as.numeric(seq_len(5)) * 1.5, tolerance = 1e-6)
})


# hive partition cols (source_id, tile_index) ####
# These cols don't live in-file — they're hive directory partitions. The
# sedonadb R binding doesn't expose `table_partition_cols`, so the sedona
# pipeline reconstructs them via per-tile-dir view registration + SQL literal
# injection. Tests lock both the presence of the reconstructed cols and the
# downstream invariants (row_index disambiguation, file-level pruning).

.make_tile_pts_sdb <- function() {
    coords <- data.frame(
        x = rep(seq(5, 45, by = 10), 5),
        y = rep(seq(5, 45, by = 10), each = 5)
    )
    terra::vect(
        data.frame(coords, id = seq_len(nrow(coords))),
        geom = c("x", "y"), crs = ""
    )
}

test_that("sedonadb: flat parquetStore injects source_id only (no tile_index)", {
    skip_if_not_installed("sedonadb")
    ps <- parquetStore() |> storeWrite(mtcars)
    sdf <- storeRead(ps, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_true("source_id" %in% names(df))
    expect_false("tile_index" %in% names(df))
    expect_equal(length(unique(df$source_id)), 1L)
})

test_that("sedonadb: parquetGeomStore injects source_id and tile_index=0", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    sdf <- storeRead(pgs, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_true(all(c("source_id", "tile_index") %in% names(df)))
    expect_equal(length(unique(df$source_id)), 1L)
    expect_equal(unique(df$tile_index), 0L)
})

test_that("sedonadb: parquetGeomTileStore unions all tiles, (tile_index, row_index) is unique", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(.make_tile_pts_sdb())
    pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
    sdf <- storeRead(pgts, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_equal(nrow(df), 25L)
    # row_index alone collides across tiles -- composite must disambiguate
    expect_lt(length(unique(df$row_index)), 25L)
    composite <- paste(df$tile_index, df$row_index)
    expect_equal(length(unique(composite)), 25L)
})

test_that("sedonadb: @tile_filter prunes at file level", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(.make_tile_pts_sdb())
    pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
    pgts_cropped <- crop(pgts, terra::ext(0, 25, 0, 25))
    expect_gt(length(pgts_cropped@tile_filter), 0L)
    sdf <- storeRead(pgts_cropped, output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_lt(nrow(df), 25L)
    expect_true(all(unique(df$tile_index) %in% pgts_cropped@tile_filter))
})

test_that("sedonadb: tile_idx arg overrides filter, selects only that tile", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(.make_tile_pts_sdb())
    pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
    sdf <- storeRead(pgts, output = "sedona", tile_idx = 1L)
    df <- sedonadb::sd_collect(sdf)
    expect_equal(unique(df$tile_index), 1L)
})

test_that("sedonadb: tile_idx with no matching tile errors clearly", {
    skip_if_not_installed("sedonadb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts_sdb(5))
    expect_error(
        storeRead(pgs, output = "sedona", tile_idx = 999L),
        "no tile directories match"
    )
})


# fields narrowing via sedona ####
# `[, j]` sets `@fields` which the sedona pipeline projects through SQL;
# direct `fields = ...` arg behaves the same.

test_that("sedonadb: [, j] column selection narrows projection", {
    skip_if_not_installed("sedonadb")
    ps <- parquetStore() |> storeWrite(mtcars)
    sdf <- ps[, c("mpg", "cyl")] |> storeRead(output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_setequal(names(df), c("mpg", "cyl"))
})

test_that("sedonadb: fields= arg narrows projection", {
    skip_if_not_installed("sedonadb")
    ps <- parquetStore() |> storeWrite(mtcars)
    sdf <- storeRead(ps, output = "sedona", fields = c("mpg", "hp"))
    df <- sedonadb::sd_collect(sdf)
    expect_setequal(names(df), c("mpg", "hp"))
})

test_that("sedonadb: col-select + filter on different col composes", {
    # Same composition as test-parquetStore.R but through the SQL path --
    # WHERE references a col not in SELECT; sedona's FROM scope covers it.
    skip_if_not_installed("sedonadb")
    ps <- parquetStore() |> storeWrite(mtcars)
    sdf <- ps[, "mpg"] |> subset(cyl == 4) |> storeRead(output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_setequal(names(df), "mpg")
    expect_equal(nrow(df), sum(mtcars$cyl == 4))
})
