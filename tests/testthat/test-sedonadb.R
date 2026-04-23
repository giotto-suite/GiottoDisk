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
