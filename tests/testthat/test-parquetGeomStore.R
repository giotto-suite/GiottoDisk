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


# polygon split_geom / part_col ####

# SpatVector with two polygons: c1 has 2 parts (multipart), c2 has 1 part.
# Built via createGiottoPolygon's data.frame path with part_col.
make_multipart_sv <- function() {
    dt <- data.table::data.table(
        poly_ID = c(rep("c1", 4), rep("c1", 4), rep("c2", 4)),
        part = c(rep(1L, 4), rep(2L, 4), rep(3L, 4)),
        x = c(0, 1, 1, 0, 2, 3, 3, 2, 10, 11, 11, 10),
        y = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
    )
    GiottoClass::createGiottoPolygon(dt, part_col = "part", verbose = FALSE)[]
}

# Vertex-level polygon data: cell_id="c1" has 2 parts (label_id 1, 2),
# cell_id="c2" has 1 part (label_id 3). Used for part_col + tile tests.
make_multipart_vertex_dt <- function() {
    data.table::data.table(
        cell_id = c(rep("c1", 4), rep("c1", 4), rep("c2", 4)),
        label_id = c(rep(1L, 4), rep(2L, 4), rep(3L, 4)),
        vertex_x = c(0, 1, 1, 0, 2, 3, 3, 2, 10, 11, 11, 10),
        vertex_y = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
    )
}

test_that("parquetGeomStore: split_geom disaggs multipart polygons", {
    sv <- make_multipart_sv()
    pgs_flat <- parquetGeomStore() |> storeWrite(sv)
    pgs_split <- parquetGeomStore() |> storeWrite(sv, split_geom = TRUE)
    expect_equal(nrow(pgs_flat), 2L) # c1 multipart, c2 single
    expect_equal(nrow(pgs_split), 3L) # disagg: 2 from c1 + 1 from c2
})

test_that("parquetGeomStore: split_geom_sourcename preserves original poly_ID", {
    sv <- make_multipart_sv()
    pgs <- parquetGeomStore() |> storeWrite(
        sv,
        split_geom = TRUE,
        split_geom_fmt = "p_%d",
        split_geom_sourcename = "orig_id"
    )
    tbl <- storeRead(pgs, output = "tibble")
    expect_true("orig_id" %in% colnames(tbl))
    # c1 appears twice (was multipart), c2 once
    expect_equal(sort(tbl$orig_id), c("c1", "c1", "c2"))
})

test_that("parquetGeomStore: part_col groups multipart vertex rows", {
    # Without part_col, .df_to_terra_poly's unique() on meta would mismatch
    # (label_id varies within cell_id). With part_col = "label_id", c1's
    # vertices are grouped into a single multipart polygon.
    dt <- make_multipart_vertex_dt()
    pgs <- parquetGeomStore() |> storeWrite(
        dt,
        type = "polygons",
        id_col = "cell_id",
        sdimx = "vertex_x",
        sdimy = "vertex_y",
        part_col = "label_id"
    )
    expect_equal(nrow(pgs), 2L) # 2 polygons: c1 multipart, c2 single
})

test_that("parquetGeomTileStore: split_geom IDs unique and tile-tagged", {
    dt <- make_multipart_vertex_dt()
    ps <- parquetStore() |> storeWrite(dt)
    pgts <- parquetGeomTileStore() |> storeWrite(
        ps,
        type = "polygons",
        id_col = "cell_id", sdimx = "vertex_x", sdimy = "vertex_y",
        group_col = "cell_id", part_col = "label_id",
        split_geom = TRUE,
        split_geom_fmt = "nucleus_%d",
        split_geom_sourcename = "cell_poly_id"
    )
    tbl <- storeRead(pgts, output = "tibble")
    # 3 disaggregated polys
    expect_equal(nrow(tbl), 3L)
    # all unique IDs (tile_idx injection prevents cross-tile collisions)
    expect_equal(length(unique(tbl$poly_ID)), 3L)
    # IDs follow `nucleus_<tile_idx>_<n>` pattern
    expect_true(all(grepl("^nucleus_[0-9]+_[0-9]+$", tbl$poly_ID)))
    # cell_poly_id preserves the pre-split cell_ids
    expect_setequal(unique(tbl$cell_poly_id), c("c1", "c2"))
})

# points: coordinates only, WKB on request ####

.pt_files <- function(pgs) {
    list.files(storePaths(pgs), pattern = "\\.parquet$", recursive = TRUE,
        full.names = TRUE)
}

test_that("parquetGeomStore: points are written without WKB or GeoParquet", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    tab <- arrow::read_parquet(.pt_files(pgs)[[1L]], as_data_frame = FALSE)
    expect_false("geom" %in% names(tab))
    expect_null(tab$metadata$geo)
    expect_equal(geomtype(pgs), "points")
    sv <- storeRead(pgs, output = "terra")
    expect_equal(unname(terra::crds(sv)[order(sv$cell_ID), ]),
        cbind(as.numeric(1:5), as.numeric(1:5)))
    expect_s3_class(storeRead(pgs, output = "sf"), "sf")
})

test_that("parquetGeomStore: data.frame points are written without WKB", {
    df <- data.frame(feat_ID = letters[1:4], x = 1:4, y = 4:1)
    pgs <- parquetGeomStore() |>
        storeWrite(df, type = "points", id_col = "feat_ID",
            sdimx = "x", sdimy = "y")
    expect_false("geom" %in% names(arrow::read_parquet(.pt_files(pgs)[[1L]])))
    sv <- storeRead(pgs, output = "terra")
    expect_equal(nrow(sv), 4L)
})

test_that("parquetGeomStore: giottodisk.points_wkb writes a GeoParquet geom", {
    pgs <- GiottoUtils::gwith_options(list(giottodisk.points_wkb = TRUE),
        parquetGeomStore() |> storeWrite(make_pts(5)))
    tab <- arrow::read_parquet(.pt_files(pgs)[[1L]], as_data_frame = FALSE)
    expect_true("geom" %in% names(tab))
    expect_false(is.null(tab$metadata$geo))
    expect_equal(nrow(storeRead(pgs, output = "terra")), 5L)
})

test_that("parquetGeomStore: a points store reopened from its path is points", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    reopened <- parquetGeomStore(pgs@path)
    expect_equal(geomtype(reopened), "points")
    expect_equal(nrow(storeRead(reopened, output = "terra")), 5L)
    # a store that never went through initialize() after its write
    unread <- parquetGeomStore(tempfile())
    storeWrite(unread, make_pts(5))
    expect_equal(nrow(storeRead(unread, output = "terra")), 5L)
})

test_that("parquetGeomStore: points as.data.table builds XY and WKB from coords", {
    pgs <- parquetGeomStore() |> storeWrite(make_pts(3))
    xy <- data.table::as.data.table(pgs, geom = "XY")
    expect_true(all(c("x", "y", "cell_ID") %in% names(xy)))
    expect_equal(sort(xy$x), c(1, 2, 3))
    wkb <- data.table::as.data.table(pgs, geom = "wkb")
    expect_true("geom" %in% names(wkb))
    expect_equal(nrow(terra::vect(wkb$geom)), 3L)
})

test_that("parquetGeomStore: arrow reads of points never see an on-disk geom", {
    pgs <- GiottoUtils::gwith_options(list(giottodisk.points_wkb = TRUE),
        parquetGeomStore() |> storeWrite(make_pts(5)))
    expect_true("geom" %in% names(arrow::read_parquet(.pt_files(pgs)[[1L]])))
    expect_false("geom" %in% names(storeRead(pgs, output = "query")))
    expect_false("geom" %in% names(storeRead(pgs, output = "tibble")))
    # `geom` is virtual on points: nameable, never read from disk
    sel <- pgs[, c("cell_ID", "geom")]
    expect_false("geom" %in% names(storeRead(sel, output = "query")))
    expect_equal(nrow(storeRead(sel, output = "terra")), 5L)
})
