# Hive layout of written geom stores ####
#
# A tile store is `source_id=<uid>/tile_index=<NNN>/part-*.parquet`, one
# level each. Flat geom stores take `tile_index=000` so they share the tile
# namespace; a tile of a tile store must not take it again (#74).

options("tilework.warn_sequential" = FALSE)

.tile_layout <- function(store) {
    files <- list.files(storePaths(store), pattern = "\\.parquet$",
        recursive = TRUE)
    sub("^source_id=[^/]+/", "", files)
}

.layout_pts_df <- function() {
    data.frame(
        x = rep(seq(5, 45, by = 10), 5),
        y = rep(seq(5, 45, by = 10), each = 5),
        feat = rep(c("a", "b", "c", "d", "e"), 5)
    )
}

.layout_poly_dt <- function() {
    # four unit squares in the corners of [0, 50]^2 so the tiler splits them
    sq <- function(id, x0, y0) {
        data.frame(
            cell_id = id,
            vertex_x = x0 + c(0, 1, 1, 0, 0),
            vertex_y = y0 + c(0, 0, 1, 1, 0)
        )
    }
    rbind(sq("c1", 2, 2), sq("c2", 47, 2), sq("c3", 2, 47), sq("c4", 47, 47))
}

.expect_tile_layout <- function(store) {
    layout <- .tile_layout(store)
    expect_gt(length(layout), 1L) # threshold forces several tiles
    expect_false(any(grepl("tile_index=[0-9]+/tile_index=", layout)))
    expect_true(all(grepl("^tile_index=[0-9]+/[^/]+\\.parquet$", layout)))
}

.expect_duckdb_reads <- function(store) {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("dbplyr")
    tbl <- storeRead(store, output = "duckdb")
    expect_equal(nrow(dplyr::collect(tbl)), nrow(store))
}

describe("tile writers do not nest a second tile_index level", {

    test_that("parquetGeomTileStore <- parquetGeomStore", {
        pts <- terra::vect(.layout_pts_df(), geom = c("x", "y"), crs = "")
        pgs <- parquetGeomStore() |> storeWrite(pts)
        pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
        .expect_tile_layout(pgts)
        .expect_duckdb_reads(pgts)
    })

    test_that("parquetGeomTileStore <- queryableStore", {
        f <- tempfile(fileext = ".parquet")
        arrow::write_parquet(.layout_pts_df(), f)
        qs <- as(fileStore(path = f,
            read_fun = function(x, ...) arrow::open_dataset(x)),
            "queryableStore")
        pgts <- parquetGeomTileStore() |> storeWrite(qs,
            type = "points", id_col = "feat", sdimx = "x", sdimy = "y",
            threshold = 4L)
        .expect_tile_layout(pgts)
        .expect_duckdb_reads(pgts)
    })

    test_that("parquetGeomTileStore <- parquetStore", {
        ps <- parquetStore() |> storeWrite(.layout_poly_dt())
        pgts <- parquetGeomTileStore() |> storeWrite(ps,
            type = "polygons", id_col = "cell_id",
            sdimx = "vertex_x", sdimy = "vertex_y",
            group_col = "cell_id", threshold = 1L)
        .expect_tile_layout(pgts)
        .expect_duckdb_reads(pgts)
    })
})

test_that("flat parquetGeomStore keeps its tile_index=000 level", {
    pts <- terra::vect(.layout_pts_df(), geom = c("x", "y"), crs = "")
    pgs <- parquetGeomStore() |> storeWrite(pts)
    expect_equal(.tile_layout(pgs), "tile_index=000/part-0.parquet")
})

test_that(".write_parquet refuses to nest a tile_index level", {
    tile_store <- storeCreate(
        path = file.path(tempfile(), "source_id=x", "tile_index=001"),
        type = "parquetGeom"
    )
    expect_error(
        .write_parquet(tile_store, .layout_pts_df(), uid_partition = FALSE,
            tile_idx = 0L),
        "already a tile_index partition"
    )
})
