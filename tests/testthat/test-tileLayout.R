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

# SQL readers over hive discovery ####

.layout_tile_store <- function() {
    pts <- terra::vect(.layout_pts_df(), geom = c("x", "y"), crs = "")
    pgs <- parquetGeomStore() |> storeWrite(pts)
    parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
}

# Rewrite a tile store into the pre-#74 layout: every tile's files moved
# down into a second `tile_index=000/` level.
.nest_tile_store <- function(store) {
    src <- file.path(storePaths(store), paste0("source_id=", store@uid))
    for (d in list.files(src, pattern = "^tile_index=", full.names = TRUE)) {
        leaf <- file.path(d, "tile_index=000")
        dir.create(leaf)
        files <- list.files(d, pattern = "\\.parquet$", full.names = TRUE)
        file.rename(files, file.path(leaf, basename(files)))
    }
    store
}

.sql_read <- function(store, engine, ...) {
    out <- storeRead(store, output = engine, ...)
    as.data.frame(if (engine == "sedona") sedonadb::sd_collect(out) else
        dplyr::collect(out))
}

.expected_tiles <- function(store) {
    tbl <- storeRead(store, output = "tibble", omit_internals = FALSE)
    sort(unique(tbl$tile_index))
}

describe("SQL readers take partition columns from hive discovery", {

    test_that("discovery is chosen unless the layout or sedonadb rules it out", {
        s <- .layout_tile_store()
        specs <- .pstore_tile_specs(s)
        expect_true(.pstore_sql_discovery(specs, "duckdb"))
        nested <- .pstore_tile_specs(.nest_tile_store(s))
        expect_true(all(vapply(nested, `[[`, logical(1L), "nested")))
        expect_false(.pstore_sql_discovery(nested, "duckdb"))
        skip_if_not_installed("sedonadb", minimum_version = "0.4.0")
        expect_true(.pstore_sql_discovery(specs, "sedona"))
        expect_false(.pstore_sql_discovery(nested, "sedona"))
    })

    for (engine in c("duckdb", "sedona")) {
        test_that(sprintf("%s: partition columns keep their types", engine), {
            if (engine == "duckdb") {
                skip_if_not_installed("duckdb")
                skip_if_not_installed("dbplyr")
            } else {
                skip_if_not_installed("sedonadb")
            }
            s <- .layout_tile_store()
            df <- .sql_read(s, engine)
            expect_type(df$source_id, "character")
            expect_type(df$tile_index, "integer")
            expect_setequal(unique(df$source_id), s@uid)
            expect_equal(sort(unique(df$tile_index)), .expected_tiles(s))
            expect_equal(nrow(df), nrow(s))
        })

        test_that(sprintf("%s: tile_idx narrows to that tile", engine), {
            if (engine == "duckdb") {
                skip_if_not_installed("duckdb")
                skip_if_not_installed("dbplyr")
            } else {
                skip_if_not_installed("sedonadb")
            }
            s <- .layout_tile_store()
            ti <- .expected_tiles(s)[2L]
            df <- .sql_read(s, engine, tile_idx = ti)
            expect_gt(nrow(df), 0L)
            expect_true(all(df$tile_index == ti))
        })

        test_that(sprintf("%s: reads a store in the pre-#74 nested layout", engine), {
            if (engine == "duckdb") {
                skip_if_not_installed("duckdb")
                skip_if_not_installed("dbplyr")
            } else {
                skip_if_not_installed("sedonadb")
            }
            s <- .layout_tile_store()
            tiles <- .expected_tiles(s)
            .nest_tile_store(s)
            df <- .sql_read(s, engine)
            expect_equal(nrow(df), nrow(s))
            expect_type(df$tile_index, "integer")
            expect_equal(sort(unique(df$tile_index)), tiles)
        })
    }
})
