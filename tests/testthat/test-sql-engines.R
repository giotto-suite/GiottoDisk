# The two SQL read paths for parquet stores: duckdb (`.pstore_to_duckdb`) and
# sedonadb. The shared behavior is written once and run on each engine, so the
# two stay parallel by construction; what only one engine has is at the end.

.sql_engines <- list(
    duckdb = list(
        output = "duckdb",
        lazy_class = "tbl_dbi",
        skip = function() skip_if_not_installed("duckdb"),
        collect = function(x) dplyr::collect(x),
        # Pull coords via SQL rather than decoding the WKB list-column R-side.
        xy = function(x) {
            q <- DBI::dbGetQuery(x$src$con, sprintf(
                'SELECT ST_X(geom) AS x, ST_Y(geom) AS y FROM "%s"',
                dbplyr::remote_name(x)))
            cbind(x = q$x, y = q$y)
        }
    ),
    sedonadb = list(
        output = "sedona",
        lazy_class = "sedonadb_dataframe",
        skip = function() skip_if_not_installed("sedonadb",
            minimum_version = "0.4.0"),
        collect = function(x) sedonadb::sd_collect(x),
        xy = function(x) {
            df <- sedonadb::sd_collect(x)
            terra::crds(terra::vect(unclass(wk::as_wkb(df$geom))))
        }
    )
)

for (eng in names(.sql_engines)) {
    e <- .sql_engines[[eng]]
    .read <- function(x, ...) e$collect(storeRead(x, output = e$output, ...))


    # output type ####

    test_that(sprintf("%s: storeRead returns the engine's lazy table", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts())
        expect_s3_class(storeRead(pgs, output = e$output), e$lazy_class)
    })

    test_that(sprintf("%s: collect returns correct row count", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        expect_equal(nrow(.read(pgs)), 5L)
    })


    # filtering ####

    test_that(sprintf("%s: crop reduces row count", eng), {
        e$skip()
        pts <- terra::vect(
            data.frame(x = 1:10, y = 1:10),
            geom = c("x", "y"), crs = ""
        )
        pgs <- parquetGeomStore() |> storeWrite(pts)
        df <- .read(crop(pgs, terra::ext(3, 7, 3, 7)))
        expect_gt(nrow(df), 0L)
        expect_lt(nrow(df), 10L)
    })

    test_that(sprintf("%s: @ops filter translates to SQL", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        df <- .read(subset(pgs, cell_ID %in% c("cell_1", "cell_2")))
        expect_equal(nrow(df), 2L)
        expect_setequal(df$cell_ID, c("cell_1", "cell_2"))
    })

    test_that(sprintf("%s: %%in%% above the threshold uses a VALUES subquery", eng), {
        # DuckDB names a VALUES column `col0`, not DataFusion's `column1`.
        e$skip()
        old <- options(giottodisk.sedona_in_subquery_threshold = 1L)
        on.exit(options(old), add = TRUE)
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        df <- .read(subset(pgs, cell_ID %in% c("cell_1", "cell_2")))
        expect_setequal(df$cell_ID, c("cell_1", "cell_2"))
    })


    # affine transform ####

    test_that(sprintf("%s: pending affine applies via ST_Affine", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        aff <- GiottoClass::affine(diag(c(2, 1)))  # scale x by 2, y unchanged
        xy <- e$xy(storeRead(affine(pgs, aff), output = e$output))
        expect_equal(sort(xy[, "x"]), seq_len(5) * 2, tolerance = 1e-6)
        expect_equal(sort(xy[, "y"]), seq_len(5),     tolerance = 1e-6)
    })

    # Shear catches the ST_Affine sedona-vs-duckdb transpose. Diagonal
    # scale is transpose-invariant so needs off-diagonal to expose the bug.
    test_that(sprintf("%s: pending affine applies via ST_Affine (shear, off-diagonal)", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        # Row-vector post-mult M = [[1, 0.5], [0, 1]]:
        #   x' = x*M[1,1] + y*M[2,1] = x
        #   y' = x*M[1,2] + y*M[2,2] = 0.5*x + y
        # So (n, n) -> (n, 1.5*n). Buggy transpose would give (1.5*n, n).
        aff <- GiottoClass::affine(matrix(c(1, 0, 0.5, 1), nrow = 2L))
        xy <- e$xy(storeRead(affine(pgs, aff), output = e$output))
        expect_equal(sort(xy[, "x"]), as.numeric(seq_len(5)),       tolerance = 1e-6)
        expect_equal(sort(xy[, "y"]), as.numeric(seq_len(5)) * 1.5, tolerance = 1e-6)
    })


    # spatial predicates ####

    test_that(sprintf("%s: spat_relate intersects matches expected ids", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        df <- pgs |>
            spatRelate("POLYGON ((0 0, 3 0, 3 3, 0 3, 0 0))", "intersects") |>
            .read()
        expect_setequal(df$cell_ID, c("cell_1", "cell_2", "cell_3"))
    })


    # hive partition cols (source_id, tile_index) ####
    # These cols don't live in-file -- they're hive directory partitions. The
    # sedonadb R binding doesn't expose `table_partition_cols`, so the sedona
    # pipeline reconstructs them via per-tile-dir view registration + SQL
    # literal injection. Tests lock both the presence of the cols and the
    # downstream invariants (row_index disambiguation, file-level pruning).

    test_that(sprintf("%s: flat parquetStore injects source_id only (no tile_index)", eng), {
        e$skip()
        ps <- parquetStore() |> storeWrite(mtcars)
        df <- .read(ps)
        expect_true("source_id" %in% names(df))
        expect_false("tile_index" %in% names(df))
        expect_equal(length(unique(df$source_id)), 1L)
    })

    test_that(sprintf("%s: parquetGeomStore injects source_id and tile_index=0", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        df <- .read(pgs)
        expect_true(all(c("source_id", "tile_index") %in% names(df)))
        expect_equal(length(unique(df$source_id)), 1L)
        expect_equal(unique(df$tile_index), 0L)
    })

    test_that(sprintf("%s: parquetGeomTileStore unions all tiles, (tile_index, row_index) is unique", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_tile_pts())
        pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
        df <- .read(pgts)
        expect_equal(nrow(df), 25L)
        # row_index alone collides across tiles -- composite must disambiguate
        expect_lt(length(unique(df$row_index)), 25L)
        composite <- paste(df$tile_index, df$row_index)
        expect_equal(length(unique(composite)), 25L)
    })

    test_that(sprintf("%s: @tile_filter prunes at file level", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_tile_pts())
        pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
        pgts_cropped <- crop(pgts, terra::ext(0, 25, 0, 25))
        expect_gt(length(pgts_cropped@tile_filter), 0L)
        df <- .read(pgts_cropped)
        expect_lt(nrow(df), 25L)
        expect_true(all(unique(df$tile_index) %in% pgts_cropped@tile_filter))
    })

    test_that(sprintf("%s: tile_idx arg overrides filter, selects only that tile", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_tile_pts())
        pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
        df <- .read(pgts, tile_idx = 1L)
        expect_equal(unique(df$tile_index), 1L)
    })

    test_that(sprintf("%s: tile_idx with no matching tile errors clearly", eng), {
        e$skip()
        pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
        expect_error(
            storeRead(pgs, output = e$output, tile_idx = 999L),
            "no tile directories match"
        )
    })


    # fields narrowing ####
    # `[, j]` sets `@fields`, which each engine projects through SQL; a direct
    # `fields = ...` arg behaves the same.

    test_that(sprintf("%s: [, j] column selection narrows projection", eng), {
        e$skip()
        ps <- parquetStore() |> storeWrite(mtcars)
        df <- .read(ps[, c("mpg", "cyl")])
        expect_setequal(names(df), c("mpg", "cyl"))
    })

    test_that(sprintf("%s: fields= arg narrows projection", eng), {
        e$skip()
        ps <- parquetStore() |> storeWrite(mtcars)
        df <- .read(ps, fields = c("mpg", "hp"))
        expect_setequal(names(df), c("mpg", "hp"))
    })

    test_that(sprintf("%s: col-select + filter on different col composes", eng), {
        # Same composition as test-parquetStore.R but through the SQL path --
        # WHERE references a col not in SELECT; the FROM scope covers it.
        e$skip()
        ps <- parquetStore() |> storeWrite(mtcars)
        df <- .read(subset(ps[, "mpg"], cyl == 4))
        expect_setequal(names(df), "mpg")
        expect_equal(nrow(df), sum(mtcars$cyl == 4))
    })
}


# duckdb only: connection handling ####

test_that("duckdb: conn / name passed directly are rejected", {
    # Same guard as the expression stores -- `...` is not read by any duckdb
    # path, so a direct `conn` would silently give the caller a tbl_dbi on an
    # ephemeral connection instead of theirs.
    skip_if_not_installed("duckdb")
    conn <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    ps <- storeWrite(parquetStore(path = tempfile(fileext = ".parquet")),
        data.frame(a = 1:5, b = letters[1:5]))
    expect_error(storeRead(ps, output = "duckdb", conn = conn),
        "must be passed inside `duckdb_params`")
    expect_error(storeRead(ps, output = "duckdb", name = "x"),
        "must be passed inside `duckdb_params`")
})

test_that("duckdb: storeRead honours duckdb_params$name", {
    skip_if_not_installed("duckdb")
    conn <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    ps <- storeWrite(parquetStore(path = tempfile(fileext = ".parquet")),
        data.frame(a = 1:5, b = letters[1:5]))
    tbl <- storeRead(ps, output = "duckdb",
        duckdb_params = list(conn = conn, name = "ps_view"))
    expect_equal(dbplyr::remote_name(tbl), "ps_view")
    expect_equal(nrow(dplyr::collect(tbl)), 5L)
})

test_that("duckdb: user-supplied connection is honored", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts(3))
    conn <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    tbl <- storeRead(pgs, output = "duckdb", duckdb_params = list(conn = conn))
    expect_identical(tbl$src$con, conn)
})

test_that("duckdb: ephemeral connection is created when none supplied", {
    skip_if_not_installed("duckdb")
    pgs <- parquetGeomStore() |> storeWrite(make_pts(3))
    tbl <- storeRead(pgs, output = "duckdb")
    # ephemeral conn lives on tbl$src$con and is valid for the duration
    df <- dplyr::collect(tbl)
    expect_equal(nrow(df), 3L)
})


# sedonadb only: geoarrow output and the registered view ####

test_that("sedonadb: geom column is geoarrow_vctr after collect", {
    skip_if_not_installed("sedonadb", minimum_version = "0.4.0")
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    df <- sedonadb::sd_collect(storeRead(pgs, output = "sedona"))
    expect_true(inherits(df$geom, "geoarrow_vctr"))
})

test_that("sedonadb: sd_view_ref returns double-quoted lowercase name", {
    skip_if_not_installed("sedonadb", minimum_version = "0.4.0")
    pgs <- parquetGeomStore() |> storeWrite(make_pts())
    sdf <- storeRead(pgs, output = "sedona")
    ref <- sd_view_ref(sdf)
    expect_type(ref, "character")
    expect_match(ref, '^"gd_[a-z0-9_]+"$')
})

test_that("sedonadb: sd_view_ref errors on object without view_name attribute", {
    expect_error(sd_view_ref(list()), "view_name")
})

test_that("sedonadb: view is queryable via sd_sql", {
    skip_if_not_installed("sedonadb", minimum_version = "0.4.0")
    pgs <- parquetGeomStore() |> storeWrite(make_pts(5))
    sdf <- storeRead(pgs, output = "sedona")
    result <- sedonadb::sd_collect(
        sedonadb::sd_sql(sprintf("SELECT COUNT(*) AS n FROM %s", sd_view_ref(sdf)))
    )
    expect_equal(result$n, 5L)
})

# sedonadb only: CRS ####

test_that("sedonadb: a store without a CRS reads in the planar engineering frame", {
    skip_if_not_installed("sedonadb", minimum_version = "0.4.0")
    # two distinct shapes: sedona reads a geom column whose statistics have
    # min == max (one row, identical shapes) as plain binary
    sq <- terra::vect(c("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))",
        "POLYGON ((20 0, 30 0, 30 10, 20 10, 20 0))"), crs = "")
    sq$poly_ID <- c("a", "b")
    for (s in list(parquetGeomStore() |> storeWrite(make_pts(5)),
            parquetGeomStore() |> storeWrite(sq))) {
        sedonadb::sd_to_view(storeRead(s, output = "sedona"), "gd_crs_t",
            overwrite = TRUE)
        crs <- sedonadb::sd_collect(sedonadb::sd_sql(
            'SELECT ST_Crs(geom) AS crs FROM "gd_crs_t" LIMIT 1'))$crs
        expect_match(crs, "EngineeringCRS", fixed = TRUE)
    }
    # planar area, in coordinate units
    area <- expanse(parquetGeomStore() |> storeWrite(sq), engine = "sedona",
        output = "vector")
    expect_equal(area, c(100, 100))
})
