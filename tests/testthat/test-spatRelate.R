# Helpers ####

.make_pts_store <- function() {
    pts <- terra::vect(
        data.frame(
            x = c(1, 3, 5, 7, 9),
            y = c(1, 3, 5, 7, 9),
            id = letters[1:5]
        ),
        geom = c("x", "y"),
        crs = ""
    )
    parquetGeomStore() |> storeWrite(pts)
}

.roi <- function(wkt = "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))") {
    terra::vect(wkt, crs = "")
}


# Queue contract ####
# spatRelate() queues a "spat_relate" op without I/O, then storeRead materializes.

test_that("spatRelate(): SpatVector input queues spat_relate op", {
    pgs <- .make_pts_store()
    s <- spatRelate(pgs, .roi(), relation = "intersects")
    expect_length(s@ops, 1L)
    expect_equal(s@ops[[1L]]$type, "spat_relate")
    expect_equal(s@ops[[1L]]$relation, "intersects")
    expect_equal(s@ops[[1L]]$form, "filter")
    expect_true(is.character(s@ops[[1L]]$y_wkt))
    expect_null(s@ops[[1L]]$y_store)
    # original store unaffected
    expect_length(pgs@ops, 0L)
})

test_that("spatRelate(): character WKT input is canonical, no normalization", {
    pgs <- .make_pts_store()
    wkt <- "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))"
    s <- spatRelate(pgs, wkt, relation = "intersects")
    expect_equal(s@ops[[1L]]$y_wkt, wkt)
})

test_that("spatRelate(): rejects invalid relation", {
    pgs <- .make_pts_store()
    expect_error(spatRelate(pgs, .roi(), relation = "nonsense"))
})

test_that("spatRelate(): rejects empty WKT", {
    pgs <- .make_pts_store()
    expect_error(spatRelate(pgs, ""))
})


# Arrow path: spat_relate routes through sedonadb internally to evaluate
# the predicate, then narrows the arrow query with a semi_join on
# surviving row_index (+ tile_index for tile stores). The rest of the
# chain stays lazy on the arrow side.

test_that("spatRelate(): arrow path -- intersects matches expected ids", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    tbl <- spatRelate(pgs, .roi(), "intersects") |>
        storeRead(output = "tibble")
    expect_setequal(tbl$id, c("a", "b"))
})

test_that("spatRelate(): arrow path -- each predicate produces sane result", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    roi <- .roi()
    results <- lapply(
        c("intersects", "within", "contains", "disjoint"),
        function(rel) {
            tbl <- spatRelate(pgs, roi, rel) |> storeRead(output = "tibble")
            sort(tbl$id)
        }
    )
    expect_setequal(results[[1L]], c("a", "b"))  # intersects
    expect_setequal(results[[2L]], c("a", "b"))  # within
    expect_setequal(results[[3L]], character(0L))  # contains
    expect_setequal(results[[4L]], c("c", "d", "e"))  # disjoint
})

test_that("spatRelate(): arrow path -- two spat_relate ops compose", {
    # Caching path: the second spat_relate's trimmed-store evaluation
    # narrows via the cached ids from the first (as an `id_filter` op),
    # not re-evaluating the first's spatial predicate.
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    roi1 <- .roi("POLYGON ((0 0, 6 0, 6 6, 0 6, 0 0))")  # hits a, b, c
    roi2 <- .roi("POLYGON ((2 2, 8 2, 8 8, 2 8, 2 2))")  # hits b, c, d
    tbl <- spatRelate(pgs, roi1, "intersects") |>
        spatRelate(roi2, "intersects") |>
        storeRead(output = "tibble")
    expect_setequal(tbl$id, c("b", "c"))
})

test_that("spatRelate(): arrow path -- three spat_relate ops compose", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    roi1 <- .roi("POLYGON ((0 0, 8 0, 8 8, 0 8, 0 0))")  # hits a, b, c, d
    roi2 <- .roi("POLYGON ((2 2, 10 2, 10 10, 2 10, 2 2))")  # hits b..e
    roi3 <- .roi("POLYGON ((4 4, 8 4, 8 8, 4 8, 4 4))")  # hits c, d
    tbl <- spatRelate(pgs, roi1, "intersects") |>
        spatRelate(roi2, "intersects") |>
        spatRelate(roi3, "intersects") |>
        storeRead(output = "tibble")
    expect_setequal(tbl$id, c("c", "d"))  # intersection of all three
})

test_that("spatRelate(): arrow path -- spat_relate interleaved with filter", {
    # An attribute filter sitting between two spat_relate ops shouldn't
    # break the cache or change the row set.
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    roi1 <- .roi("POLYGON ((0 0, 6 0, 6 6, 0 6, 0 0))")  # hits a, b, c
    roi2 <- .roi("POLYGON ((2 2, 8 2, 8 8, 2 8, 2 2))")  # hits b, c, d
    tbl <- spatRelate(pgs, roi1, "intersects") |>
        subset(id != "c") |>
        spatRelate(roi2, "intersects") |>
        storeRead(output = "tibble")
    expect_setequal(tbl$id, "b")
})

test_that("spatRelate(): arrow path -- different predicates compose", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    big_box <- .roi("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")  # all 5
    small_box <- .roi("POLYGON ((2 2, 4 2, 4 4, 2 4, 2 2))")    # contains b only
    # intersects(big) AND disjoint(small) = all 5 minus b = a, c, d, e
    tbl <- spatRelate(pgs, big_box, "intersects") |>
        spatRelate(small_box, "disjoint") |>
        storeRead(output = "tibble")
    expect_setequal(tbl$id, c("a", "c", "d", "e"))
})

test_that("spatRelate(): arrow + sedona paths agree on two-op chain", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    roi1 <- .roi("POLYGON ((0 0, 6 0, 6 6, 0 6, 0 0))")
    roi2 <- .roi("POLYGON ((2 2, 8 2, 8 8, 2 8, 2 2))")
    arrow_tbl <- spatRelate(pgs, roi1, "intersects") |>
        spatRelate(roi2, "intersects") |>
        storeRead(output = "tibble")
    sedona_tbl <- sedonadb::sd_collect(
        spatRelate(pgs, roi1, "intersects") |>
            spatRelate(roi2, "intersects") |>
            storeRead(output = "sedona")
    )
    expect_setequal(arrow_tbl$id, sedona_tbl$id)
})


# Sedona path correctness ####
# Five points on the diagonal at (1,1), (3,3), (5,5), (7,7), (9,9).
# ROI is the square (0,0)-(4,4). Expected hits: a, b.

test_that("spatRelate(): sedona path -- intersects matches expected ids", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    sdf <- spatRelate(pgs, .roi(), "intersects") |>
        storeRead(output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_setequal(df$id, c("a", "b"))
})

test_that("spatRelate(): sedona path -- each predicate produces sane result", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    roi <- .roi()
    results <- lapply(
        c("intersects", "within", "contains", "disjoint"),
        function(rel) {
            sdf <- spatRelate(pgs, roi, rel) |> storeRead(output = "sedona")
            sort(sedonadb::sd_collect(sdf)$id)
        }
    )
    expect_setequal(results[[1L]], c("a", "b"))  # intersects
    expect_setequal(results[[2L]], c("a", "b"))  # within
    expect_setequal(results[[3L]], character(0L))  # contains
    expect_setequal(results[[4L]], c("c", "d", "e"))  # disjoint
})


# Composition ####
# spat_relate composes with the existing op queue (subset, [, j], head, etc.).

test_that("spatRelate(): composes with subset() on a different col", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    tbl <- spatRelate(pgs, .roi(), "intersects") |>
        subset(id != "a") |>
        storeRead(output = "tibble")
    expect_equal(tbl$id, "b")
})

test_that("spatRelate(): composes with [, j] narrowing", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    tbl <- spatRelate(pgs, .roi(), "intersects")[, "id"] |>
        storeRead(output = "tibble")
    expect_setequal(names(tbl), "id")
    expect_setequal(tbl$id, c("a", "b"))
})

test_that("spatRelate(): composes with head()", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    tbl <- spatRelate(pgs, .roi(), "intersects") |>
        head(1L) |>
        storeRead(output = "tibble")
    expect_equal(nrow(tbl), 1L)
})


# Input-type dispatch ####

test_that("spatRelate(): single-feature SpatVector is passed through", {
    pgs <- .make_pts_store()
    s <- spatRelate(pgs, .roi(), "intersects")
    # ROI is a single-feature poly -- WKT should be the same input string
    expect_match(s@ops[[1L]]$y_wkt, "^POLYGON")
})

test_that("spatRelate(): multi-feature SpatVector is unioned (single WKT)", {
    pgs <- .make_pts_store()
    # Two ROIs -- terra::aggregate should union them into a MULTIPOLYGON
    multi_roi <- terra::vect(c(
        "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))",
        "POLYGON ((10 10, 12 10, 12 12, 10 12, 10 10))"
    ), crs = "")
    s <- spatRelate(pgs, multi_roi, "intersects")
    expect_length(s@ops, 1L)
    expect_match(s@ops[[1L]]$y_wkt, "MULTIPOLYGON|GEOMETRYCOLLECTION")
})

test_that("spatRelate(): >1000-feature SpatVector errors with nudge", {
    pgs <- .make_pts_store()
    big <- terra::vect(
        data.frame(x = seq_len(1500), y = 1, id = seq_len(1500)),
        geom = c("x", "y"), crs = ""
    )
    expect_error(
        spatRelate(pgs, big, "intersects"),
        "parquetGeomStore"
    )
})

test_that("spatRelate(): store/store filter form queues with y_store ref", {
    pgs <- .make_pts_store()
    # Build a small y store from a single ROI poly
    roi_sv <- .roi()
    y_store <- parquetGeomStore() |> storeWrite(roi_sv)
    s <- spatRelate(pgs, y_store, relation = "intersects", form = "filter")
    expect_length(s@ops, 1L)
    expect_equal(s@ops[[1L]]$form, "filter")
    expect_null(s@ops[[1L]]$y_wkt)
    expect_true(inherits(s@ops[[1L]]$y_store, "parquetGeomBase"))
})

test_that("spatRelate(): store/store form='join' errors (Phase 5)", {
    pgs <- .make_pts_store()
    y_store <- parquetGeomStore() |> storeWrite(.roi())
    expect_error(
        spatRelate(pgs, y_store, relation = "intersects", form = "join"),
        "not yet implemented"
    )
})


# Rejection: in-mem x against stored y ####

test_that("spatRelate(): (SpatVector, parquetGeomBase) is rejected", {
    pgs <- .make_pts_store()
    sv <- .roi()
    expect_error(spatRelate(sv, pgs, "intersects"), "stored y is unsupported")
})


# saveRDS roundtrip ####

test_that("spatRelate(): saveRDS roundtrip preserves op + result", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    s <- spatRelate(pgs, .roi(), "intersects")
    pre <- storeRead(s, output = "tibble")

    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(s, tmp)
    rt <- readRDS(tmp)
    expect_length(rt@ops, 1L)
    expect_equal(rt@ops[[1L]]$type, "spat_relate")
    expect_equal(rt@ops[[1L]]$y_wkt, s@ops[[1L]]$y_wkt)
    post <- storeRead(rt, output = "tibble")
    expect_setequal(pre$id, post$id)
})


# Engine selection ####
# The narrow path resolves an engine per call:
#   1. `giottodisk.spatial_query_engine` option (specific or "auto")
#   2. auto: sedona > duckdb > terra
# Tests pin via GiottoUtils::gwith_options so they don't depend on what's
# globally installed.

test_that("spatRelate(): explicit duckdb engine produces same ids as sedona", {
    skip_if_not_installed("sedonadb")
    skip_if_not_installed("duckdb")
    pgs <- .make_pts_store()
    sedona_tbl <- GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = "sedona"),
        spatRelate(pgs, .roi(), "intersects") |> storeRead(output = "tibble")
    )
    duckdb_tbl <- GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = "duckdb"),
        spatRelate(pgs, .roi(), "intersects") |> storeRead(output = "tibble")
    )
    expect_setequal(sedona_tbl$id, duckdb_tbl$id)
})

test_that("spatRelate(): duckdb engine handles two-op spat_relate chain", {
    # Cached id_filter handling works against the duckdb compile path too.
    skip_if_not_installed("duckdb")
    pgs <- .make_pts_store()
    roi1 <- .roi("POLYGON ((0 0, 6 0, 6 6, 0 6, 0 0))")
    roi2 <- .roi("POLYGON ((2 2, 8 2, 8 8, 2 8, 2 2))")
    tbl <- GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = "duckdb"),
        spatRelate(pgs, roi1, "intersects") |>
            spatRelate(roi2, "intersects") |>
            storeRead(output = "tibble")
    )
    expect_setequal(tbl$id, c("b", "c"))
})

test_that("spatRelate(): terra engine produces correct ids", {
    pgs <- .make_pts_store()
    tbl <- GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = "terra"),
        spatRelate(pgs, .roi(), "intersects") |> storeRead(output = "tibble")
    )
    expect_setequal(tbl$id, c("a", "b"))
})

test_that("spatRelate(): terra engine handles two-op spat_relate chain", {
    pgs <- .make_pts_store()
    roi1 <- .roi("POLYGON ((0 0, 6 0, 6 6, 0 6, 0 0))")  # a, b, c
    roi2 <- .roi("POLYGON ((2 2, 8 2, 8 8, 2 8, 2 2))")  # b, c, d
    tbl <- GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = "terra"),
        spatRelate(pgs, roi1, "intersects") |>
            spatRelate(roi2, "intersects") |>
            storeRead(output = "tibble")
    )
    expect_setequal(tbl$id, c("b", "c"))
})

test_that("spatRelate(): terra engine streams parquetGeomTileStore via tileApply", {
    # Tile-store path: terra engine uses tilework::tileApply per-tile
    # rather than materializing the whole trim. Verify correctness on a
    # multi-tile fixture.
    coords <- data.frame(
        x = rep(seq(5, 45, by = 10), 5),
        y = rep(seq(5, 45, by = 10), each = 5),
        id = sprintf("p%02d", seq_len(25))
    )
    pts <- terra::vect(coords, geom = c("x", "y"), crs = "")
    pgs <- parquetGeomStore() |> storeWrite(pts)
    pgts <- parquetGeomTileStore() |> storeWrite(pgs, threshold = 4L)
    # Quarter-box: 0..25 x 0..25 -- expect points whose centroid lies
    # in that AABB (rows 1..3 in x, rows 1..3 in y after the binning).
    tbl <- GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = "terra"),
        spatRelate(pgts,
            "POLYGON ((0 0, 25 0, 25 25, 0 25, 0 0))",
            "intersects") |>
            storeRead(output = "tibble")
    )
    # Centroids at multiples of 10 starting at 5 -- (5,5), (15,5), (5,15),
    # (15,15), (25,5), (5,25), (25,15), (15,25), (25,25) all intersect
    # the closed [0,25] box; the points at 35+ on either axis do not.
    expect_gt(nrow(tbl), 0L)
    expect_lt(nrow(tbl), 25L)
})

test_that("spatRelate(): terra/sedona/duckdb engines agree on results", {
    skip_if_not_installed("sedonadb")
    skip_if_not_installed("duckdb")
    pgs <- .make_pts_store()
    by_engine <- function(eng) {
        sort(GiottoUtils::gwith_options(
            list(giottodisk.spatial_query_engine = eng),
            spatRelate(pgs, .roi(), "intersects") |>
                storeRead(output = "tibble"))$id)
    }
    expect_equal(by_engine("terra"),  c("a", "b"))
    expect_equal(by_engine("sedona"), c("a", "b"))
    expect_equal(by_engine("duckdb"), c("a", "b"))
})

test_that("spatRelate(): auto fallback to terra nudges via inform", {
    pgs <- .make_pts_store()
    # Force the resolver to fall through to terra by mocking the
    # availability check. local_mocked_bindings is testthat 3.
    testthat::local_mocked_bindings(
        .spat_engine_available = function(pkg) FALSE
    )
    expect_message(
        GiottoUtils::gwith_options(
            list(giottodisk.spatial_query_engine = "auto"),
            spatRelate(pgs, .roi(), "intersects") |>
                storeRead(output = "tibble")
        ),
        "sedonadb"
    )
})

test_that("spatRelate(): auto engine resolves to an installed backend", {
    skip_if_not(
        requireNamespace("sedonadb", quietly = TRUE) ||
        requireNamespace("duckdb", quietly = TRUE),
        "no SQL spatial engine installed"
    )
    pgs <- .make_pts_store()
    # Default (auto) should produce results, not error out
    tbl <- GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = "auto"),
        spatRelate(pgs, .roi(), "intersects") |> storeRead(output = "tibble")
    )
    expect_setequal(tbl$id, c("a", "b"))
})
