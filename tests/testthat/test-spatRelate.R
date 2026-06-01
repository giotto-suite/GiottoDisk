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


# Arrow path: spatial predicates are unsupported on arrow.
# A correct arrow implementation would tile + stream the predicate;
# the previous one-shot collect() approach was not safe at atlas scale.

test_that("spatRelate(): arrow storeRead errors with sedona nudge", {
    pgs <- .make_pts_store()
    expect_error(
        spatRelate(pgs, .roi(), "intersects") |>
            storeRead(output = "tibble"),
        "sedona"
    )
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
    sdf <- spatRelate(pgs, .roi(), "intersects") |>
        subset(id != "a") |>
        storeRead(output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_equal(df$id, "b")
})

test_that("spatRelate(): composes with [, j] narrowing", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    sdf <- spatRelate(pgs, .roi(), "intersects")[, "id"] |>
        storeRead(output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    # sedona carries internal index cols; what matters is the user col is
    # present and filters/values are correct.
    expect_true("id" %in% names(df))
    expect_setequal(df$id, c("a", "b"))
})

test_that("spatRelate(): composes with head()", {
    skip_if_not_installed("sedonadb")
    pgs <- .make_pts_store()
    sdf <- spatRelate(pgs, .roi(), "intersects") |>
        head(1L) |>
        storeRead(output = "sedona")
    df <- sedonadb::sd_collect(sdf)
    expect_equal(nrow(df), 1L)
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
    pre <- sedonadb::sd_collect(storeRead(s, output = "sedona"))

    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(s, tmp)
    rt <- readRDS(tmp)
    expect_length(rt@ops, 1L)
    expect_equal(rt@ops[[1L]]$type, "spat_relate")
    expect_equal(rt@ops[[1L]]$y_wkt, s@ops[[1L]]$y_wkt)
    post <- sedonadb::sd_collect(storeRead(rt, output = "sedona"))
    expect_setequal(pre$id, post$id)
})
