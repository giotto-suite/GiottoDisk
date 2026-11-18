# End-to-end CosMx disk reader tests. Fixture: helper-cosmx-fixture.R
# (2 FOVs, a `1e+05` vertex after > 1 MB of integer vertices, an orphan
# polygon, a cell without a polygon).

.cosmx_skip <- function() {
    skip_if_not_installed("arrow")
    skip_if_not_installed("Giotto")
}

.cosmx_obj <- function(d, ...) {
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    suppressMessages(
        importCosMxDisk(d, backend = p, version = "v6")$create_gobject(
            ..., verbose = FALSE)
    )
}

test_that("polygon CSV with coordinates written as 1e+05 is read", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir()
    g <- suppressWarnings(.cosmx_obj(d))

    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    # c_1_2_1 is centered at (5100, 99990); its top edge is the 1e+05 vertex
    expect_equal(sl[cell_ID == "c_1_2_1", c(sdimx, sdimy)], c(5100, 99990))
})

test_that("polygons of cells with no expression are dropped, with a warning", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir()
    expect_warning(g <- .cosmx_obj(d),
                   "1 of 15003 polygons belong to cells")

    poly_ids <- GiottoClass::spatIDs(GiottoClass::getPolygonInfo(g,
        return_giottoPolygon = TRUE))
    expr_ids <- colnames(GiottoClass::getExpression(g, feat_type = "rna",
        values = "raw")[])
    expect_false("c_1_2_4" %in% poly_ids)
    expect_setequal(poly_ids, expr_ids)
    # the filter must not reorder vertices: every fixture cell is a 20 x 20
    # square, which a scrambled ring would not be
    areas <- terra::expanse(GiottoClass::getPolygonInfo(g,
        return_giottoPolygon = TRUE))
    area <- areas[[which(vapply(areas, is.numeric, logical(1L)))[[1L]]]]
    expect_true(all(abs(area - 400) < 1e-6))
})

test_that("cells with expression but no polygon are dropped", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir()
    g <- suppressWarnings(.cosmx_obj(d))

    ids <- GiottoClass::spatIDs(g)
    expect_length(ids, 15002L)
    expect_false("c_1_2_3" %in% ids)
    expect_setequal(GiottoClass::pDataDT(g)$cell_ID, ids)
})

test_that("an export with no orphan polygons loads without a warning", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir()
    pf <- file.path(d, "Fixture-polygons.csv")
    writeLines(grep("c_1_2_4", readLines(pf), value = TRUE, invert = TRUE), pf)
    expect_no_warning(.cosmx_obj(d), message = "polygons belong to cells")
})

test_that("control probes are split out by the default split_keyword", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir()
    g <- suppressWarnings(.cosmx_obj(d))
    expect_setequal(GiottoClass::list_expression(g)$feat_type,
                    c("rna", "negprobes", "falsecode"))
})
