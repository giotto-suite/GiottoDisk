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

# transcripts ####

.tx_count <- function(g, feat_type) {
    gp <- GiottoClass::getFeatureInfo(g, feat_type = feat_type,
        return_giottoPoints = TRUE)
    nrow(storeRead(gp@spatVector, output = "tibble"))
}

test_that("load_transcripts writes every detection, split by feat_type", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    expected <- .mk_cosmx_tx(d)
    withr::local_options(giottodisk.artifact_dump = tempfile("dump_"))
    g <- suppressWarnings(.cosmx_obj(d, load_transcripts = TRUE))

    expect_equal(.tx_count(g, "rna"),
                 sum(expected[c("GENEA", "GENEB")]))
    expect_equal(.tx_count(g, "negprobes"), expected[["Negative1"]])
    expect_equal(.tx_count(g, "falsecode"), expected[["SystemControl1"]])
})

test_that("transcript coordinates and vendor cell assignment are kept", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    .mk_cosmx_tx(d, gz = TRUE)
    withr::local_options(giottodisk.artifact_dump = tempfile("dump_"))
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    r <- suppressWarnings(importCosMxDisk(d, backend = p, version = "v6"))
    st <- r$load_transcripts(output = "store", verbose = FALSE)
    expect_s4_class(st, "parquetGeomTileStore")

    tx <- storeRead(st, output = "tibble")
    hit <- tx[tx$cell == "c_1_2_1", ]
    # the 1e+05 detection, untouched; CellComp survives the blank lead-in
    # the tile store keeps the coordinates as x_index / y_index
    expect_setequal(hit$y_index, c(99990, 100000))
    expect_true("Membrane" %in% tx$CellComp)
})

test_that("the tx_file is converted to parquet once and reused", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    .mk_cosmx_tx(d)
    dump <- tempfile("dump_")
    withr::local_options(giottodisk.artifact_dump = dump)
    tx <- file.path(d, "Fixture_tx_file.csv")

    first <- GiottoDisk:::.cosmx_ensure_parquet(tx)
    stamp <- file.mtime(list.files(first, full.names = TRUE, recursive = TRUE))
    expect_message(second <- GiottoDisk:::.cosmx_ensure_parquet(tx,
        verbose = TRUE), "using cached conversion")
    expect_identical(first, second)
    expect_equal(file.mtime(list.files(second, full.names = TRUE,
        recursive = TRUE)), stamp)
    expect_true(startsWith(first, dump))
})

test_that("load_transcripts restricts to the reader's fovs", {
    .cosmx_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    .mk_cosmx_tx(d)
    withr::local_options(giottodisk.artifact_dump = tempfile("dump_"))
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    r <- suppressWarnings(importCosMxDisk(d, backend = p, version = "v6",
        fovs = 2))
    tx <- storeRead(r$load_transcripts(output = "store", verbose = FALSE),
        output = "tibble")
    expect_setequal(unique(tx$fov), 2L)
    expect_equal(nrow(tx), 3L)
})
