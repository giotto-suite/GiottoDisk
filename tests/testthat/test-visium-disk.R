# VisiumDiskReader. Every check runs on a synthetic Space Ranger `outs`
# directory written by helper-visium-fixture.R, and again on a real one when
# GIOTTODISK_VISIUM_DATA points at it. Ground truth is the in-memory reader.

skip_if_not_installed("Giotto")
skip_if_not_installed("hdf5r")
skip_if_not("VisiumReader" %in% methods::getClasses(asNamespace("Giotto")))

.vis_fixture <- .write_visium_fixture()
.vis_sources <- list(fixture = .vis_fixture$root,
                     dataset = Sys.getenv("GIOTTODISK_VISIUM_DATA", ""))

.vis_skip <- function(dir) {
    if (!nzchar(dir) || !dir.exists(file.path(dir, "spatial"))) {
        skip("GIOTTODISK_VISIUM_DATA not set")
    }
}

.q <- function(x) suppressWarnings(suppressMessages(x))

test_that("importVisiumDisk requires a backend", {
    expect_error(importVisiumDisk(), "backend")
})

test_that("the disk reader extends the in-memory one", {
    x <- importVisiumDisk(.vis_fixture$root, backend = withr::local_tempdir())
    expect_s4_class(x, "VisiumDiskReader")
    expect_true(methods::extends("VisiumDiskReader", "VisiumReader"))
    expect_identical(names(formals(Giotto::importVisium(.vis_fixture$root)$
        create_gobject)), setdiff(names(formals(x$create_gobject)), "gobject"))
})

for (src in names(.vis_sources)) local({
    dir <- .vis_sources[[src]]
    tag <- sprintf("[%s] ", src)

    test_that(paste0(tag, "disk object matches the in-memory reader"), {
        .vis_skip(dir)
        withr::local_options(list(giottodisk.dgc_max_rows = Inf,
                                  giottodisk.dgc_max_cols = Inf))
        gd <- .q(importVisiumDisk(dir, backend = withr::local_tempdir())$
            create_gobject(load_image = FALSE, verbose = FALSE))
        gm <- .q(Giotto::importVisium(dir)$
            create_gobject(load_image = FALSE, verbose = FALSE))

        expect_identical(GiottoClass::list_expression(gd)$feat_type,
                         GiottoClass::list_expression(gm)$feat_type)
        for (ft in GiottoClass::list_expression(gm)$feat_type) {
            ed <- GiottoClass::getExpression(gd, feat_type = ft,
                                             output = "exprObj")[]
            expect_s4_class(ed, "parquetExprStore")
            expect_identical(GiottoClass::featIDs(gd, feat_type = ft),
                             GiottoClass::featIDs(gm, feat_type = ft))
            mm <- GiottoClass::getExpression(gm, feat_type = ft,
                                             output = "matrix")
            expect_equal(max(abs(storeRead(ed, output = "dgcmatrix") - mm)), 0)
            expect_equal(GiottoClass::pDataDT(gd, feat_type = ft),
                         GiottoClass::pDataDT(gm, feat_type = ft),
                         ignore_attr = TRUE)
        }
        expect_identical(GiottoClass::spatIDs(gd), GiottoClass::spatIDs(gm))
        expect_identical(
            as.data.frame(GiottoClass::getSpatialLocations(gd, output = "data.table")),
            as.data.frame(GiottoClass::getSpatialLocations(gm, output = "data.table")))
        # disk polygons are written to the vault by the backed setter
        pd <- GiottoClass::getPolygonInfo(gd, name = "cell",
                                          return_giottoPolygon = FALSE)
        expect_s4_class(pd, "parquetGeomStore")
        pd <- storeRead(pd, output = "terra")
        pm <- GiottoClass::getPolygonInfo(gm, name = "cell",
                                          return_giottoPolygon = FALSE)
        expect_identical(pd$poly_ID, pm$poly_ID)
        expect_equal(terra::geom(pd), terra::geom(pm))
        expect_identical(GiottoClass::instructions(gd, "micron_scale"),
                         GiottoClass::instructions(gm, "micron_scale"))
    })
})

test_that("[fixture] createGiottoVisiumObject(backend =) builds a backed object", {
    proj <- withr::local_tempdir()
    g <- .q(Giotto::createGiottoVisiumObject(.vis_fixture$root,
        gene_column_index = 2, backend = proj, verbose = FALSE))
    expect_false(is.null(g@source))
    expect_s4_class(GiottoClass::getExpression(g, output = "exprObj")[],
                    "parquetExprStore")
    expect_identical(GiottoClass::spatIDs(g), .vis_fixture$grid$barcode)
    expect_identical(GiottoClass::list_images_names(g), "image")
})

test_that("[fixture] a Visium disk object survives a snapshot and reload", {
    # loadGiotto() otherwise looks for a python environment to reconnect
    withr::local_options(list(giotto.use_conda = FALSE,
                              tilework.warn_sequential = FALSE))
    proj <- withr::local_tempdir()
    g <- .q(importVisiumDisk(.vis_fixture$root, backend = proj)$
        create_gobject(load_image = FALSE, verbose = FALSE))
    g <- .q(GiottoClass::saveGiotto(g, name = "visium"))
    g2 <- .q(GiottoClass::loadGiotto(proj, name = "visium"))
    expect_identical(GiottoClass::spatIDs(g2), GiottoClass::spatIDs(g))
    expect_identical(GiottoClass::list_spatial_info_names(g2), "cell")
    expect_equal(
        storeRead(GiottoClass::getExpression(g2, output = "exprObj")[],
                  output = "dgcmatrix"),
        storeRead(GiottoClass::getExpression(g, output = "exprObj")[],
                  output = "dgcmatrix"))
})
