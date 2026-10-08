# VisiumHDDiskReader. Every check runs on a synthetic Space Ranger `outs`
# directory written by helper-visiumhd-fixture.R, and again on a real one
# when GIOTTODISK_VISIUMHD_DATA points at it (10x's "Visium HD Tiny 3'
# Dataset" takes about a minute). Ground truth is the in-memory reader for
# parity, and the shipped 10x matrices, read directly with hdf5r, for
# aggregation.

skip_if_not_installed("Giotto")
skip_if_not_installed("hdf5r")
skip_if_not_installed("tilework")

.vhd_fixture <- .write_visiumhd_fixture()
.vhd_sources <- list(fixture = .vhd_fixture$root,
                     dataset = Sys.getenv("GIOTTODISK_VISIUMHD_DATA", ""))

.vhd_skip <- function(dir) {
    if (!nzchar(dir) || !dir.exists(file.path(dir, "binned_outputs"))) {
        skip("GIOTTODISK_VISIUMHD_DATA not set")
    }
}

# 10x feature-barcode h5 -> genes x barcodes, named as the disk reader names
.vhd_h5 <- function(f) {
    h <- hdf5r::H5File$new(f, mode = "r")
    on.exit(h$close_all())
    Matrix::sparseMatrix(
        i = h[["matrix/indices"]][] + 1L, p = h[["matrix/indptr"]][],
        x = as.numeric(h[["matrix/data"]][]), dims = h[["matrix/shape"]][],
        dimnames = list(.disambiguate_feat_ids(h[["matrix/features/name"]][]),
                        h[["matrix/barcodes"]][])
    )
}

.vhd_dense <- function(m, rows, cols) {
    out <- Matrix::Matrix(0, length(rows), length(cols), sparse = TRUE,
                          dimnames = list(rows, cols))
    r <- intersect(rownames(m), rows)
    out[r, ] <- m[r, cols, drop = FALSE]
    out
}

.q <- function(x) suppressWarnings(suppressMessages(x))

test_that("importVisiumHDDisk requires a backend", {
    expect_error(importVisiumHDDisk(), "backend")
})

test_that("[fixture] disk objects hold the fixture's own counts and hierarchy", {
    fx <- .vhd_fixture
    withr::local_options(list(giottodisk.dgc_max_rows = Inf,
                              giottodisk.dgc_max_cols = Inf))
    g <- .q(importVisiumHDDisk(file.path(fx$root, "binned_outputs"),
        backend = withr::local_tempdir(), bin = 2)$
        create_gobject(load_bin_mapping = TRUE, load_image = FALSE,
                       verbose = FALSE))

    m <- storeRead(GiottoClass::getExpression(g, output = "exprObj")[],
                   output = "dgcmatrix")
    keep <- Matrix::rowSums(fx$m2) > 0          # G5 has no counts
    expect_identical(rownames(m), rownames(fx$m2)[keep])
    expect_equal(as.matrix(m), as.matrix(fx$m2[keep, colnames(m)]))

    cm <- GiottoClass::pDataDT(g, spat_unit = "bin002")
    exp <- fx$grid[match(cm$cell_ID, barcode)]
    expect_identical(cm$bin008, exp$bc8)
    expect_identical(cm$cell, ifelse(is.na(exp$cell), NA_character_,
                                     sprintf("cellid_%09d-1", exp$cell)))
    expect_identical(cm$in_cell, !is.na(exp$cell))
})

for (src in names(.vhd_sources)) local({
    dir <- .vhd_sources[[src]]
    tag <- sprintf("[%s] ", src)

    test_that(paste0(tag, "binned object matches the in-memory reader"), {
        .vhd_skip(dir)
        withr::local_options(list(giottodisk.dgc_max_rows = Inf,
                                  giottodisk.dgc_max_cols = Inf))
        bdir <- file.path(dir, "binned_outputs")
        gd <- .q(importVisiumHDDisk(bdir, backend = withr::local_tempdir(),
            bin = 8)$create_gobject(load_image = FALSE, verbose = FALSE))
        gm <- .q(Giotto::importVisiumHD(bdir, bin = 8)$
            create_gobject(load_image = FALSE, verbose = FALSE))

        ed <- GiottoClass::getExpression(gd, output = "exprObj")[]
        expect_s4_class(ed, "parquetExprStore")
        expect_identical(GiottoClass::spatIDs(gd), GiottoClass::spatIDs(gm))
        expect_identical(GiottoClass::featIDs(gd), GiottoClass::featIDs(gm))
        expect_identical(
            as.data.frame(GiottoClass::getSpatialLocations(gd, output = "data.table")),
            as.data.frame(GiottoClass::getSpatialLocations(gm, output = "data.table")))
        mm <- GiottoClass::getExpression(gm, output = "matrix")
        expect_equal(max(abs(storeRead(ed, output = "dgcmatrix") - mm)), 0)

        # same rows; the disk reader also puts them in expression order
        cd <- GiottoClass::pDataDT(gd)
        expect_identical(cd$cell_ID, GiottoClass::spatIDs(gd))
        expect_equal(cd[order(cell_ID)], GiottoClass::pDataDT(gm)[order(cell_ID)],
                     ignore_attr = TRUE)
    })

    test_that(paste0(tag, "segmented object matches the in-memory reader"), {
        .vhd_skip(dir)
        withr::local_options(list(giottodisk.dgc_max_rows = Inf,
                                  giottodisk.dgc_max_cols = Inf))
        sdir <- file.path(dir, "segmented_outputs")
        gd <- .q(importVisiumHDDisk(sdir, backend = withr::local_tempdir())$
            create_gobject(load_image = FALSE, verbose = FALSE))
        gm <- .q(Giotto::createGiottoVisiumHDObjectCell(sdir,
            load_image = FALSE, verbose = FALSE))

        expect_identical(GiottoClass::spatIDs(gd), GiottoClass::spatIDs(gm))
        md <- storeRead(GiottoClass::getExpression(gd, output = "exprObj")[],
                        output = "dgcmatrix")
        expect_equal(
            max(abs(md - GiottoClass::getExpression(gm, output = "matrix"))), 0)

        for (nm in c("cell", "nucleus")) {
            pd <- GiottoClass::getPolygonInfo(gd, nm, return_giottoPolygon = TRUE)
            pm <- GiottoClass::getPolygonInfo(gm, nm, return_giottoPolygon = TRUE)
            expect_s4_class(pd[], "parquetGeomStore")
            sv <- as.terra(pd[])
            am <- stats::setNames(suppressWarnings(terra::expanse(pm[])),
                                  pm$poly_ID)
            expect_setequal(sv$poly_ID, names(am))
            expect_equal(suppressWarnings(terra::expanse(sv)),
                         unname(am[sv$poly_ID]))
        }
        sd <- GiottoClass::getSpatialLocations(gd, spat_unit = "cell",
                                               output = "data.table")
        sm <- GiottoClass::getSpatialLocations(gm, spat_unit = "cell",
                                               output = "data.table")
        expect_equal(sd[match(sm$cell_ID, cell_ID), .(sdimx, sdimy)],
                     sm[, .(sdimx, sdimy)], ignore_attr = TRUE)
    })

    test_that(paste0(tag, "2 um bin points aggregate into the disk cells as ",
                     "the shipped cell matrix"), {
        .vhd_skip(dir)
        proj <- withr::local_tempdir()
        gp <- .q(importVisiumHDDisk(file.path(dir, "binned_outputs"),
            backend = proj, bin = 2)$load_transcripts(verbose = FALSE))[[1L]]
        expect_s4_class(gp, "giottoBinPoints")
        raw2 <- .vhd_h5(file.path(dir,
            "binned_outputs/square_002um/raw_feature_bc_matrix.h5"))
        expect_equal(sum(gp@counts$x), sum(raw2@x))

        sdir <- file.path(dir, "segmented_outputs")
        gc <- .q(importVisiumHDDisk(sdir, backend = proj)$
            create_gobject(load_expression = FALSE, load_polygons = "cell",
                           load_image = FALSE, verbose = FALSE))
        gpoly <- GiottoClass::getPolygonInfo(gc, "cell",
                                             return_giottoPolygon = TRUE)
        expect_s4_class(gpoly[], "parquetGeomStore")
        ov <- calculateOverlap(gpoly, gp, return_gpolygon = FALSE,
                               verbose = FALSE)
        agg <- overlapToMatrix(ov, feat_count_column = "count")
        shipped <- .vhd_h5(file.path(sdir, "raw_feature_cell_matrix.h5"))

        rows <- union(rownames(shipped), rownames(agg))
        cols <- colnames(shipped)
        expect_setequal(colnames(agg), cols)
        expect_equal(.vhd_dense(agg, rows, cols), .vhd_dense(shipped, rows, cols))
    })

    test_that(paste0(tag, "several units build into one object, linked by ",
                     "parent units"), {
        .vhd_skip(dir)
        proj <- withr::local_tempdir()
        g <- .q(importVisiumHDDisk(file.path(dir, "binned_outputs"),
            backend = proj, bin = 8, tissue_only = TRUE)$
            create_gobject(bin = c(2, 8), load_bin_mapping = TRUE,
                           load_image = FALSE, verbose = FALSE))
        g <- .q(importVisiumHDDisk(file.path(dir, "segmented_outputs"),
            backend = proj)$create_gobject(gobject = g, load_image = FALSE,
                                           verbose = FALSE))

        expect_setequal(names(g@expression), c("bin002", "bin008", "cell"))
        c2 <- GiottoClass::pDataDT(g, spat_unit = "bin002")
        c8 <- GiottoClass::pDataDT(g, spat_unit = "bin008")
        expect_true(all(c("bin008", "cell", "in_cell") %in% names(c2)))
        expect_false("cell" %in% names(c8))   # an 8 um bin spans two cells
        expect_true(all(stats::na.omit(c2$cell) %in%
                        GiottoClass::spatIDs(g, "cell")))
        expect_true(all(c2$bin008 %in% GiottoClass::spatIDs(g, "bin008")))
    })

    test_that(paste0(tag, "loading the 2 um bin points keeps the bin's cell ",
                     "metadata"), {
        # the parent reader replaces its barcodes with the 2 um ones after
        # this step, which leaves an 8 um object with no cell metadata rows
        .vhd_skip(dir)
        g <- .q(importVisiumHDDisk(file.path(dir, "binned_outputs"),
            backend = withr::local_tempdir(), bin = 8)$
            create_gobject(load_transcripts = TRUE, load_image = FALSE,
                           verbose = FALSE))
        expect_s4_class(GiottoClass::getFeatureInfo(g, return_giottoPoints = TRUE),
                        "giottoBinPoints")
        expect_identical(GiottoClass::pDataDT(g, spat_unit = "bin008")$cell_ID,
                         GiottoClass::spatIDs(g, "bin008"))
    })

    test_that(paste0(tag, "the Giotto wrappers build the disk object with ",
                     "`backend =`"), {
        .vhd_skip(dir)
        if (!"backend" %in% names(formals(Giotto::createGiottoVisiumHDObjectCell))) {
            skip("installed Giotto predates `backend` on the VisiumHD wrappers")
        }
        gb <- .q(Giotto::createGiottoVisiumHDObjectBin(
            file.path(dir, "binned_outputs"), bin = 8, load_image = FALSE,
            verbose = FALSE, backend = withr::local_tempdir()))
        expect_false(is.null(gb@source))
        expect_s4_class(GiottoClass::getExpression(gb, output = "exprObj")[],
                        "parquetExprStore")

        gc <- .q(Giotto::createGiottoVisiumHDObjectCell(
            file.path(dir, "segmented_outputs"), load_image = FALSE,
            load_transcripts = TRUE, verbose = FALSE,
            backend = withr::local_tempdir()))
        expect_s4_class(GiottoClass::getPolygonInfo(gc, "cell",
            return_giottoPolygon = TRUE)[], "parquetGeomStore")
        expect_s4_class(GiottoClass::getFeatureInfo(gc, return_giottoPoints = TRUE),
                        "giottoBinPoints")
    })
})
