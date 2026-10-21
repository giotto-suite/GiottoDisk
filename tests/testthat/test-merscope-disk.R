# End-to-end MERSCOPE disk reader tests. Fixture: helper-merscope-fixture.R
# (6 cells over 2 FOVs, 19-digit cell IDs, csv or parquet transcripts, parquet
# or HDF5 boundaries).

.mer_skip <- function() {
    skip_if_not_installed("arrow")
    skip_if_not_installed("wk")
}

.proj <- function() {
    p <- tempfile("merscope_proj_")
    dir.create(p)
    p
}

.disk_obj <- function(d, ..., polygon_format = "parquet") {
    suppressMessages(suppressWarnings(
        importMerscopeDisk(d, backend = .proj(),
            polygon_format = polygon_format)$create_gobject(..., verbose = FALSE)
    ))
}

.expr_cols <- function(g, feat_type = "rna") {
    withr::local_options(list(giottodisk.dgc_max_rows = Inf,
        giottodisk.dgc_max_cols = Inf))
    ed <- GiottoClass::getExpression(g, feat_type = feat_type,
        output = "exprObj")[]
    storeRead(ed, output = "dgcmatrix")
}

.poly_ids <- function(g) {
    GiottoClass::spatIDs(GiottoClass::getPolygonInfo(g,
        return_giottoPolygon = TRUE))
}


test_that("disk object matches the in-memory reader (csv and parquet tx)", {
    .mer_skip()
    for (fmt in c("csv", "parquet")) {
        d <- .mk_merscope_disk_dir(tx_format = fmt)
        gd <- .disk_obj(d)
        gm <- suppressWarnings(Giotto::importMerscope(d)$create_gobject(
            verbose = FALSE))

        for (ft in c("rna", "neg_probe")) {
            md <- .expr_cols(gd, ft)
            mm <- GiottoClass::getExpression(gm, feat_type = ft,
                output = "matrix")
            expect_setequal(colnames(md), colnames(mm))
            expect_equal(as.matrix(md[sort(rownames(md)), sort(colnames(md))]),
                as.matrix(mm[sort(rownames(mm)), sort(colnames(mm))]))

            pd <- data.table::as.data.table(GiottoClass::getFeatureInfo(gd,
                feat_type = ft, return_giottoPoints = TRUE)[], geom = "XY")
            pm <- terra::crds(GiottoClass::getFeatureInfo(gm, feat_type = ft,
                return_giottoPoints = TRUE)[])
            expect_equal(sort(pd$x), sort(pm[, 1]))
            expect_equal(sort(pd$y), sort(pm[, 2]))
        }
        expect_setequal(.poly_ids(gd), .poly_ids(gm))
        sd <- GiottoClass::getSpatialLocations(gd, output = "data.table")
        sm <- GiottoClass::getSpatialLocations(gm, output = "data.table")
        data.table::setorder(sd, cell_ID)
        data.table::setorder(sm, cell_ID)
        expect_identical(sd$cell_ID, sm$cell_ID)
        expect_equal(sd$sdimy, sm$sdimy)
    }
})

test_that("19-digit cell IDs are exact in every data type", {
    .mer_skip()
    g <- .disk_obj(.mk_merscope_disk_dir(), load_cellmeta = TRUE)
    want <- .mer_ids()
    expect_setequal(.poly_ids(g), want)
    expect_setequal(colnames(.expr_cols(g)), want)
    expect_setequal(GiottoClass::pDataDT(g)$cell_ID, want)
})

test_that("FOVs subsets every data type, the expression matrix included", {
    .mer_skip()
    g <- .disk_obj(.mk_merscope_disk_dir(), FOVs = 1, load_cellmeta = TRUE)
    want <- .mer_ids()[4:6]
    expect_setequal(colnames(.expr_cols(g)), want)
    expect_setequal(colnames(.expr_cols(g, "neg_probe")), want)
    expect_setequal(.poly_ids(g), want)
    expect_setequal(GiottoClass::pDataDT(g)$cell_ID, want)
    expect_setequal(GiottoClass::getSpatialLocations(g,
        output = "data.table")$cell_ID, want)
})

test_that("FOVs with csv transcripts keeps only those FOVs' points", {
    .mer_skip()
    g <- .disk_obj(.mk_merscope_disk_dir(tx_format = "csv"), FOVs = 0)
    pd <- data.table::as.data.table(GiottoClass::getFeatureInfo(g,
        feat_type = "rna", return_giottoPoints = TRUE)[])
    expect_setequal(unique(pd$cell_id), .mer_ids()[1:3])
})

test_that("replicated z-planes load one plane", {
    .mer_skip()
    g <- .disk_obj(.mk_merscope_disk_dir(z_planes = 3L))
    expect_length(.poly_ids(g), 6L)
})

test_that("HDF5 boundaries are read on disk", {
    .mer_skip()
    skip_if_not_installed("rhdf5")
    withr::local_options(list(giotto.warn_sequential = FALSE))
    d <- .mk_merscope_disk_dir(bounds = "hdf5")
    g <- .disk_obj(d, polygon_format = "hdf5")
    expect_setequal(.poly_ids(g), .mer_ids())
    gp <- .disk_obj(.mk_merscope_disk_dir())
    sd <- GiottoClass::getSpatialLocations(g, output = "data.table")
    sp <- GiottoClass::getSpatialLocations(gp, output = "data.table")
    data.table::setorder(sd, cell_ID)
    data.table::setorder(sp, cell_ID)
    # same squares as the parquet boundaries, so the same centroids
    expect_equal(sd$sdimx, sp$sdimx)
    expect_equal(sd$sdimy, sp$sdimy)
})

test_that("load_image() works without an explicit transform_path", {
    .mer_skip()
    d <- .mk_merscope_disk_dir()
    dir.create(file.path(d, "images"))
    r <- terra::rast(nrows = 40, ncols = 80, xmin = 0, xmax = 80,
        ymin = 0, ymax = 40, vals = seq_len(3200))
    terra::writeRaster(r, file.path(d, "images", "mosaic_DAPI_z0.tif"))
    # 2 pixels per micron, offset 4 / 6 pixels
    utils::write.table(matrix(c(2, 0, 4, 0, 2, 6, 0, 0, 1), 3, byrow = TRUE),
        file.path(d, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep = " ", row.names = FALSE, col.names = FALSE)
    x <- importMerscopeDisk(d, backend = .proj())
    img <- x$load_image(path = file.path(d, "images", "mosaic_DAPI_z0.tif"))
    e <- as.vector(terra::ext(img))
    expect_equal(unname(e), c(-2, 38, -17, 3))
})

test_that("csv transcripts are converted to parquet once and cached", {
    .mer_skip()
    withr::local_options(list(giottodisk.artifact_dump = withr::local_tempdir()))
    f <- file.path(.mk_merscope_disk_dir(tx_format = "csv"),
        "detected_transcripts.csv")
    a <- .merscope_csv_to_parquet(f)
    mt <- file.mtime(list.files(a, full.names = TRUE))
    b <- .merscope_csv_to_parquet(f)
    expect_identical(a, b)
    expect_identical(file.mtime(list.files(b, full.names = TRUE)), mt)
    tx <- dplyr::collect(arrow::open_dataset(a))
    expect_type(tx$cell_id, "character")
    expect_false(any(grepl("unnamed", names(tx))))
})
