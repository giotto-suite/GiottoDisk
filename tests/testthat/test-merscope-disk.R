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

# A mosaic like MERSCOPE's: too large to read whole, uncompressed, in strips,
# with no geotransform (GDAL's baseline profile, its .aux.xml sidecar
# removed); every row and column different, so a flip shows. Transform: 2
# pixels per micron, offset 4 / 6 pixels.
.mk_mosaic <- function(d, ncol = 4200L) {
    dir.create(file.path(d, "images"), showWarnings = FALSE)
    v <- outer(seq_len(64L) * 50L, seq_len(ncol) %% 37L, `+`)
    r <- terra::rast(nrows = 64, ncols = ncol, xmin = 0, xmax = ncol,
        ymin = 0, ymax = 64, vals = as.vector(t(v)))
    f <- file.path(d, "images", "mosaic_DAPI_z0.tif")
    terra::writeRaster(r, f, datatype = "INT2U", overwrite = TRUE,
        gdal = c("PROFILE=BASELINE", "COMPRESS=NONE"))
    unlink(paste0(f, ".aux.xml"))
    utils::write.table(matrix(c(2, 0, 4, 0, 2, 6, 0, 0, 1), 3, byrow = TRUE),
        file.path(d, "images", "micron_to_mosaic_pixel_transform.csv"),
        sep = " ", row.names = FALSE, col.names = FALSE)
    f
}

test_that("only large images without tiles or overviews need converting", {
    .mer_skip()
    d <- withr::local_tempdir()
    f <- .mk_mosaic(d)
    expect_false(.tif_tile_readable(f))
    small <- file.path(d, "small.tif")
    terra::writeRaster(terra::rast(nrows = 40, ncols = 80, vals = 1), small)
    expect_true(.tif_tile_readable(small))
    cog <- file.path(d, "cog.tif")
    terra::writeRaster(suppressWarnings(terra::rast(f)), cog, filetype = "COG",
        datatype = "INT2U", gdal = "BLOCKSIZE=512")
    expect_true(.tif_tile_readable(cog))
    expect_true(.tif_tile_readable(file.path(d, "not-an-image.tif")))
})

test_that("load_image() converts a mosaic once to a COG in the vault", {
    .mer_skip()
    d <- .mk_merscope_disk_dir()
    f <- .mk_mosaic(d)
    proj <- .proj()
    x <- importMerscopeDisk(d, backend = proj)
    img <- x$load_image(path = f)
    cog <- img@file_path
    expect_true(startsWith(cog, normalizePath(.gdsrc_vault_dir(proj))))
    expect_true(.tif_tile_readable(cog))
    expect_match(paste(terra::describe(cog), collapse = "\n"), "Overviews:")

    # Placed as the mosaic itself: the same value at every point Giotto reads
    plain <- x$load_image(path = f, convert = FALSE)
    expect_identical(plain@file_path, f)
    expect_equal(as.vector(terra::ext(img)), as.vector(terra::ext(plain)))
    e <- as.vector(terra::ext(plain))
    xy <- as.matrix(expand.grid(seq(e[1] + 1, e[2] - 1, length.out = 40),
                                seq(e[3] + 0.1, e[4] - 0.1, length.out = 16)))
    expect_identical(terra::extract(img@raster_object, xy)[, 1],
                     terra::extract(plain@raster_object, xy)[, 1])
    # and in the source's row order, which tile readers such as GiottoLens
    # use: terra shows the source (no geotransform) upside down, the COG not
    a <- suppressWarnings(terra::rast(f))
    b <- terra::rast(cog)
    expect_identical(terra::values(b, row = 1, nrows = 1),
                     terra::values(a, row = nrow(a), nrows = 1))

    # Reused while the mosaic is unchanged; converted again once it changes
    mt <- file.mtime(cog)
    again <- x$load_image(path = f)
    expect_identical(again@file_path, cog)
    expect_identical(file.mtime(cog), mt)
    Sys.setFileTime(f, Sys.time() + 60)
    expect_false(identical(x$load_image(path = f)@file_path, cog))
})

test_that("a georeferenced source keeps its orientation in the COG", {
    .mer_skip()
    d <- withr::local_tempdir()
    v <- outer(seq_len(64L) * 50L, seq_len(4200L) %% 37L, `+`)
    r <- terra::rast(nrows = 64, ncols = 4200, xmin = 0, xmax = 4200,
        ymin = 0, ymax = 64, vals = as.vector(t(v)))
    f <- file.path(d, "geo.tif")
    terra::writeRaster(r, f, datatype = "INT2U")
    expect_false(.tif_tile_readable(f))
    cog <- .image_cog(f, gDirSource(path = .proj()))
    a <- terra::rast(f)
    b <- terra::rast(cog)
    expect_identical(terra::datatype(b), "INT2U")
    expect_equal(as.vector(terra::ext(b)), as.vector(terra::ext(a)))
    expect_identical(terra::values(b), terra::values(a))
})

test_that("a snapshot keeps the converted image; convert_images = FALSE skips it", {
    .mer_skip()
    d <- .mk_merscope_disk_dir()
    f <- .mk_mosaic(d)
    x <- importMerscopeDisk(d, backend = .proj())
    g <- suppressMessages(suppressWarnings(
        x$create_gobject(load_images = "DAPI_z0", verbose = FALSE)))
    cog <- g@images[[1L]]@file_path
    expect_false(identical(normalizePath(cog), normalizePath(f)))
    placed <- as.vector(terra::ext(g@images[[1L]]@raster_object))
    suppressWarnings(snapshotSave(g@source, g, name = "with_image", verbose = FALSE))
    sourcePrune(g@source)
    expect_true(file.exists(cog))
    # Loading rebuilds the raster from the file: it keeps its placement
    loaded <- suppressMessages(suppressWarnings(snapshotLoad(g@source)))
    img <- GiottoClass::reconnect(loaded@images[[1L]])
    expect_equal(as.vector(terra::ext(img@raster_object)), placed)

    g2 <- suppressMessages(suppressWarnings(
        importMerscopeDisk(d, backend = .proj())$create_gobject(
            load_images = "DAPI_z0", convert_images = FALSE, verbose = FALSE)))
    expect_identical(normalizePath(g2@images[[1L]]@file_path), normalizePath(f))
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
