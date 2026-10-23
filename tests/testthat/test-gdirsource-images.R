# Adopting images into a gDirSource vault. GiottoClass reads a JPEG-2000 TIFF
# through a VRT under tempdir() (GiottoClass adr/0005), so terra::sources()
# names the VRT, not the image: the vault must copy @file_path instead.

# a JPEG-2000 giottoLargeImage over a private copy of the GiottoClass fixture,
# so a test can delete the original
.jp2k_image <- function(name = "jp2k") {
    fx <- system.file("extdata", "toy_jp2k.ome.tif", package = "GiottoClass")
    if (!nzchar(fx)) testthat::skip("GiottoClass JPEG-2000 fixture not found")
    d <- tempfile("jp2k_")
    dir.create(d)
    f <- file.path(d, basename(fx))
    file.copy(fx, f)
    img <- GiottoClass::createGiottoLargeImage(f, name = name, verbose = FALSE)
    if (!grepl("\\.vrt$", terra::sources(img@raster_object))) {
        testthat::skip("GiottoClass does not read JPEG-2000 through a VRT")
    }
    img
}

.img_src <- function() {
    td <- tempfile("imgsrc_")
    list(td = td, src = gDirSource(td))
}

.img_uid <- function(img) basename(dirname(img@file_path))

describe("sourceAdopt(gDirSource, giottoLargeImage)", {
    it("copies the JPEG-2000 file, not the VRT it is read through", {
        x <- .img_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        img <- .jp2k_image()
        orig <- img@file_path

        out <- sourceAdopt(x$src, img)
        expect_true(startsWith(normalizePath(out@file_path),
            normalizePath(.gdsrc_vault_dir(x$td))))
        expect_identical(unname(tools::md5sum(out@file_path)),
            unname(tools::md5sum(orig)))
        expect_identical(unname(terra::values(out@raster_object)),
            unname(terra::values(img@raster_object)))
    })

    it("leaves an image that reads after the original is gone", {
        x <- .img_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        img <- .jp2k_image()
        vals <- unname(terra::values(img@raster_object))

        out <- sourceAdopt(x$src, img)
        unlink(dirname(img@file_path), recursive = TRUE)
        out <- GiottoClass::reconnect(out)
        expect_identical(unname(terra::values(out@raster_object)), vals)
    })

    it("is seen by sourceContains, so a second adopt copies nothing", {
        x <- .img_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        img <- .jp2k_image()
        expect_false(sourceContains(x$src, img))

        out <- sourceAdopt(x$src, img)
        expect_true(sourceContains(x$src, out))
        n_vault <- length(list.dirs(.gdsrc_vault_dir(x$td),
            recursive = FALSE))
        again <- sourceAdopt(x$src, out)
        expect_identical(again@file_path, out@file_path)
        expect_length(list.dirs(.gdsrc_vault_dir(x$td), recursive = FALSE),
            n_vault)
    })

    it("records a content fingerprint, not a per-object hash", {
        x <- .img_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        a <- sourceAdopt(x$src, .jp2k_image("a"))
        tif <- system.file("extdata", "toy_intensity.tif",
            package = "GiottoClass")
        skip_if_not(nzchar(tif))
        b <- sourceAdopt(x$src,
            GiottoClass::createGiottoLargeImage(tif, verbose = FALSE))

        ha <- x$src[.img_uid(a), "hash"]
        expect_identical(ha, GiottoClass::fingerprint_data(a@raster_object))
        expect_false(identical(ha, x$src[.img_uid(b), "hash"]))
    })
})

describe("sourceAdopt(gDirSource, SpatRaster)", {
    it("refuses a raster read through a VRT", {
        x <- .img_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        img <- .jp2k_image()
        expect_error(sourceAdopt(x$src, img@raster_object), "VRT")
        expect_length(list.dirs(.gdsrc_vault_dir(x$td), recursive = FALSE), 0L)
    })
})

describe("snapshotSave with a JPEG-2000 image", {
    skip_if_not_installed("GiottoData")
    g0 <- GiottoUtils::gwith_options(list(giotto.use_conda = FALSE),
        GiottoData::loadGiottoMini("visium", verbose = FALSE)
    )

    it("copies the image into the vault when exporting", {
        x <- .img_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        img <- .jp2k_image()
        g <- GiottoClass::setGiotto(g, img, verbose = FALSE)

        g <- snapshotSave(x$src, g, name = "v1", verbose = FALSE)
        out <- GiottoClass::getGiottoImage(g, name = "jp2k")
        expect_identical(unname(tools::md5sum(out@file_path)),
            unname(tools::md5sum(img@file_path)))
        expect_true(.img_uid(out) %in% .ss_gdsrc_detect_uid_images(g))
        expect_true(.snap_tag(x$td, "v1") %in%
            unlist(x$src[.img_uid(out), "giottosave"]))
    })

    it("keeps a reference to the original when not exporting", {
        x <- .img_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        img <- .jp2k_image()
        g <- GiottoClass::setGiotto(g, img, verbose = FALSE)

        g <- snapshotSave(x$src, g, name = "v1", export_image = FALSE,
            verbose = FALSE)
        out <- GiottoClass::getGiottoImage(g, name = "jp2k")
        expect_identical(out@file_path, img@file_path)
        expect_false(sourceContains(x$src, out))
    })
})
