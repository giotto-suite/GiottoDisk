# A SpatVector written straight to a parquetGeomTileStore, and through
# sourceWrite() with the store type option, as GiottoLens needs for an object
# whose points and polygons are in memory.

.tsv_polys <- function() {
    sq <- function(cx, cy, r = 3) sprintf("POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))",
        cx - r, cy - r, cx + r, cy - r, cx + r, cy + r, cx - r, cy + r, cx - r, cy - r)
    cen <- expand.grid(x = seq(5, 95, 10), y = seq(-95, -5, 10))
    pv <- terra::vect(sq(cen$x, cen$y))
    pv$poly_ID <- sprintf("cell_%03d", seq_len(nrow(cen)))
    pv
}

.tsv_points <- function() {
    GiottoUtils::local_seed(1)
    terra::vect(data.frame(x = stats::runif(400, 0, 100), y = stats::runif(400, -100, 0),
        feat_ID = sample(c("A", "B", "C"), 400, TRUE)), geom = c("x", "y"))
}

.tsv_quiet <- function(expr) GiottoUtils::gwith_options(
    list(tilework.warn_sequential = FALSE), suppressMessages(expr))

describe("storeWrite(parquetGeomTileStore, SpatVector)", {
    it("tiles polygons and points, and leaves no temporary store", {
        dump <- getArtifactDumpDir()
        before <- list.files(dump)
        polys <- .tsv_quiet(storeWrite(parquetGeomTileStore(), .tsv_polys(), threshold = 30L))
        expect_s4_class(polys, "parquetGeomTileStore")
        expect_gt(nrow(polys@tiles@metadata), 1L)
        back <- storeRead(polys, output = "terra")
        expect_setequal(back$poly_ID, .tsv_polys()$poly_ID)
        pts <- .tsv_quiet(storeWrite(parquetGeomTileStore(), .tsv_points(), threshold = 100L))
        expect_gt(nrow(pts@tiles@metadata), 1L)
        expect_equal(nrow(storeRead(pts, output = "terra")), 400)
        expect_setequal(setdiff(list.files(dump), before),
            basename(c(polys@path, pts@path)))
    })

    it("writes nothing for an empty SpatVector", {
        expect_null(storeWrite(parquetGeomTileStore(), .tsv_polys()[0]))
    })
})

describe("sourceWrite(gDirSource, giotto) with tiled geometries", {
    it("writes the polygons and points as tile stores in the vault", {
        skip_if_not_installed("GiottoClass")
        td <- tempfile("tilesv_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        src <- suppressMessages(gDirSource(td))
        g <- suppressWarnings(suppressMessages(GiottoClass::giotto()))
        g <- GiottoClass::setGiotto(g, GiottoClass::createGiottoPolygon(.tsv_polys(),
            name = "cell", verbose = FALSE), verbose = FALSE)
        g <- GiottoClass::setGiotto(g, GiottoClass::createGiottoPoints(.tsv_points(),
            verbose = FALSE), verbose = FALSE)
        g2 <- GiottoUtils::gwith_options(
            list(giotto.gdsrc_spatvector_format = "parquetGeomTile"),
            .tsv_quiet(sourceWrite(src, g)))
        expect_s4_class(g2@spatial_info$cell@spatVector, "parquetGeomTileStore")
        expect_s4_class(g2@feat_info$rna@spatVector, "parquetGeomTileStore")
        reg <- as.data.frame(src)
        expect_setequal(reg$uid, c(g2@spatial_info$cell@spatVector@uid,
            g2@feat_info$rna@spatVector@uid))
        expect_identical(unique(reg$store), "parquetGeomTile")
        expect_setequal(list.files(file.path(td, "artifacts")), reg$uid)
    })
})
