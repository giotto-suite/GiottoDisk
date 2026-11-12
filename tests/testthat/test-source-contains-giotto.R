# sourceContains(gDirSource, giotto): whether a snapshotSave() would move any
# of the object's stores into the vault.

.contains_polys <- function() {
    sq <- function(cx, cy, r = 3) sprintf("POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))",
        cx - r, cy - r, cx + r, cy - r, cx + r, cy + r, cx - r, cy + r, cx - r, cy - r)
    cen <- expand.grid(x = seq(5, 95, 10), y = seq(-95, -5, 10))
    pv <- terra::vect(sq(cen$x, cen$y))
    pv$poly_ID <- sprintf("cell_%03d", seq_len(nrow(cen)))
    GiottoClass::createGiottoPolygon(pv, name = "cell", verbose = FALSE)
}

.contains_points <- function() {
    GiottoUtils::local_seed(1)
    GiottoClass::createGiottoPoints(terra::vect(data.frame(x = stats::runif(50, 0, 100),
        y = stats::runif(50, -100, 0), feat_ID = sample(c("A", "B"), 50, TRUE)),
        geom = c("x", "y")), verbose = FALSE)
}

# a giotto object in project `src` holding `objs`
.contains_object <- function(src, objs) {
    g <- suppressWarnings(suppressMessages(GiottoClass::giotto()))
    for (o in objs) g <- GiottoClass::setGiotto(g, o, verbose = FALSE)
    g@source <- src
    g
}

describe("sourceContains(gDirSource, giotto)", {
    skip_if_not_installed("GiottoClass")

    it("names the stores a save would move, and is TRUE after the save", {
        td <- tempfile("contains_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        src <- suppressMessages(gDirSource(td))
        gp <- .contains_polys()
        gp@spatVector <- suppressMessages(storeWrite(parquetGeomStore(), gp[]))
        pts <- .contains_points()
        pts@spatVector <- suppressMessages(storeWrite(parquetGeomStore(), pts[]))
        g <- .contains_object(src, list(gp, pts))
        res <- sourceContains(src, g)
        expect_false(res)
        expect_setequal(attr(res, "external"), c("polygons 'cell'", "points 'rna'"))
        g <- suppressMessages(snapshotSave(src, g, name = "s1", verbose = FALSE))
        expect_identical(sourceContains(src, g), TRUE)
    })

    it("does not count data held in memory", {
        td <- tempfile("contains_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        src <- suppressMessages(gDirSource(td))
        expect_identical(sourceContains(src, .contains_object(src, list(.contains_polys()))), TRUE)
    })
})
