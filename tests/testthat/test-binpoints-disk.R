# Store-backed giottoBinPoints: 2 um bins on disk, aggregated into polygons. Ground truth is
# the fixture's 2 um matrix summed over the bins inside each polygon.

skip_if_not_installed("Giotto")
skip_if_not_installed("hdf5r")
skip_if_not_installed("tilework")
skip_if_not_installed("duckdb")

.q <- function(x) suppressWarnings(suppressMessages(x))

# rectangle polygons with their edges between bin centres
.bp_rects <- function(rects) {
    wkt <- sprintf("POLYGON ((%1$s %3$s, %2$s %3$s, %2$s %4$s, %1$s %4$s, %1$s %3$s))",
        rects$xmin, rects$xmax, rects$ymin, rects$ymax)
    sv <- terra::vect(wkt)
    sv$poly_ID <- rects$poly_ID
    sv
}

test_that("disk bin points sum the bins inside each polygon", {
    fx <- .write_visiumhd_fixture()
    withr::local_options(list(giottodisk.dgc_max_rows = Inf,
                              giottodisk.dgc_max_cols = Inf,
                              giottodisk.spatial_query_engine = "duckdb"))
    proj <- withr::local_tempdir()
    g <- .q(importVisiumHDDisk(file.path(fx$root, "binned_outputs"),
        backend = proj, bin = 2)$create_gobject(load_image = FALSE,
                                                verbose = FALSE))
    pe <- GiottoClass::getExpression(g, spat_unit = "bin002",
                                     output = "exprObj")[]
    sl <- GiottoClass::getSpatialLocations(g, spat_unit = "bin002",
                                           output = "data.table")
    bp <- .q(.binpoints_disk(pe, sl, gDirSource(path = proj)))
    expect_s4_class(bp, "giottoBinPoints")
    expect_true(.gbp_store_backed(bp))
    expect_identical(featIDs(bp), pe@feat_ids)
    # bins with counts: a scan of the store agrees with its cached marginals
    rows <- seq_len(pe@n_cells)
    no_stats <- pe
    no_stats@stats <- list()
    expect_identical(.pe_rows_with_counts(no_stats, rows),
                     .pe_rows_with_counts(pe, rows))
    # only bins with counts
    expect_equal(nrow(bp), sum(Matrix::colSums(fx$m2[, sl$cell_ID]) > 0))

    # two overlapping rectangles and one apart; edges on even coordinates
    r <- lapply(sl[, c("sdimx", "sdimy")], function(v) floor(range(v)))
    mid <- function(v, f) 2 * round((v[1] + f * diff(v)) / 2)
    rects <- data.frame(
        poly_ID = c("a", "b", "c"),
        xmin = c(mid(r$sdimx, 0), mid(r$sdimx, 0.3), mid(r$sdimx, 0.7)),
        xmax = c(mid(r$sdimx, 0.5), mid(r$sdimx, 0.6), mid(r$sdimx, 1) + 2),
        ymin = c(mid(r$sdimy, 0), mid(r$sdimy, 0.2), mid(r$sdimy, 0)),
        ymax = c(mid(r$sdimy, 0.6), mid(r$sdimy, 0.8), mid(r$sdimy, 0.4))
    )
    g <- .q(GiottoClass::setPolygonInfo(g,
        GiottoClass::createGiottoPolygon(.bp_rects(rects), name = "cell")))
    g <- GiottoClass::setFeatureInfo(g, bp, feat_type = "rna",
                                     verbose = FALSE)
    expect_s4_class(GiottoClass::getFeatureInfo(g, feat_type = "rna",
        return_giottoPoints = TRUE), "giottoBinPoints")
    expect_true(.gbp_store_backed(GiottoClass::getFeatureInfo(g, feat_type = "rna",
        return_giottoPoints = TRUE)))

    g <- .q(GiottoClass::calculateOverlap(g, spat_info = "cell",
                                          feat_info = "rna"))
    g <- .q(GiottoClass::overlapToMatrix(g, poly_info = "cell",
                                         feat_info = "rna", name = "raw"))
    got <- storeRead(GiottoClass::getExpression(g, spat_unit = "cell",
        feat_type = "rna", output = "exprObj")[], output = "dgcmatrix")

    m2 <- fx$m2[pe@feat_ids, ]
    inside <- function(k) {
        sl$cell_ID[sl$sdimx > rects$xmin[k] & sl$sdimx < rects$xmax[k] &
                   sl$sdimy > rects$ymin[k] & sl$sdimy < rects$ymax[k]]
    }
    exp <- vapply(seq_len(nrow(rects)), function(k) {
        Matrix::rowSums(m2[, inside(k), drop = FALSE])
    }, numeric(nrow(m2)))
    dimnames(exp) <- list(rownames(m2), rects$poly_ID)
    expect_true(all(lengths(lapply(seq_len(3), inside)) > 0L))
    expect_gt(length(intersect(inside(1), inside(2))), 0L)
    expect_identical(colnames(got), rects$poly_ID)
    expect_equal(as.matrix(got[rownames(exp), ]), exp)
})

test_that("disk bin points take the bin size", {
    fx <- .write_visiumhd_fixture()
    withr::local_options(list(tilework.warn_sequential = FALSE))
    b8 <- file.path(fx$root, "binned_outputs", "square_008um")
    bp <- .q(importVisiumHDDisk(file.path(fx$root, "binned_outputs"),
        backend = withr::local_tempdir(), bin = 8)$load_transcripts(
        bin = 8, expr_path = b8, tissue_positions_path = b8,
        scalefactors_path = b8, verbose = FALSE))[[1L]]
    expect_s4_class(bp, "giottoBinPoints")
    expect_true(.gbp_store_backed(bp))
    expect_equal(nrow(bp), sum(Matrix::colSums(fx$m8) > 0))
    tot <- dplyr::collect(dplyr::summarise(storeRead(bp@counts),
                                           x = sum(value)))$x
    expect_equal(tot, sum(fx$m8))
})

test_that("a crop view and feature selection narrow disk bin points", {
    fx <- .write_visiumhd_fixture()
    withr::local_options(list(tilework.warn_sequential = FALSE))
    proj <- withr::local_tempdir()
    g <- .q(importVisiumHDDisk(file.path(fx$root, "binned_outputs"),
        backend = proj, bin = 2)$create_gobject(load_transcripts = TRUE,
                                                load_image = FALSE,
                                                verbose = FALSE))
    sl <- GiottoClass::getSpatialLocations(g, spat_unit = "bin002",
                                           output = "data.table")
    r <- lapply(sl[, c("sdimx", "sdimy")], range)
    e <- terra::ext(r$sdimx[1] - 0.5, mean(r$sdimx), r$sdimy[1] - 0.5,
                    mean(r$sdimy))
    feats <- rownames(fx$m2)[1:2]
    g <- GiottoClass::crop(g, e, geom = "poly", view = "roi")
    got <- GiottoClass::combineFeatureData(g, feat_type = "rna",
        spat_unit = NULL, sel_feats = list(rna = feats), view = "roi")[[1]]

    inside <- sl[sdimx >= e$xmin & sdimx <= e$xmax &
                 sdimy >= e$ymin & sdimy <= e$ymax]
    m <- Matrix::summary(fx$m2[feats, inside$cell_ID, drop = FALSE])
    exp <- data.table::data.table(feat_ID = feats[m$i], count = m$x,
        x = inside$sdimx[m$j], y = inside$sdimy[m$j])
    key <- function(d) {
        d <- d[, list(feat_ID, count = as.numeric(count), x, y)]
        data.table::setorder(d, feat_ID, x, y)[]
    }
    expect_gt(nrow(exp), 0L)
    expect_equal(key(got), key(exp), ignore_attr = TRUE)
})

test_that("subsetGiotto() narrows disk bin points by feature and extent", {
    fx <- .write_visiumhd_fixture()
    withr::local_options(list(tilework.warn_sequential = FALSE))
    g <- .q(importVisiumHDDisk(file.path(fx$root, "binned_outputs"),
        backend = withr::local_tempdir(), bin = 2)$create_gobject(
        load_transcripts = TRUE, load_image = FALSE, verbose = FALSE))
    bp <- GiottoClass::getFeatureInfo(g, return_giottoPoints = TRUE)
    e <- terra::ext(bp)
    xmid <- mean(c(e$xmin, e$xmax))
    feats <- featIDs(bp)[1:2]
    g2 <- .q(GiottoClass::subsetGiotto(g, spat_unit = "bin002",
        feat_ids = feats))
    g2 <- .q(GiottoClass::subsetGiottoLocs(g2, spat_unit = "bin002",
        x_max = xmid))
    bp2 <- GiottoClass::getFeatureInfo(g2, return_giottoPoints = TRUE)
    expect_s4_class(bp2, "giottoBinPoints")
    expect_true(.gbp_store_backed(bp2))
    d <- data.table::as.data.table(bp2)
    expect_setequal(unique(d$feat_ID), intersect(feats, unique(d$feat_ID)))
    expect_true(all(d$feat_ID %in% feats))
    expect_true(all(d$x <= xmid))
})
