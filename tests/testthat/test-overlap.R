# calculateOverlap / overlapToMatrix on disk-backed stores. The overlap is a
# queryableStore over parquet with the fixed schema poly_ID, feat_ID,
# pt_tile_index, pt_row_index; every engine must find the same pairs.

options("tilework.warn_sequential" = FALSE)

# 3x3 square cells, side 8. Centres sit 0.3 off the point-tile edges at 20:
# the terra point-tile path double counts a cell whose centroid lies exactly
# on a tile edge, which is a separate bug. Points on a .5-offset unit lattice
# with five features, so none lies on a cell edge.
.ov_fx <- local({
    fx <- NULL
    function() {
        if (!is.null(fx)) return(fx)
        sq <- expand.grid(cx = c(10, 20, 30) + 0.3, cy = c(10, 20, 30) + 0.3)
        wkt <- sprintf(
            "POLYGON ((%1$s %2$s, %3$s %2$s, %3$s %4$s, %1$s %4$s, %1$s %2$s))",
            sq$cx - 4, sq$cy - 4, sq$cx + 4, sq$cy + 4)
        polys <- terra::vect(wkt, crs = "")
        polys$poly_ID <- paste0("p", seq_len(nrow(sq)))
        g <- expand.grid(x = seq(0.5, 39.5, by = 1), y = seq(0.5, 39.5, by = 1))
        g$feat_ID <- c("a", "b", "c", "d", "e")[(seq_len(nrow(g)) %% 5) + 1]
        pts_flat <- parquetGeomStore() |>
            storeWrite(terra::vect(g, geom = c("x", "y"), crs = ""))
        truth <- do.call(rbind, lapply(seq_len(nrow(sq)), function(k) {
            inb <- abs(g$x - sq$cx[k]) < 4 & abs(g$y - sq$cy[k]) < 4
            data.frame(poly_ID = paste0("p", k), feat_ID = g$feat_ID[inb])
        }))
        fx <<- list(
            polys = parquetGeomStore() |> storeWrite(polys),
            pts_flat = pts_flat,
            pts_tile = parquetGeomTileStore() |>
                storeWrite(pts_flat, threshold = 400L),
            truth = truth,
            cells = paste0("p", seq_len(nrow(sq))),
            feats = c("a", "b", "c", "d", "e")
        )
        fx
    }
})

# feature x cell counts from overlap rows, over the given universes
.ov_counts <- function(d, feats, cells) {
    d <- as.data.frame(d)
    d <- d[d$feat_ID %in% feats & d$poly_ID %in% cells, ]
    m <- matrix(0, length(feats), length(cells), dimnames = list(feats, cells))
    agg <- stats::aggregate(list(n = rep(1, nrow(d))),
        by = list(f = d$feat_ID, p = d$poly_ID), FUN = sum)
    m[cbind(agg$f, agg$p)] <- agg$n
    m
}

.ov_dense <- function(pe) {
    as.matrix(storeRead(pe, output = "dgcmatrix",
        max_rows = Inf, max_cols = Inf))
}

.skip_engine <- function(engine) {
    if (engine == "duckdb") {
        skip_if_not_installed("duckdb")
        skip_if_not_installed("dbplyr")
    }
    if (engine == "sedona") {
        skip_if_not_installed("sedonadb", minimum_version = "0.4.0")
    }
}

.engines <- list(
    terra_flat = list(engine = "terra", pts = "pts_flat"),
    terra_tile = list(engine = "terra", pts = "pts_tile"),
    duckdb     = list(engine = "duckdb", pts = "pts_tile"),
    sedona     = list(engine = "sedona", pts = "pts_tile")
)

describe("calculateOverlap engines", {
    for (nm in names(.engines)) {
        cfg <- .engines[[nm]]
        test_that(sprintf("%s: the carrier holds exactly the true pairs", nm), {
            .skip_engine(cfg$engine)
            fx <- .ov_fx()
            ov <- calculateOverlap(fx$polys, fx[[cfg$pts]],
                engine = cfg$engine)
            expect_s4_class(ov, "overlapPointDisk")
            expect_identical(class(ov@data)[[1L]], "queryableStore")
            d <- storeRead(ov@data, output = "tibble")
            expect_setequal(names(d),
                c("poly_ID", "feat_ID", "pt_tile_index", "pt_row_index"))
            expect_type(d$poly_ID, "character")
            expect_type(d$pt_row_index, "integer")
            # the key is unique: each point lands in at most one cell here
            expect_false(anyDuplicated(
                d[c("pt_tile_index", "pt_row_index", "poly_ID")]) > 0L)
            expect_equal(.ov_counts(d, fx$feats, fx$cells),
                .ov_counts(fx$truth, fx$feats, fx$cells))
        })
    }

    for (engine in c("duckdb", "sedona")) {
        test_that(sprintf("%s: pending filters on both stores apply", engine), {
            .skip_engine(engine)
            fx <- .ov_fx()
            polys <- subset(fx$polys, poly_ID != "p5")
            pts <- subset(fx$pts_tile, feat_ID != "a")
            d <- storeRead(calculateOverlap(polys, pts, engine = engine)@data,
                output = "tibble")
            expect_false("p5" %in% d$poly_ID)
            expect_false("a" %in% d$feat_ID)
            keep <- fx$truth$poly_ID != "p5" & fx$truth$feat_ID != "a"
            expect_equal(nrow(d), sum(keep))
        })
    }

    test_that("an overlap with no hits still reads, as zero rows", {
        fx <- .ov_fx()
        # two cells: the terra path cannot plan tiles for a single polygon
        far <- terra::vect(c(
            "POLYGON ((100 100, 101 100, 101 101, 100 101, 100 100))",
            "POLYGON ((200 200, 201 200, 201 201, 200 201, 200 200))"),
            crs = "")
        far$poly_ID <- c("far1", "far2")
        ov <- calculateOverlap(parquetGeomStore() |> storeWrite(far),
            fx$pts_flat)
        expect_equal(nrow(storeRead(ov@data, output = "tibble")), 0L)
    })
})

describe("calculateOverlap engine resolution", {
    auto_pick <- if (.spat_engine_available("sedonadb")) "sedona" else
        if (.spat_engine_available("duckdb")) "duckdb" else "terra"
    with_opt <- function(val, expr) GiottoUtils::gwith_options(
        list(giottodisk.spatial_query_engine = val), expr)

    test_that("an explicit engine wins, then the option, then auto", {
        expect_identical(with_opt("duckdb", .resolve_overlap_engine("terra")),
            "terra")
        expect_identical(with_opt("duckdb", .resolve_overlap_engine()), "duckdb")
        expect_identical(with_opt("auto", suppressMessages(
            .resolve_overlap_engine())), auto_pick)
        expect_identical(with_opt(NULL, suppressMessages(
            .resolve_overlap_engine())), auto_pick)
        expect_error(.resolve_overlap_engine("geos"), "should be one of")
    })

    test_that("auto picks terra when a tiling param is supplied", {
        expect_identical(with_opt("auto",
            .resolve_overlap_engine(tiling = TRUE)), "terra")
        # a named engine is honoured; it ignores the tiling params
        expect_identical(with_opt("auto",
            .resolve_overlap_engine("duckdb", tiling = TRUE)), "duckdb")
        expect_identical(with_opt("sedona",
            .resolve_overlap_engine(tiling = TRUE)), "sedona")
    })

    test_that("a tiling param reaches the resolver from the method", {
        skip_if_not_installed("duckdb")
        skip_if_not_installed("dbplyr")
        fx <- .ov_fx()
        files <- function(ov) basename(list.files(ov@data@path,
            pattern = "[.]parquet$", recursive = TRUE))
        # terra writes tile_NNNN.parquet files, duckdb one overlap.parquet
        with_opt("auto", {
            tiled <- calculateOverlap(fx$polys, fx$pts_tile, prune_tiles = TRUE)
        })
        expect_true(all(grepl("^tile_[0-9]+[.]parquet$", files(tiled))))
        with_opt("duckdb", {
            untiled <- calculateOverlap(fx$polys, fx$pts_tile)
        })
        expect_identical(files(untiled), "overlap.parquet")
    })
})

describe("overlapToMatrix", {
    test_that("the parquetExpr build matches the truth, cell-major", {
        skip_if_not_installed("duckdb")
        fx <- .ov_fx()
        ov <- calculateOverlap(fx$polys, fx$pts_flat)
        pe <- overlapToMatrix(ov, path = tempfile())
        expect_s4_class(pe, "parquetExprStore")
        .expect_cell_major(pe@path)
        M <- .ov_dense(pe)
        expect_equal(M, .ov_counts(fx$truth, fx$feats, fx$cells)[
            rownames(M), colnames(M)])
    })

    test_that("the duckdb and arrow builds give the same store", {
        skip_if_not_installed("duckdb")
        fx <- .ov_fx()
        ov <- calculateOverlap(fx$polys, fx$pts_flat)
        via_duckdb <- overlapToMatrix(ov, path = tempfile())
        via_arrow <- overlapToMatrix(ov@data, path = tempfile(),
            all_feat_ids = ov@feat_ids, all_cell_ids = ov@spat_ids)
        expect_identical(via_duckdb@cell_ids, via_arrow@cell_ids)
        expect_identical(via_duckdb@feat_ids, via_arrow@feat_ids)
        expect_identical(.ov_dense(via_duckdb), .ov_dense(via_arrow))
    })

    test_that("narrowed ID universes drop the rest of the overlap", {
        skip_if_not_installed("duckdb")
        fx <- .ov_fx()
        ov <- calculateOverlap(fx$polys, fx$pts_flat)
        ov@spat_ids <- c("p1", "p2")
        ov@feat_ids <- c("b", "c")
        pe <- overlapToMatrix(ov, path = tempfile())
        expect_identical(pe@cell_ids, c("p1", "p2"))
        expect_identical(pe@feat_ids, c("b", "c"))
        M <- .ov_dense(pe)
        expect_equal(M, .ov_counts(fx$truth, c("b", "c"), c("p1", "p2"))[
            rownames(M), colnames(M)])
        # the arrow build narrows the same way
        M2 <- .ov_dense(overlapToMatrix(ov@data, path = tempfile(),
            all_feat_ids = ov@feat_ids, all_cell_ids = ov@spat_ids))
        expect_identical(M2, M)
    })
})

test_that("snapshotLoad drops overlaps saved with a parquetStore carrier", {
    v <- terra::vect(c("POLYGON ((0 0, 2 0, 2 2, 0 2, 0 0))",
        "POLYGON ((3 3, 5 3, 5 5, 3 5, 3 3))"), crs = "")
    v$poly_ID <- c("p1", "p2")
    gp <- GiottoClass::createGiottoPolygon(v, name = "cell", verbose = FALSE)
    mk <- function(data) new("overlapPointDisk", data = data,
        poly_id_col = "poly_ID", feat_id_col = "feat_ID",
        spat_ids = c("p1", "p2"), feat_ids = "a")
    gp@overlaps <- list(
        rna = mk(parquetStore() |>
            storeWrite(data.frame(poly_ID = "p1", feat_ID = "a"))),
        protein = mk(.overlap_store(tempfile()))
    )
    g <- GiottoClass::setGiotto(GiottoClass::giotto(), gp, verbose = FALSE)

    td <- tempfile("ovload_")
    dir.create(file.path(td, "giottosave"), recursive = TRUE)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    saveRDS(g, file.path(td, "giottosave", "snap.rds"))

    expect_warning(loaded <- snapshotLoad(gDirSource(td)),
        "dropped overlaps saved by an older GiottoDisk: cell/rna")
    expect_identical(names(loaded@spatial_info$cell@overlaps), "protein")
})

test_that("the carrier adopts into a project vault and still reads", {
    fx <- .ov_fx()
    ov <- calculateOverlap(fx$polys, fx$pts_flat)
    n <- nrow(storeRead(ov@data, output = "tibble"))
    td <- tempfile("ovsrc_")
    dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    src <- gDirSource(td)
    old_path <- ov@data@path
    adopted <- sourceAdopt(src, ov@data,
        depends = c(ov@poly_uids, ov@feat_uids))
    expect_false(dir.exists(old_path))
    expect_true(sourceContains(src, adopted))
    expect_identical(adopted@uid, ov@data@uid)
    expect_equal(nrow(storeRead(adopted, output = "tibble")), n)
})
