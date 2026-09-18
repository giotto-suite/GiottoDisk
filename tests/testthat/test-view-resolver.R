# Tests for parquetCoordinator-side view resolution helpers and dispatch.
#
# Covers:
#   .find_store_with_cols
#   .narrow_dt_via_arrow      (same-store, cross-store DT owner)
#   .narrow_store_by_predicate (same-store, cross-store DT owner)
#   resolveSubobject(cellMetaObj/spatLocsObj, parquetCoordinator)
#
# Fixture: GiottoData::loadGiottoMini("visium") provides a gobject whose
# cell_metadata has `leiden_clus` (in_tissue, nr_feats, perc_feats,
# total_expr, custom_leiden) and whose spatial_locs has sdimx/sdimy
# keyed by cell_ID. Loaded once per file via setup chunk.

skip_if_no_mini <- function() {
    skip_if_not_installed("GiottoData")
    skip_if_not_installed("arrow")
}

.mini_g <- local({
    g_cache <- NULL
    function() {
        if (!is.null(g_cache)) return(g_cache)
        g_cache <<- GiottoData::loadGiottoMini("visium", verbose = FALSE)
        g_cache
    }
})


# Recipe builders ---------------------------------------------------------
#
# A resolver test needs a recipe with no gobject in scope, and recording
# normally happens on a gobject under a name. These assemble one
# directly -- but every STEP is built by the same public builder verb the
# recorder routes through, so a step-shape change upstream breaks these
# tests instead of silently drifting past them.

.mk_view <- function(...) {
    methods::new("giottoView", steps = list(...))
}

# A frame with no sample identity -- one set of steps applied to whatever
# subobject it is handed. That is a `perSampleSpace`; the `":default:"`
# sample key it used to be built with is gone, because the case it stood
# for is now the class rather than a magic key.
.mk_space <- function(...) {
    methods::new("perSampleSpace", name = "s", steps = list(...))
}

# filter step from an unevaluated predicate
.vfilter <- function(pred, ...) {
    subset(.mk_view(), subset = substitute(pred), quote = FALSE, ...)[[1L]]
}

# filter step from an already-built language object
.vfilter_lang <- function(pred, ...) {
    subset(.mk_view(), subset = pred, quote = FALSE, ...)[[1L]]
}

# crop step. `region` goes through the same WKT normalisation the
# recorder applies, and `geom` is the declared cell representation --
# "centroid" reduces to a cell_ID set, "poly" is evaluated on the cell
# polygon.
.vcrop <- function(region, relation = "intersects", geom = "centroid",
                   space = NULL) {
    crop(.mk_view(), region, relation = relation, geom = geom,
        space = space)[[1L]]
}

.stransform <- function(op, ...) {
    sp <- methods::new("perSampleSpace", name = "s")
    do.call(op, c(list(x = sp), list(...)))[[1L]]
}

# Cell_IDs inside `box`, namespaced `<sample>::<local>` the way the joint
# cell vocabulary of a giottoMulti is.
.ns_ids <- function(sl_list, box) {
    unique(unlist(lapply(names(sl_list), function(nm) {
        sd <- sl_list[[nm]]@coordinates
        ids <- sd$cell_ID[sd$sdimx >= box[1] & sd$sdimx <= box[2] &
                          sd$sdimy >= box[3] & sd$sdimy <= box[4]]
        if (length(ids) == 0L) return(character())
        paste(nm, ids, sep = "::")
    })))
}


# .find_store_with_cols ---------------------------------------------------

test_that(".find_store_with_cols: locates cellMeta col as in-mem data.table", {
    skip_if_no_mini()
    g <- .mini_g()
    hit <- .find_store_with_cols(g, "leiden_clus")
    expect_false(is.null(hit))
    expect_identical(hit$kind, "data.table")
    expect_identical(hit$key, "cell_ID")
    expect_true(data.table::is.data.table(hit$source))
    expect_true("leiden_clus" %in% names(hit$source))
})

test_that(".find_store_with_cols: unknown column returns NULL", {
    skip_if_no_mini()
    g <- .mini_g()
    expect_null(.find_store_with_cols(g, "does_not_exist"))
})

test_that(".find_store_with_cols: prefers cellMeta over spatLocs when both have cell_ID", {
    skip_if_no_mini()
    g <- .mini_g()
    # cell_ID exists on both cellMeta and spatLocs; cellMeta wins per
    # spatValues precedence order.
    hit <- .find_store_with_cols(g, "cell_ID")
    expect_identical(hit$key, "cell_ID")
    expect_true("nr_feats" %in% names(hit$source))  # cellMeta cols
})


# .narrow_dt_via_arrow ----------------------------------------------------

test_that(".narrow_dt_via_arrow: same-store filter on cellMeta", {
    skip_if_no_mini()
    g <- .mini_g()
    cm_dt <- GiottoClass::getCellMetadata(g, output = "data.table")
    pred  <- quote(leiden_clus == 1)
    out   <- .narrow_dt_via_arrow(cm_dt, pred, gobject = g, key = "cell_ID")
    expect_true(data.table::is.data.table(out))
    expect_true(nrow(out) > 0L)
    expect_true(all(out$leiden_clus == 1))
    expect_setequal(out$cell_ID, cm_dt$cell_ID[cm_dt$leiden_clus == 1])
})

test_that(".narrow_dt_via_arrow: cross-store DT-owner narrows spatLocs by cellMeta col", {
    skip_if_no_mini()
    g <- .mini_g()
    sl_dt <- GiottoClass::getSpatialLocations(g, output = "data.table")
    cm_dt <- GiottoClass::getCellMetadata(g, output = "data.table")
    pred  <- quote(leiden_clus == 1)
    out   <- .narrow_dt_via_arrow(sl_dt, pred, gobject = g, key = "cell_ID")
    expect_true(data.table::is.data.table(out))
    expect_setequal(names(out), names(sl_dt))
    expect_setequal(out$cell_ID, cm_dt$cell_ID[cm_dt$leiden_clus == 1])
})

test_that(".narrow_dt_via_arrow: errors when predicate cols are unknown", {
    skip_if_no_mini()
    g <- .mini_g()
    sl_dt <- GiottoClass::getSpatialLocations(g, output = "data.table")
    pred  <- quote(does_not_exist == 1)
    expect_error(
        .narrow_dt_via_arrow(sl_dt, pred, gobject = g, key = "cell_ID"),
        "no subobject covers predicate cols"
    )
})


# .narrow_store_by_predicate ----------------------------------------------

test_that(".narrow_store_by_predicate: same-store filter queues a filter op", {
    dt <- data.table::data.table(
        cell_ID = paste0("c", 1:10),
        cluster = rep(c("A", "B"), each = 5L),
        score   = seq_len(10L)
    )
    ps <- parquetStore() |> storeWrite(dt)
    pred <- quote(cluster == "A")
    out  <- .narrow_store_by_predicate(ps, pred, gobject = NULL)

    expect_true(inherits(out, "parquetBase"))
    types <- vapply(out@ops, `[[`, character(1L), "type")
    expect_true("filter" %in% types)

    res <- storeRead(out, output = "tibble")
    expect_equal(nrow(res), 5L)
    expect_true(all(res$cluster == "A"))
})

test_that(".narrow_store_by_predicate: cross-store DT-owner queues id_filter", {
    skip_if_no_mini()
    g <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")

    target_dt <- data.table::data.table(
        cell_ID = cmeta$cell_ID,
        value = seq_len(nrow(cmeta))
    )
    ps <- parquetStore() |> storeWrite(target_dt)

    pred <- quote(leiden_clus == 1)
    out  <- .narrow_store_by_predicate(ps, pred, gobject = g)

    expect_true(inherits(out, "parquetBase"))
    types <- vapply(out@ops, `[[`, character(1L), "type")
    expect_true("id_filter" %in% types)

    res <- storeRead(out, output = "tibble")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(res$cell_ID, expected)
})

test_that(".narrow_store_by_predicate: predicate cols across multiple stores fails cleanly", {
    skip_if_no_mini()
    g <- .mini_g()
    target_dt <- data.table::data.table(cell_ID = paste0("c", 1:3))
    ps <- parquetStore() |> storeWrite(target_dt)
    # leiden_clus is on cellMeta, sdimx on spatLocs -- different owners
    pred <- quote(leiden_clus == 1 & sdimx > 0)
    expect_error(
        .narrow_store_by_predicate(ps, pred, gobject = g),
        "no subobject covers"
    )
})


# resolveSubobject (parquetCoordinator) -----------------------------------

test_that("resolveSubobject(cellMetaObj, parquetCoordinator): empty view returns unchanged", {
    skip_if_no_mini()
    g <- .mini_g()
    cm <- GiottoClass::getCellMetadata(g, output = "cellMetaObj", copy_obj = TRUE)
    out <- resolveSubobject(cm, gobject = g,
        view = NULL, space = NULL,
        coordinator = parquetCoordinator())
    expect_identical(out@metaDT, cm@metaDT)
})

test_that("resolveSubobject(cellMetaObj, parquetCoordinator): same-store filter narrows DT", {
    skip_if_no_mini()
    g <- .mini_g()
    cm <- GiottoClass::getCellMetadata(g, output = "cellMetaObj", copy_obj = TRUE)
    v  <- .mk_view(.vfilter(leiden_clus == 1))
    out <- resolveSubobject(cm, gobject = g,
        view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_true(data.table::is.data.table(out@metaDT))
    expect_true(all(out@metaDT$leiden_clus == 1))
})

test_that("resolveSubobject(spatLocsObj, parquetCoordinator): cross-store filter via cellMeta", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "spatLocsObj",
        copy_obj = TRUE)
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    v  <- .mk_view(.vfilter(leiden_clus == 1))
    out <- resolveSubobject(sl, gobject = g,
        view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_true(data.table::is.data.table(out@coordinates))
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(out@coordinates$cell_ID, expected)
})


# Backed giottoPolygon / giottoPoints fixtures ----------------------------

.mk_backed_polygon <- function() {
    # poly_ID + region attribute. No geom needed for filter-op queueing —
    # the slot is ANY so a plain parquetStore is accepted by direct construction.
    dt <- data.table::data.table(
        poly_ID = paste0("c", 1:10),
        region  = rep(c("tumor", "stroma"), 5L),
        score   = seq_len(10L)
    )
    ps <- parquetStore() |> storeWrite(dt)
    new("giottoPolygon", spatVector = ps, name = "cell")
}

.mk_backed_points <- function(n = 20L) {
    dt <- data.table::data.table(
        cell_ID = rep(paste0("c", seq_len(n / 2L)), each = 2L),
        feature_name = rep(c("GENE1", "GENE2"), times = n / 2L),
        qv = stats::runif(n, 10, 40)
    )
    ps <- parquetStore() |> storeWrite(dt)
    new("giottoPoints", spatVector = ps, feat_type = "rna")
}


# resolveSubobject(giottoPolygon, parquetCoordinator) ---------------------

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): empty view returns unchanged", {
    gp <- .mk_backed_polygon()
    out <- resolveSubobject(gp, gobject = NULL,
        view = NULL, space = NULL,
        coordinator = parquetCoordinator())
    expect_identical(out@spatVector, gp@spatVector)
})

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): same-store filter queues filter op", {
    gp <- .mk_backed_polygon()
    v  <- .mk_view(.vfilter(region == "tumor"))
    out <- resolveSubobject(gp, gobject = NULL,
        view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_true(inherits(out@spatVector, "parquetBase"))
    types <- vapply(out@spatVector@ops, `[[`, character(1L), "type")
    expect_true("filter" %in% types)

    res <- storeRead(out@spatVector, output = "tibble")
    expect_equal(nrow(res), 5L)
    expect_true(all(res$region == "tumor"))
})

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): cross-store DT-owner queues id_filter with named-by", {
    skip_if_no_mini()
    g  <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")

    # Build a backed polygon whose poly_IDs match cell_IDs from the mini fixture
    poly_dt <- data.table::data.table(
        poly_ID = cmeta$cell_ID,
        score   = seq_len(nrow(cmeta))
    )
    ps <- parquetStore() |> storeWrite(poly_dt)
    gp <- new("giottoPolygon", spatVector = ps, name = "cell")

    v  <- .mk_view(.vfilter(leiden_clus == 1))
    out <- resolveSubobject(gp, gobject = g,
        view = v, space = NULL,
        coordinator = parquetCoordinator())

    types <- vapply(out@spatVector@ops, `[[`, character(1L), "type")
    expect_true("id_filter" %in% types)

    # The queued op should have named-by mapping poly_ID -> cell_ID
    id_op <- out@spatVector@ops[[which(types == "id_filter")[1L]]]
    expect_named(id_op$by, "poly_ID")
    expect_identical(unname(id_op$by), "cell_ID")

    res <- storeRead(out@spatVector, output = "tibble")
    expected_ids <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(res$poly_ID, expected_ids)
})


# resolveSubobject(giottoPoints, parquetCoordinator) ----------------------

test_that("resolveSubobject(giottoPoints, parquetCoordinator): same-store filter on points queues filter op", {
    gp <- .mk_backed_points(20L)
    v  <- .mk_view(.vfilter(feature_name == "GENE1"))
    out <- resolveSubobject(gp, gobject = NULL,
        view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_true(inherits(out@spatVector, "parquetBase"))
    types <- vapply(out@spatVector@ops, `[[`, character(1L), "type")
    expect_true("filter" %in% types)

    res <- storeRead(out@spatVector, output = "tibble")
    expect_true(all(res$feature_name == "GENE1"))
})


# viewCrop on backed parquetGeomStore targets -----------------------------

.mk_backed_geom_points <- function(n = 10L) {
    # 10 points at (1,1)..(10,10), each with cell_ID + feature_name
    sv <- terra::vect(
        data.frame(
            x = seq_len(n), y = seq_len(n),
            cell_ID = paste0("c", seq_len(n)),
            feature_name = rep(c("GENE1", "GENE2"), length.out = n)
        ),
        geom = c("x", "y"), crs = ""
    )
    ps <- parquetGeomStore() |> storeWrite(sv)
    new("giottoPoints", spatVector = ps, feat_type = "rna")
}

.mk_backed_geom_polygon <- function(n = 10L) {
    # use the same coords as backed_geom_points so cells/polys align
    sv <- terra::vect(
        data.frame(
            x = seq_len(n), y = seq_len(n),
            poly_ID = paste0("c", seq_len(n)),
            region = rep(c("tumor", "stroma"), length.out = n)
        ),
        geom = c("x", "y"), crs = ""
    )
    ps <- parquetGeomStore() |> storeWrite(sv)
    new("giottoPolygon", spatVector = ps, name = "cell")
}

test_that("resolveSubobject(giottoPoints, parquetCoordinator): viewCrop queues spat_relate; results narrowed by AABB", {
    gp <- .mk_backed_geom_points(10L)
    # Crop to points (3,3)..(6,6) — should keep c3, c4, c5, c6
    v  <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5)))
    out <- resolveSubobject(gp, gobject = NULL,
        view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_true(inherits(out@spatVector, "parquetGeomBase"))
    types <- vapply(out@spatVector@ops, `[[`, character(1L), "type")
    expect_true("spat_relate" %in% types)
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("cell_ID", "feature_name"))
    expect_setequal(res$cell_ID, c("c3", "c4", "c5", "c6"))
})

# A gobject whose polygon source IS the backed store, which is what a
# `geom = "poly"` crop evaluates against. Also the fixture the pre-rework
# suite never had: a backed polygon source paired with a tabular target.
.mk_backed_poly_gobject <- function(n = 10L) {
    gp <- .mk_backed_geom_polygon(n)
    g <- GiottoClass::giotto()
    g <- GiottoClass::setPolygonInfo(g, gp, name = "cell",
        centroids_to_spatlocs = TRUE, verbose = FALSE, initialize = FALSE)
    g
}

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): geom = 'poly' crop narrows by cell_ID like every cell-keyed target", {
    # adr/0015: a crop resolves to a surviving cell_ID set, and the
    # polygon store consumes it as one `id_filter` -- the same way cell
    # metadata and expression do. It does NOT get its own pushdown arm:
    # that was a second evaluation path for one predicate, and a recipe
    # must not mean different things depending on which slot reads it.
    g  <- .mk_backed_poly_gobject(10L)
    gp <- GiottoClass::getPolygonInfo(g, return_giottoPolygon = TRUE)
    v  <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"))
    out <- resolveSubobject(gp, gobject = g,
        view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_true(inherits(out@spatVector, "parquetGeomBase"))
    types <- vapply(out@spatVector@ops, `[[`, character(1L), "type")
    expect_true("id_filter" %in% types)
    expect_false("spat_relate" %in% types)
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("poly_ID", "region"))
    expect_setequal(res$poly_ID, c("c3", "c4", "c5", "c6"))
})

test_that("a poly crop gives a tabular target the same ID set as the polygon store", {
    # The regression the pre-rework suite could not catch: no fixture
    # paired a backed polygon source with a tabular target, so the
    # `terra::relate(<store>, ...)` dispatch failure never fired.
    g <- .mk_backed_poly_gobject(10L)
    target <- data.table::data.table(
        cell_ID = paste0("c", seq_len(10L)), value = seq_len(10L))
    ps <- parquetStore() |> storeWrite(target)
    v  <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"))

    out <- .push_view_to_pstore(ps, v, gobject = g,
        coordinator = parquetCoordinator())
    res <- storeRead(out, output = "tibble")
    expect_setequal(res$cell_ID, c("c3", "c4", "c5", "c6"))
})

test_that("the poly arm gives the same ID set on every engine", {
    # `geom` picks the geometry; `engine` picks the evaluator. They are
    # independent, so pinning the engine must not move the answer.
    g <- .mk_backed_poly_gobject(10L)
    v <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"))
    engines <- c("terra",
        if (requireNamespace("duckdb", quietly = TRUE)) "duckdb",
        if (requireNamespace("sedonadb", quietly = TRUE)) "sedona")

    answers <- lapply(engines, function(e) {
        GiottoUtils::gwith_options(
            list(giottodisk.spatial_query_engine = e),
            sort(dplyr::collect(.surviving_cell_ids_arrow(v, g,
                new.env(parent = emptyenv()),
                parquetCoordinator()))$cell_ID)
        )
    })
    expect_setequal(answers[[1L]], c("c3", "c4", "c5", "c6"))
    for (a in answers[-1L]) expect_identical(a, answers[[1L]])
})

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): a centroid crop with no centroid source warns", {
    gp <- .mk_backed_geom_polygon(10L)
    # `geom = "centroid"` asks about the cell's spatial_locs row. With no
    # gobject there is no centroid source, so it warns rather than
    # answering the polygon question in its place.
    v  <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5)))
    expect_warning(
        out <- resolveSubobject(gp, gobject = NULL,
            view = v, space = NULL,
            coordinator = parquetCoordinator()),
        "no spatial locations available"
    )
    types <- vapply(out@spatVector@ops, `[[`, character(1L), "type")
    expect_false("spat_relate" %in% types)
})

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): filter + geom crop compose", {
    g  <- .mk_backed_poly_gobject(10L)
    gp <- GiottoClass::getPolygonInfo(g, return_giottoPolygon = TRUE)
    # Crop to (2.5, 6.5, 2.5, 6.5) keeps c3-c6; further filter to "tumor"
    # (every other index: c1, c3, c5, c7, c9) -> intersection: c3, c5
    v  <- .mk_view(
        .vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"),
        .vfilter(region == "tumor"))
    out <- resolveSubobject(gp, gobject = g,
        view = v, space = NULL,
        coordinator = parquetCoordinator())
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("poly_ID", "region"))
    expect_setequal(res$poly_ID, c("c3", "c5"))
    expect_true(all(res$region == "tumor"))
})

test_that("push_view_to_pstore: a crop on a non-geom store reduces to an id_filter", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    xr <- range(sl$sdimx); yr <- range(sl$sdimy)
    box <- c(mean(xr) - 200, mean(xr) + 200, mean(yr) - 200, mean(yr) + 200)

    target <- data.table::data.table(cell_ID = sl$cell_ID, v = 1L)
    ps <- parquetStore() |> storeWrite(target)
    v  <- .mk_view(.vcrop(box))
    out <- .push_view_to_pstore(ps, v, gobject = g,
        coordinator = parquetCoordinator())
    types <- vapply(out@ops, `[[`, character(1L), "type")
    expect_true("id_filter" %in% types)

    res <- storeRead(out, output = "tibble")
    expected <- sl$cell_ID[sl$sdimx >= box[1] & sl$sdimx <= box[2] &
                           sl$sdimy >= box[3] & sl$sdimy <= box[4]]
    expect_setequal(res$cell_ID, expected)
})

test_that("push_view_to_pstore: a crop with nothing to evaluate it against warns", {
    dt <- data.table::data.table(
        cell_ID = paste0("c", 1:5), value = seq_len(5L)
    )
    ps <- parquetStore() |> storeWrite(dt)  # plain, no geom
    v  <- .mk_view(.vcrop(c(0, 100, 0, 100)))
    expect_warning(
        out <- .push_view_to_pstore(ps, v, gobject = NULL,
            coordinator = parquetCoordinator()),
        "no spatial locations available"
    )
    # Store unchanged
    expect_equal(length(out@ops), 0L)
})

test_that("push_view_to_pstore: non-intersects relation queues spat_relate op", {
    gp <- .mk_backed_geom_points(5L)
    v  <- .mk_view(.vcrop(c(0, 100, 0, 100), relation = "within"))
    # `cell_keyed = FALSE`: one row is one transcript, so a crop means
    # "clip these points" and is evaluated on the points' own geometry
    # regardless of the step's `geom` (which describes a CELL).
    out <- .push_view_to_pstore(gp@spatVector, v, gobject = NULL,
        coordinator = parquetCoordinator(), cell_keyed = FALSE)
    types <- vapply(out@ops, `[[`, character(1L), "type")
    expect_true("spat_relate" %in% types)
    relate_op <- out@ops[[which(types == "spat_relate")[1L]]]
    expect_identical(relate_op$relation, "within")
})

test_that("push_view_to_pstore: SpatVector region queues spat_relate carrying WKT", {
    gp <- .mk_backed_geom_points(10L)
    region <- terra::vect(
        "POLYGON((2.5 2.5, 6.5 2.5, 6.5 6.5, 2.5 6.5, 2.5 2.5))",
        crs = ""
    )
    v <- .mk_view(.vcrop(region))
    out <- .push_view_to_pstore(gp@spatVector, v, gobject = NULL,
        coordinator = parquetCoordinator(), cell_keyed = FALSE)
    types <- vapply(out@ops, `[[`, character(1L), "type")
    expect_true("spat_relate" %in% types)
    relate_op <- out@ops[[which(types == "spat_relate")[1L]]]
    expect_identical(relate_op$relation, "intersects")
    expect_type(relate_op$y_wkt, "character")
})

# spaceTransform pushdown -------------------------------------------------

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): space transforms compose into @post_ops", {
    gp <- .mk_backed_geom_polygon(5L)
    # Translate +5 +5; then spin 30
    sp <- .mk_space(
        .stransform("spatShift", dx = 5, dy = 5),
        .stransform("spin", 30))
    out <- resolveSubobject(gp, gobject = NULL,
        view = NULL, space = sp,
        coordinator = parquetCoordinator())
    expect_true(inherits(out@spatVector, "parquetGeomBase"))
    aff <- .pgeom_pending_transform(out@spatVector)
    expect_s4_class(aff, "affine2d")
    # After spatShift + spin, the affine should be non-identity.
    expect_false(isTRUE(all.equal(aff@affine, diag(3L))))
})

test_that("resolveSubobject(giottoPolygon, parquetCoordinator): space + view compose; both queue", {
    g  <- .mk_backed_poly_gobject(10L)
    gp <- GiottoClass::getPolygonInfo(g, return_giottoPolygon = TRUE)
    sp <- .mk_space(.stransform("spatShift", dx = 100, dy = 100))
    # The crop region is drawn in the view's OWN frame, which here is the
    # native one -- `space` is the output frame, applied to the returned
    # geometry, and the two are deliberately not the same thing.
    v  <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"))
    out <- resolveSubobject(gp, gobject = g,
        view = v, space = sp,
        coordinator = parquetCoordinator())
    # @post_ops carries the affine
    expect_s4_class(.pgeom_pending_transform(out@spatVector), "affine2d")
    # the crop queues an id_filter, like every cell-keyed target
    types <- vapply(out@spatVector@ops, `[[`, character(1L), "type")
    expect_true("id_filter" %in% types)
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("poly_ID", "region"))
    expect_setequal(res$poly_ID, c("c3", "c4", "c5", "c6"))
})

test_that("resolveSubobject(giottoPoints, parquetCoordinator): space transforms apply to backed @spatVector", {
    gp <- .mk_backed_geom_points(5L)
    sp <- .mk_space(.stransform("spatShift", dx = 10, dy = 10))
    out <- resolveSubobject(gp, gobject = NULL,
        view = NULL, space = sp,
        coordinator = parquetCoordinator())
    expect_true(inherits(out@spatVector, "parquetGeomBase"))
    expect_s4_class(.pgeom_pending_transform(out@spatVector), "affine2d")
})

test_that("resolveSubobject(spatLocsObj, parquetCoordinator): space transforms apply to coordinates", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "spatLocsObj",
        copy_obj = TRUE)
    sp <- .mk_space(.stransform("spatShift", dx = 100, dy = 100))
    out <- resolveSubobject(sl, gobject = g,
        view = NULL, space = sp,
        coordinator = parquetCoordinator())
    # spatLocsObj's spatShift method mutates @coordinates eagerly
    expect_true(data.table::is.data.table(out@coordinates))
    expect_equal(out@coordinates$sdimx, sl@coordinates$sdimx + 100)
    expect_equal(out@coordinates$sdimy, sl@coordinates$sdimy + 100)
})

test_that(".apply_space_to_subobj: NULL space is a no-op", {
    gp <- .mk_backed_geom_polygon(3L)
    out <- .apply_space_to_subobj(gp, gobject = NULL, space = NULL)
    expect_identical(out, gp)
})

# Cache thread-through -----------------------------------------------------

test_that(".surviving_cell_ids_arrow: stores one entry per view, intersection across steps", {
    skip_if_no_mini()
    g  <- .mini_g()
    cache <- new.env(parent = emptyenv())
    v  <- .mk_view(.vfilter(leiden_clus == 1), .vfilter(in_tissue == 1))
    surv <- .surviving_cell_ids_arrow(v, g, cache, parquetCoordinator())
    expect_s3_class(surv, "Table")
    # ONE target-independent quantity (adr/0015): every cell-keyed slot
    # narrows by the same set, so there is nothing to split.
    expect_setequal(ls(cache), "surviving_cell_ids")

    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1 & cmeta$in_tissue == 1]
    expect_setequal(dplyr::collect(surv)$cell_ID, expected)
})

test_that(".surviving_cell_ids_arrow: cache hit returns the same object without recompute", {
    skip_if_no_mini()
    g  <- .mini_g()
    cache <- new.env(parent = emptyenv())
    v  <- .mk_view(.vfilter(leiden_clus == 1))
    surv1 <- .surviving_cell_ids_arrow(v, g, cache, parquetCoordinator())
    surv2 <- .surviving_cell_ids_arrow(v, g, cache, parquetCoordinator())
    expect_identical(surv1, surv2)
})

test_that("cache path: multiple targets sharing one view share one id_filter table", {
    skip_if_no_mini()
    g  <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    target_dt1 <- data.table::data.table(
        cell_ID = cmeta$cell_ID, value1 = seq_len(nrow(cmeta))
    )
    target_dt2 <- data.table::data.table(
        cell_ID = cmeta$cell_ID, value2 = seq_len(nrow(cmeta))
    )
    ps1 <- parquetStore() |> storeWrite(target_dt1)
    ps2 <- parquetStore() |> storeWrite(target_dt2)

    cache <- new.env(parent = emptyenv())
    v <- .mk_view(.vfilter(leiden_clus == 1))
    out1 <- .push_view_to_pstore(ps1, v, gobject = g,
        coordinator = parquetCoordinator(), .cache = cache)
    out2 <- .push_view_to_pstore(ps2, v, gobject = g,
        coordinator = parquetCoordinator(), .cache = cache)
    expect_setequal(ls(cache), "surviving_cell_ids")

    t1 <- vapply(out1@ops, `[[`, character(1L), "type")
    t2 <- vapply(out2@ops, `[[`, character(1L), "type")
    op1 <- out1@ops[[which(t1 == "id_filter")[1L]]]
    op2 <- out2@ops[[which(t2 == "id_filter")[1L]]]
    expect_identical(op1$ids_tab, op2$ids_tab)

    res1 <- storeRead(out1, output = "tibble")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(res1$cell_ID, expected)
})

test_that("cache path: viewCrop folds into surviving_cell_ids via spatial_locs", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    # Crop region around the central spot — picks cells whose centroid
    # is in the box.
    xrange <- range(sl$sdimx); yrange <- range(sl$sdimy)
    box <- c(mean(xrange) - 200, mean(xrange) + 200,
        mean(yrange) - 200, mean(yrange) + 200)

    v  <- .mk_view(.vcrop(box))
    cache <- new.env(parent = emptyenv())
    surv <- .surviving_cell_ids_arrow(v, g, cache, parquetCoordinator())
    expect_s3_class(surv, "Table")
    surv_ids <- dplyr::collect(surv)$cell_ID

    expected <- sl$cell_ID[sl$sdimx >= box[1] & sl$sdimx <= box[2] &
                           sl$sdimy >= box[3] & sl$sdimy <= box[4]]
    expect_setequal(surv_ids, expected)
})

test_that("cache path: viewFilter + viewCrop intersection in surviving_cell_ids", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    xrange <- range(sl$sdimx); yrange <- range(sl$sdimy)
    box <- c(mean(xrange) - 500, mean(xrange) + 500,
        mean(yrange) - 500, mean(yrange) + 500)

    v <- .mk_view(.vfilter(leiden_clus == 1), .vcrop(box))
    cache <- new.env(parent = emptyenv())
    surv <- .surviving_cell_ids_arrow(v, g, cache, parquetCoordinator())
    surv_ids <- dplyr::collect(surv)$cell_ID

    in_box <- sl$cell_ID[sl$sdimx >= box[1] & sl$sdimx <= box[2] &
                         sl$sdimy >= box[3] & sl$sdimy <= box[4]]
    in_clus <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(surv_ids, intersect(in_box, in_clus))
})

test_that("cache path: multi-filter view produces ONE id_filter op per target", {
    skip_if_no_mini()
    g  <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    target_dt <- data.table::data.table(cell_ID = cmeta$cell_ID, v = 1L)
    ps <- parquetStore() |> storeWrite(target_dt)
    v <- .mk_view(.vfilter(leiden_clus == 1), .vfilter(in_tissue == 1))

    cache <- new.env(parent = emptyenv())
    out <- .push_view_to_pstore(ps, v, gobject = g,
        coordinator = parquetCoordinator(), .cache = cache)
    types <- vapply(out@ops, `[[`, character(1L), "type")
    # exactly one id_filter representing the full intersection
    expect_equal(sum(types == "id_filter"), 1L)

    res <- storeRead(out, output = "tibble")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1 & cmeta$in_tissue == 1]
    expect_setequal(res$cell_ID, expected)
})

test_that(".cache: NULL cache preserves the non-cache branching path", {
    # With no cache, .narrow_store_by_predicate uses its old branches —
    # same-store fast path queues a `filter` op (not id_filter).
    dt <- data.table::data.table(
        cell_ID = paste0("c", 1:5), cluster = c("A","B","A","B","A")
    )
    ps <- parquetStore() |> storeWrite(dt)
    out <- .narrow_store_by_predicate(ps, quote(cluster == "A"),
        gobject = NULL)
    types <- vapply(out@ops, `[[`, character(1L), "type")
    expect_true("filter" %in% types)
    expect_false("id_filter" %in% types)
})

test_that("resolveSubobject(spatLocsObj) + cellMeta predicate: cache fills via push_view_to_dt", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "spatLocsObj",
        copy_obj = TRUE)
    v  <- .mk_view(.vfilter(leiden_clus == 1))
    cache <- new.env(parent = emptyenv())
    out <- resolveSubobject(sl, gobject = g, view = v, space = NULL,
        coordinator = parquetCoordinator(), .cache = cache)
    expect_true(data.table::is.data.table(out@coordinates))
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(out@coordinates$cell_ID, expected)
    expect_true("surviving_cell_ids" %in% ls(cache))
})


# exprObj / dimObj narrowing -----------------------------------------------

test_that("resolveSubobject(exprObj, parquetCoordinator): backed exprMat narrows by cell_ID via [,j]", {
    skip_if_no_mini()
    g  <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    cell_ids <- cmeta$cell_ID

    pe <- parquetExprStore(
        path     = tempfile(fileext = ".parquet"),
        cell_ids = cell_ids,
        feat_ids = c("g1", "g2", "g3")
    )
    eo <- GiottoClass::createExprObj(expression_data = pe, name = "raw",
        spat_unit = "cell", feat_type = "rna")

    v <- .mk_view(.vfilter(leiden_clus == 1))
    out <- resolveSubobject(eo, gobject = g, view = v, space = NULL,
        coordinator = parquetCoordinator())

    expect_s4_class(out@exprMat, "parquetExprStore")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(out@exprMat@cell_ids, expected)
    # @cell_idx records the original positions of the surviving cells
    expect_equal(length(out@exprMat@cell_idx), length(expected))
})

test_that("resolveSubobject(exprObj, parquetCoordinator): empty narrow yields zero-cell store", {
    skip_if_no_mini()
    g  <- .mini_g()
    pe <- parquetExprStore(
        path = tempfile(fileext = ".parquet"),
        cell_ids = c("not_in_mini_1", "not_in_mini_2"),
        feat_ids = "g1"
    )
    eo <- GiottoClass::createExprObj(expression_data = pe, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    v <- .mk_view(.vfilter(leiden_clus == 1))
    out <- resolveSubobject(eo, gobject = g, view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_s4_class(out@exprMat, "parquetExprStore")
    expect_equal(out@exprMat@n_cells, 0)
})

test_that("resolveSubobject(exprObj, parquetCoordinator): in-mem exprMat falls through", {
    skip_if_no_mini()
    g  <- .mini_g()
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    # The mini's expression is in-mem (matrix / dgCMatrix) — should
    # callNextMethod into dataTableCoordinator's narrow path.
    v <- .mk_view(.vfilter(leiden_clus == 1))
    out <- resolveSubobject(eo, gobject = g, view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_s4_class(out, "exprObj")
    # dataTableCoordinator's narrow on exprObj does column subsetting by
    # surviving cell_IDs (via .cached_surviving_cell_ids inside that method).
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(colnames(out@exprMat), expected)
})

test_that("resolveSubobject(exprObj, parquetCoordinator): empty view returns unchanged", {
    pe <- parquetExprStore(
        path = tempfile(fileext = ".parquet"),
        cell_ids = paste0("c", 1:5),
        feat_ids = c("g1", "g2")
    )
    eo <- GiottoClass::createExprObj(expression_data = pe, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    out <- resolveSubobject(eo, gobject = NULL,
        view = NULL, space = NULL,
        coordinator = parquetCoordinator())
    expect_identical(out@exprMat@cell_ids, eo@exprMat@cell_ids)
})

# Piece B: expression-value predicates ------------------------------------

.mk_g_with_backed_expr <- function() {
    # Wrap a tiny labeled dgCMatrix in a parquetExprStore and inject it
    # into a fresh giotto fixture's @expression slot.
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    cells <- cmeta$cell_ID
    set.seed(1L)
    n_genes <- 4L
    m <- Matrix::rsparsematrix(n_genes, length(cells), density = 0.5,
        rand.x = function(n) as.double(stats::rpois(n, 5L) + 1L))
    rownames(m) <- c("MYC", "IFN_R1", "TP53", "ACTB")
    colnames(m) <- cells

    pe <- storeWrite(parquetExprStore(
        path = tempfile(fileext = ".parquet")), m)
    eo <- GiottoClass::createExprObj(expression_data = pe, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    list(g = GiottoClass::setExpression(g, eo, verbose = FALSE),
        gene_mat = m)
}

test_that(".find_store_with_cols: locates gene names on parquetExprStore", {
    skip_if_no_mini()
    fx <- .mk_g_with_backed_expr()
    hit <- .find_store_with_cols(fx$g, "MYC")
    expect_false(is.null(hit))
    expect_identical(hit$kind, "data.table")
    expect_identical(hit$key, "cell_ID")
    expect_true("MYC" %in% names(hit$source))
    expect_true("cell_ID" %in% names(hit$source))
    # One row per cell in the fixture (gobject has N cells, slice has N rows)
    cmeta <- GiottoClass::getCellMetadata(fx$g, output = "data.table")
    expect_equal(nrow(hit$source), nrow(cmeta))
})

test_that(".find_store_with_cols: multi-gene predicate cols on parquetExprStore", {
    skip_if_no_mini()
    fx <- .mk_g_with_backed_expr()
    hit <- .find_store_with_cols(fx$g, c("MYC", "TP53"))
    expect_false(is.null(hit))
    expect_true(all(c("MYC", "TP53", "cell_ID") %in% names(hit$source)))
})

test_that("view filter on expression value: surviving cell_IDs match the manual filter", {
    skip_if_no_mini()
    fx <- .mk_g_with_backed_expr()
    # Compute the expected surviving cell_IDs manually
    myc_vals <- as.numeric(fx$gene_mat["MYC", ])
    expected_ids <- colnames(fx$gene_mat)[myc_vals > 0]

    v <- .mk_view(.vfilter(MYC > 0))
    cache <- new.env(parent = emptyenv())
    surv <- .surviving_cell_ids_arrow(v, fx$g, cache, parquetCoordinator())
    surv_ids <- dplyr::collect(surv)$cell_ID
    expect_setequal(surv_ids, expected_ids)
})

test_that(".find_store_with_cols: generic matrix-rownames branch picks up dgCMatrix expression", {
    skip_if_no_mini()
    g  <- .mini_g()
    # visium mini expression is a dgCMatrix in-mem; pick a real gene
    e  <- GiottoClass::getExpression(g, output = "exprObj")
    gene <- rownames(e@exprMat)[[1L]]
    hit <- .find_store_with_cols(g, gene)
    expect_false(is.null(hit))
    expect_identical(hit$kind, "data.table")
    expect_identical(hit$key, "cell_ID")
    expect_true(gene %in% names(hit$source))
    # Row count = n_cells of the underlying matrix
    expect_equal(nrow(hit$source), ncol(e@exprMat))
})

test_that("view filter on gene expression: in-mem dgCMatrix path matches manual subset", {
    skip_if_no_mini()
    g  <- .mini_g()
    e  <- GiottoClass::getExpression(g, output = "exprObj")
    gene <- rownames(e@exprMat)[[1L]]
    vals <- as.numeric(e@exprMat[gene, ])
    names(vals) <- colnames(e@exprMat)
    expected <- names(vals)[vals > 0]

    pred <- bquote(.(as.name(gene)) > 0)
    v <- .mk_view(.vfilter_lang(pred))
    cache <- new.env(parent = emptyenv())
    surv <- .surviving_cell_ids_arrow(v, g, cache, parquetCoordinator())
    surv_ids <- dplyr::collect(surv)$cell_ID
    expect_setequal(surv_ids, expected)
})

test_that("view filter on gene expression intersects with cellMeta predicate", {
    skip_if_no_mini()
    fx <- .mk_g_with_backed_expr()
    cmeta <- GiottoClass::getCellMetadata(fx$g, output = "data.table")
    myc_vals <- as.numeric(fx$gene_mat["MYC", ])
    expressing <- colnames(fx$gene_mat)[myc_vals > 0]
    clus1 <- cmeta$cell_ID[cmeta$leiden_clus == 1]

    v <- .mk_view(.vfilter(leiden_clus == 1), .vfilter(MYC > 0))
    cache <- new.env(parent = emptyenv())
    surv <- .surviving_cell_ids_arrow(v, fx$g, cache, parquetCoordinator())
    surv_ids <- dplyr::collect(surv)$cell_ID
    expect_setequal(surv_ids, intersect(expressing, clus1))
})


test_that("resolveSubobject(dimObj, parquetCoordinator): in-mem coordinates fall through", {
    skip_if_no_mini()
    g  <- .mini_g()
    dr <- GiottoClass::getDimReduction(g, reduction = "cells",
        reduction_method = "pca", name = "pca", output = "dimObj")
    if (is.null(dr)) skip("no PCA in mini fixture")
    v <- .mk_view(.vfilter(leiden_clus == 1))
    out <- resolveSubobject(dr, gobject = g, view = v, space = NULL,
        coordinator = parquetCoordinator())
    expect_s4_class(out, "dimObj")
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    # dataTableCoordinator narrows by intersect(rownames, keep)
    expect_setequal(rownames(out@coordinates), expected)
})


# giottoMulti verification --------------------------------------------------

.mk_multi <- function() {
    g1 <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    g2 <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    GiottoClass::createGiottoMulti(list(s1 = g1, s2 = g2))
}

test_that("multi: spatIDs works as the bootstrap for surviving_cell_ids", {
    skip_if_no_mini()
    mg <- .mk_multi()
    ids <- GiottoClass::spatIDs(mg)
    expect_type(ids, "character")
    expect_gt(length(ids), 0L)
})

test_that("multi: .find_store_with_cols hits joint cellMeta with list_ID column", {
    skip_if_no_mini()
    mg <- .mk_multi()
    hit <- .find_store_with_cols(mg, "leiden_clus")
    expect_identical(hit$kind, "data.table")
    expect_true("list_ID" %in% names(hit$source))
    expect_true("leiden_clus" %in% names(hit$source))
})

test_that("multi: viewFilter via cache surviving_cell_ids matches joint cellMeta filter", {
    skip_if_no_mini()
    mg <- .mk_multi()
    v <- .mk_view(.vfilter(leiden_clus == 1))
    cache <- new.env(parent = emptyenv())
    surv <- .surviving_cell_ids_arrow(v, mg, cache, parquetCoordinator())
    surv_ids <- dplyr::collect(surv)$cell_ID

    cm <- GiottoClass::getCellMetadata(mg, output = "data.table")
    expected <- unique(cm$cell_ID[cm$leiden_clus == 1])
    expect_setequal(surv_ids, expected)
})

test_that("multi: viewCrop unions cell_IDs across per-sample spatLocs", {
    skip_if_no_mini()
    mg <- .mk_multi()
    sl_list <- GiottoClass::getSpatialLocations(mg, output = "spatLocsObj")
    sd1 <- sl_list[[1L]]@coordinates
    box <- c(min(sd1$sdimx) + 1000, min(sd1$sdimx) + 3000,
        min(sd1$sdimy) + 1000, min(sd1$sdimy) + 3000)
    v <- .mk_view(.vcrop(box))
    surv <- .surviving_cell_ids_arrow(v, mg, new.env(parent = emptyenv()),
        parquetCoordinator())
    surv_ids <- dplyr::collect(surv)$cell_ID

    expected <- .ns_ids(sl_list, box)
    expect_setequal(surv_ids, expected)
})

test_that("multi: viewFilter + viewCrop intersect via the cache surviving set", {
    skip_if_no_mini()
    mg <- .mk_multi()
    sl_list <- GiottoClass::getSpatialLocations(mg, output = "spatLocsObj")
    sd1 <- sl_list[[1L]]@coordinates
    box <- c(min(sd1$sdimx), min(sd1$sdimx) + 4000,
        min(sd1$sdimy), min(sd1$sdimy) + 4000)
    v <- .mk_view(.vfilter(leiden_clus == 1), .vcrop(box))
    surv <- .surviving_cell_ids_arrow(v, mg, new.env(parent = emptyenv()),
        parquetCoordinator())
    surv_ids <- dplyr::collect(surv)$cell_ID

    in_box <- .ns_ids(sl_list, box)
    cm <- GiottoClass::getCellMetadata(mg, output = "data.table")
    clus1 <- unique(cm$cell_ID[cm$leiden_clus == 1])
    expect_setequal(surv_ids, intersect(in_box, clus1))
})


test_that(".apply_space_to_subobj: multi-step recipe composes via accumulating matrix", {
    gp <- .mk_backed_geom_polygon(3L)
    sp <- .mk_space(
        .stransform("spatShift", dx = 1, dy = 0),
        .stransform("spatShift", dx = 0, dy = 1),
        .stransform("spatShift", dx = 2, dy = 2))
    out <- .apply_space_to_subobj(gp, gobject = NULL, space = sp)
    aff <- .pgeom_pending_transform(out@spatVector)
    expect_s4_class(aff, "affine2d")
    # Net translation should be (3, 3). The translation entries of the
    # 3x3 affine are at [1,3] (x-shift) and [2,3] (y-shift).
    expect_equal(aff@affine[1L, 3L], 3)
    expect_equal(aff@affine[2L, 3L], 3)
})
