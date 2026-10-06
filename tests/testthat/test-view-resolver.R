# Tests for parquetCoordinator-side view resolution.
#
# Covers:
#   .find_store_with_cols
#   GiottoClass::resolveKeep(parquetCoordinator): the surviving set as a lazy plan
#   GiottoClass::resolveRecipe(<leaf>, parquetCoordinator): backed leaves queue it; in-memory
#     leaves fall through to dataTableCoordinator and read `vector`
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


# Resolve one leaf the way the container does: evaluate the view into the
# op's surviving set once, then hand the leaf that set, the recipe and the
# frames -- never the gobject.
.rk <- function(x, g, v, space = NULL, spaces = list(),
                co = parquetCoordinator()) {
    keep <- if (is.null(v)) NULL else GiottoClass::resolveKeep(co, g, v)
    GiottoClass::resolveRecipe(x, co, keep = keep, view = v, space = space, spaces = spaces)
}

.op_types <- function(store) vapply(store@ops, `[[`, character(1L), "type")

.keep_ids <- function(g, v) {
    sort(dplyr::collect(GiottoClass::resolveKeep(parquetCoordinator(), g, v)$arrow)$cell_ID)
}


# Backed giottoPolygon / giottoPoints fixtures ----------------------------

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

# A gobject whose polygon source IS the backed store, which is what a
# `geom = "poly"` crop evaluates against -- and a backed polygon source
# paired with a tabular target.
.mk_backed_poly_gobject <- function(n = 10L) {
    gp <- .mk_backed_geom_polygon(n)
    g <- GiottoClass::giotto()
    g <- GiottoClass::setPolygonInfo(g, gp, name = "cell",
        centroids_to_spatlocs = TRUE, verbose = FALSE, initialize = FALSE)
    g
}


# resolveKeep (parquetCoordinator) -----------------------------------------

test_that("resolveKeep: NULL when the view narrows nothing", {
    expect_null(GiottoClass::resolveKeep(parquetCoordinator(), NULL, NULL))
    expect_null(GiottoClass::resolveKeep(parquetCoordinator(), NULL, .mk_view()))
})

test_that("resolveKeep: the arrow form is a lazy plan, vector a promise", {
    # the backed polygon store is the filter owner here, so its surviving
    # key column is contributed as a query rather than read
    g <- .mk_backed_poly_gobject(10L)
    v <- .mk_view(.vfilter(region == "tumor"))
    expected <- paste0("c", c(1, 3, 5, 7, 9))

    k <- GiottoClass::resolveKeep(parquetCoordinator(), g, v)
    expect_s3_class(k, "viewKeep")
    expect_s3_class(k$arrow, "arrow_dplyr_query")
    expect_setequal(dplyr::collect(k$arrow)$cell_ID, expected)
    expect_setequal(k$vector, expected)

    # Nothing is read until a form is used: with the owner's files gone,
    # building the set still succeeds, and only reading it fails.
    g2 <- .mk_backed_poly_gobject(10L)
    k2 <- GiottoClass::resolveKeep(parquetCoordinator(), g2, v)
    gp2 <- GiottoClass::getPolygonInfo(g2, return_giottoPolygon = TRUE)
    unlink(storePaths(gp2@spatVector), recursive = TRUE)
    expect_error(k2$vector)
})

test_that("resolveKeep: filters and crops intersect, in step order", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    xr <- range(sl$sdimx); yr <- range(sl$sdimy)
    box <- c(mean(xr) - 500, mean(xr) + 500, mean(yr) - 500, mean(yr) + 500)

    v <- .mk_view(.vfilter(leiden_clus == 1), .vfilter(in_tissue == 1),
        .vcrop(box))
    in_box <- sl$cell_ID[sl$sdimx >= box[1] & sl$sdimx <= box[2] &
                         sl$sdimy >= box[3] & sl$sdimy <= box[4]]
    in_clus <- cmeta$cell_ID[cmeta$leiden_clus == 1 & cmeta$in_tissue == 1]
    expect_setequal(.keep_ids(g, v), intersect(in_box, in_clus))
})

test_that("resolveKeep: a filter nothing covers is an error", {
    skip_if_no_mini()
    expect_error(GiottoClass::resolveKeep(parquetCoordinator(), .mini_g(),
        .mk_view(.vfilter(does_not_exist == 1))),
        "no subobject covers predicate cols")
    # leiden_clus is on cellMeta, sdimx on spatLocs -- different owners
    expect_error(GiottoClass::resolveKeep(parquetCoordinator(), .mini_g(),
        .mk_view(.vfilter(leiden_clus == 1 & sdimx > 0))),
        "no subobject covers")
})

test_that("resolveKeep: a crop with nothing to evaluate it against warns", {
    v <- .mk_view(.vcrop(c(0, 100, 0, 100)))
    expect_warning(k <- GiottoClass::resolveKeep(parquetCoordinator(),
        GiottoClass::giotto(), v), "no spatial locations available")
    expect_null(k)
})


# GiottoClass::resolveRecipe(cellMetaObj / spatLocsObj / spatEnrObj) ---------------------------

test_that("GiottoClass::resolveRecipe(cellMetaObj): keep = NULL returns it unchanged", {
    skip_if_no_mini()
    cm <- GiottoClass::getCellMetadata(.mini_g(), output = "cellMetaObj",
        copy_obj = TRUE)
    out <- GiottoClass::resolveRecipe(cm, parquetCoordinator())
    expect_identical(out@metaDT, cm@metaDT)
})

test_that("GiottoClass::resolveRecipe(cellMetaObj): in-memory metadata falls through and narrows", {
    skip_if_no_mini()
    g  <- .mini_g()
    cm <- GiottoClass::getCellMetadata(g, output = "cellMetaObj",
        copy_obj = TRUE)
    out <- .rk(cm, g, .mk_view(.vfilter(leiden_clus == 1)))
    expect_true(data.table::is.data.table(out@metaDT))
    expect_true(all(out@metaDT$leiden_clus == 1))
})

test_that("resolve: a multi-step view queues ONE id_filter on a backed target", {
    skip_if_no_mini()
    g <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    ps <- parquetStore() |> storeWrite(
        data.table::data.table(cell_ID = cmeta$cell_ID, v = 1L))
    v <- .mk_view(.vfilter(leiden_clus == 1), .vfilter(in_tissue == 1))
    out <- .queue_keep(ps, GiottoClass::resolveKeep(parquetCoordinator(), g, v))
    expect_identical(.op_types(out), "id_filter")
    expect_setequal(storeRead(out, output = "tibble")$cell_ID,
        cmeta$cell_ID[cmeta$leiden_clus == 1 & cmeta$in_tissue == 1])
})

test_that("GiottoClass::resolveRecipe(spatLocsObj): in-memory coordinates narrow by a cellMeta filter", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "spatLocsObj",
        copy_obj = TRUE)
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    out <- .rk(sl, g, .mk_view(.vfilter(leiden_clus == 1)))
    expect_true(data.table::is.data.table(out@coordinates))
    expect_setequal(out@coordinates$cell_ID,
        cmeta$cell_ID[cmeta$leiden_clus == 1])
})

test_that("resolve: one keep narrows every backed target to the same set", {
    skip_if_no_mini()
    g <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    ps1 <- parquetStore() |> storeWrite(data.table::data.table(
        cell_ID = cmeta$cell_ID, v1 = seq_len(nrow(cmeta))))
    ps2 <- parquetStore() |> storeWrite(data.table::data.table(
        cell_ID = cmeta$cell_ID, v2 = seq_len(nrow(cmeta))))
    k <- GiottoClass::resolveKeep(parquetCoordinator(), g,
        .mk_view(.vfilter(leiden_clus == 1)))
    out1 <- .queue_keep(ps1, k)
    out2 <- .queue_keep(ps2, k)
    expect_identical(out1@ops[[1L]]$ids_tab, out2@ops[[1L]]$ids_tab)
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(storeRead(out1, output = "tibble")$cell_ID, expected)
    expect_setequal(storeRead(out2, output = "tibble")$cell_ID, expected)
})

test_that("resolve: a crop reduces to an id_filter on a non-geom target", {
    skip_if_no_mini()
    g  <- .mini_g()
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    xr <- range(sl$sdimx); yr <- range(sl$sdimy)
    box <- c(mean(xr) - 200, mean(xr) + 200, mean(yr) - 200, mean(yr) + 200)
    ps <- parquetStore() |> storeWrite(
        data.table::data.table(cell_ID = sl$cell_ID, v = 1L))
    out <- .queue_keep(ps, GiottoClass::resolveKeep(parquetCoordinator(), g,
        .mk_view(.vcrop(box))))
    expect_identical(.op_types(out), "id_filter")
    expected <- sl$cell_ID[sl$sdimx >= box[1] & sl$sdimx <= box[2] &
                           sl$sdimy >= box[3] & sl$sdimy <= box[4]]
    expect_setequal(storeRead(out, output = "tibble")$cell_ID, expected)
})


# GiottoClass::resolveRecipe(giottoPolygon) -----------------------------------------------------

test_that("GiottoClass::resolveRecipe(giottoPolygon): keep = NULL returns it unchanged", {
    gp <- .mk_backed_geom_polygon()
    out <- GiottoClass::resolveRecipe(gp, parquetCoordinator())
    expect_identical(out@spatVector, gp@spatVector)
})

test_that("GiottoClass::resolveRecipe(giottoPolygon): a filter on its own attributes narrows it", {
    # The polygon store is its own filter owner: poly_ID answers for the
    # cell axis, so `region` resolves there like any metadata column.
    g  <- .mk_backed_poly_gobject(10L)
    gp <- GiottoClass::getPolygonInfo(g, return_giottoPolygon = TRUE)
    out <- .rk(gp, g, .mk_view(.vfilter(region == "tumor")))
    expect_identical(.op_types(out@spatVector), "id_filter")
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("poly_ID", "region"))
    expect_setequal(res$poly_ID, paste0("c", c(1, 3, 5, 7, 9)))
})

test_that("GiottoClass::resolveRecipe(giottoPolygon): a cellMeta filter queues id_filter mapping poly_ID onto cell_ID", {
    skip_if_no_mini()
    g  <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    ps <- parquetStore() |> storeWrite(data.table::data.table(
        poly_ID = cmeta$cell_ID, score = seq_len(nrow(cmeta))))
    gp <- new("giottoPolygon", spatVector = ps, name = "cell")

    out <- .rk(gp, g, .mk_view(.vfilter(leiden_clus == 1)))
    id_op <- out@spatVector@ops[[1L]]
    expect_identical(id_op$type, "id_filter")
    expect_named(id_op$by, "poly_ID")
    expect_identical(unname(id_op$by), "cell_ID")
    res <- storeRead(out@spatVector, output = "tibble")
    expect_setequal(res$poly_ID, cmeta$cell_ID[cmeta$leiden_clus == 1])
})

test_that("GiottoClass::resolveRecipe(giottoPolygon): geom = 'poly' crop narrows by cell_ID like every cell-keyed target", {
    # adr/0015: a crop resolves to a surviving cell_ID set, and the polygon
    # store consumes it as one `id_filter` -- the same way cell metadata and
    # expression do. No separate pushdown arm: a recipe must not mean
    # different things depending on which slot reads it.
    g  <- .mk_backed_poly_gobject(10L)
    gp <- GiottoClass::getPolygonInfo(g, return_giottoPolygon = TRUE)
    out <- .rk(gp, g, .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly")))
    expect_identical(.op_types(out@spatVector), "id_filter")
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("poly_ID", "region"))
    expect_setequal(res$poly_ID, c("c3", "c4", "c5", "c6"))
})

test_that("a poly crop gives a tabular target the same ID set as the polygon store", {
    g <- .mk_backed_poly_gobject(10L)
    ps <- parquetStore() |> storeWrite(data.table::data.table(
        cell_ID = paste0("c", seq_len(10L)), value = seq_len(10L)))
    k <- GiottoClass::resolveKeep(parquetCoordinator(), g,
        .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly")))
    res <- storeRead(.queue_keep(ps, k), output = "tibble")
    expect_setequal(res$cell_ID, c("c3", "c4", "c5", "c6"))
})

test_that("the poly arm gives the same ID set on every engine", {
    # `geom` picks the geometry; `engine` picks the evaluator. They are
    # independent, so pinning the engine must not move the answer.
    g <- .mk_backed_poly_gobject(10L)
    v <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"))
    engines <- c("terra",
        if (requireNamespace("duckdb", quietly = TRUE)) "duckdb",
        if (.spat_engine_available("sedonadb")) "sedona")
    answers <- lapply(engines, function(e) {
        GiottoUtils::gwith_options(
            list(giottodisk.spatial_query_engine = e), .keep_ids(g, v))
    })
    expect_setequal(answers[[1L]], c("c3", "c4", "c5", "c6"))
    for (a in answers[-1L]) expect_identical(a, answers[[1L]])
})

test_that("GiottoClass::resolveRecipe(giottoPolygon): filter + geom crop compose", {
    g  <- .mk_backed_poly_gobject(10L)
    gp <- GiottoClass::getPolygonInfo(g, return_giottoPolygon = TRUE)
    # crop keeps c3-c6; "tumor" is every other index -> c3, c5
    v  <- .mk_view(
        .vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"),
        .vfilter(region == "tumor"))
    res <- storeRead(.rk(gp, g, v)@spatVector, output = "tibble",
        fields = c("poly_ID", "region"))
    expect_setequal(res$poly_ID, c("c3", "c5"))
})


# GiottoClass::resolveRecipe(giottoPoints) ------------------------------------------------------

test_that("GiottoClass::resolveRecipe(giottoPoints): a crop clips the store's own geometry", {
    gp <- .mk_backed_geom_points(10L)
    out <- .rk(gp, NULL, .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5))))
    expect_identical(.op_types(out@spatVector), "spat_relate")
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("cell_ID", "feature_name"))
    expect_setequal(res$cell_ID, c("c3", "c4", "c5", "c6"))
})

test_that("GiottoClass::resolveRecipe(giottoPoints): the step's relation and WKT region are carried", {
    gp <- .mk_backed_geom_points(5L)
    out <- .rk(gp, NULL, .mk_view(.vcrop(c(0, 100, 0, 100),
        relation = "within")))
    op <- out@spatVector@ops[[1L]]
    expect_identical(op$relation, "within")
    expect_type(op$y_wkt, "character")
})

test_that("GiottoClass::resolveRecipe(giottoPoints): filter steps are skipped, as in memory", {
    # filters are cell-centric; one row here is one transcript
    gp <- .mk_backed_geom_points(10L)
    out <- GiottoClass::resolveRecipe(gp, parquetCoordinator(),
        view = .mk_view(.vfilter(feature_name == "GENE1")))
    expect_length(out@spatVector@ops, 0L)
})

test_that("GiottoClass::resolveRecipe(giottoPoints): a crop drawn in a named space clips in that frame", {
    # The crop is drawn in "shifted", +100 in x. The points come back
    # native, so the region is projected back before it clips.
    gp <- .mk_backed_geom_points(10L)
    shifted <- .mk_space(.stransform("spatShift", dx = 100, dy = 0))
    v <- .mk_view(.vcrop(c(102.5, 106.5, 2.5, 6.5), space = "shifted"))
    out <- .rk(gp, NULL, v, spaces = list(shifted = shifted))
    res <- storeRead(out@spatVector, output = "tibble", fields = "cell_ID")
    expect_setequal(res$cell_ID, c("c3", "c4", "c5", "c6"))
    expect_error(.rk(gp, NULL, v), "drawn in space 'shifted'")
})


# spaceTransform pushdown -------------------------------------------------

test_that("GiottoClass::resolveRecipe(giottoPolygon): space transforms compose into @post_ops", {
    gp <- .mk_backed_geom_polygon(5L)
    sp <- .mk_space(
        .stransform("spatShift", dx = 5, dy = 5),
        .stransform("spin", 30))
    out <- GiottoClass::resolveRecipe(gp, parquetCoordinator(), space = sp)
    aff <- .pgeom_pending_transform(out@spatVector)
    expect_s4_class(aff, "affine2d")
    expect_false(isTRUE(all.equal(aff@affine, diag(3L))))
})

test_that("GiottoClass::resolveRecipe(giottoPolygon): space + view compose; both queue", {
    g  <- .mk_backed_poly_gobject(10L)
    gp <- GiottoClass::getPolygonInfo(g, return_giottoPolygon = TRUE)
    sp <- .mk_space(.stransform("spatShift", dx = 100, dy = 100))
    # The crop region is drawn in the view's own frame, the native one here;
    # `space` is the output frame applied to the returned geometry.
    v  <- .mk_view(.vcrop(c(2.5, 6.5, 2.5, 6.5), geom = "poly"))
    out <- .rk(gp, g, v, space = sp)
    expect_s4_class(.pgeom_pending_transform(out@spatVector), "affine2d")
    expect_identical(.op_types(out@spatVector), "id_filter")
    res <- storeRead(out@spatVector, output = "tibble",
        fields = c("poly_ID", "region"))
    expect_setequal(res$poly_ID, c("c3", "c4", "c5", "c6"))
})

test_that("GiottoClass::resolveRecipe(giottoPoints): space transforms apply to backed @spatVector", {
    gp <- .mk_backed_geom_points(5L)
    sp <- .mk_space(.stransform("spatShift", dx = 10, dy = 10))
    out <- GiottoClass::resolveRecipe(gp, parquetCoordinator(), space = sp)
    expect_s4_class(.pgeom_pending_transform(out@spatVector), "affine2d")
})

test_that("GiottoClass::resolveRecipe(spatLocsObj): space transforms apply to in-memory coordinates", {
    skip_if_no_mini()
    sl <- GiottoClass::getSpatialLocations(.mini_g(), output = "spatLocsObj",
        copy_obj = TRUE)
    sp <- .mk_space(.stransform("spatShift", dx = 100, dy = 100))
    out <- GiottoClass::resolveRecipe(sl, parquetCoordinator(), space = sp)
    expect_equal(out@coordinates$sdimx, sl@coordinates$sdimx + 100)
    expect_equal(out@coordinates$sdimy, sl@coordinates$sdimy + 100)
})

test_that(".apply_space_to_subobj: NULL space is a no-op", {
    gp <- .mk_backed_geom_polygon(3L)
    expect_identical(.apply_space_to_subobj(gp, NULL), gp)
})

test_that(".apply_space_to_subobj: multi-step recipe composes via accumulating matrix", {
    gp <- .mk_backed_geom_polygon(3L)
    sp <- .mk_space(
        .stransform("spatShift", dx = 1, dy = 0),
        .stransform("spatShift", dx = 0, dy = 1),
        .stransform("spatShift", dx = 2, dy = 2))
    aff <- .pgeom_pending_transform(.apply_space_to_subobj(gp, sp)@spatVector)
    expect_s4_class(aff, "affine2d")
    # net translation (3, 3), at [1,3] (x-shift) and [2,3] (y-shift)
    expect_equal(aff@affine[1L, 3L], 3)
    expect_equal(aff@affine[2L, 3L], 3)
})


# exprObj / dimObj narrowing -----------------------------------------------

test_that("GiottoClass::resolveRecipe(exprObj): backed exprMat narrows by cell_ID via [,j]", {
    skip_if_no_mini()
    g  <- .mini_g()
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    pe <- parquetExprStore(path = tempfile(fileext = ".parquet"),
        cell_ids = cmeta$cell_ID, feat_ids = c("g1", "g2", "g3"))
    eo <- GiottoClass::createExprObj(expression_data = pe, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    out <- .rk(eo, g, .mk_view(.vfilter(leiden_clus == 1)))
    expect_s4_class(out@exprMat, "parquetExprStore")
    expected <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    expect_setequal(out@exprMat@cell_ids, expected)
    expect_equal(length(out@exprMat@cell_idx), length(expected))
})

test_that("GiottoClass::resolveRecipe(exprObj): an empty narrow yields a zero-cell store", {
    skip_if_no_mini()
    pe <- parquetExprStore(path = tempfile(fileext = ".parquet"),
        cell_ids = c("not_in_mini_1", "not_in_mini_2"), feat_ids = "g1")
    eo <- GiottoClass::createExprObj(expression_data = pe, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    out <- .rk(eo, .mini_g(), .mk_view(.vfilter(leiden_clus == 1)))
    expect_equal(out@exprMat@n_cells, 0)
})

test_that("GiottoClass::resolveRecipe(exprObj): in-memory exprMat falls through and reads vector", {
    skip_if_no_mini()
    g  <- .mini_g()
    eo <- GiottoClass::getExpression(g, output = "exprObj")
    out <- .rk(eo, g, .mk_view(.vfilter(leiden_clus == 1)))
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    expect_setequal(colnames(out@exprMat),
        cmeta$cell_ID[cmeta$leiden_clus == 1])
})

test_that("GiottoClass::resolveRecipe(exprObj): keep = NULL returns it unchanged", {
    pe <- parquetExprStore(path = tempfile(fileext = ".parquet"),
        cell_ids = paste0("c", 1:5), feat_ids = c("g1", "g2"))
    eo <- GiottoClass::createExprObj(expression_data = pe, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    out <- GiottoClass::resolveRecipe(eo, parquetCoordinator())
    expect_identical(out@exprMat@cell_ids, eo@exprMat@cell_ids)
})

test_that("GiottoClass::resolveRecipe(dimObj): inherits the in-memory leaf", {
    skip_if_no_mini()
    g  <- .mini_g()
    dr <- GiottoClass::getDimReduction(g, reduction = "cells",
        reduction_method = "pca", name = "pca", output = "dimObj")
    if (is.null(dr)) skip("no PCA in mini fixture")
    out <- .rk(dr, g, .mk_view(.vfilter(leiden_clus == 1)))
    cmeta <- GiottoClass::getCellMetadata(g, output = "data.table")
    expect_setequal(rownames(out@coordinates),
        cmeta$cell_ID[cmeta$leiden_clus == 1])
})


# Expression-value predicates --------------------------------------------

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
    expect_identical(hit$kind, "data.table")
    expect_identical(hit$key, "cell_ID")
    expect_true(all(c("MYC", "cell_ID") %in% names(hit$source)))
    cmeta <- GiottoClass::getCellMetadata(fx$g, output = "data.table")
    expect_equal(nrow(hit$source), nrow(cmeta))
})

test_that(".find_store_with_cols: multi-gene predicate cols on parquetExprStore", {
    skip_if_no_mini()
    fx <- .mk_g_with_backed_expr()
    hit <- .find_store_with_cols(fx$g, c("MYC", "TP53"))
    expect_true(all(c("MYC", "TP53", "cell_ID") %in% names(hit$source)))
})

test_that(".find_store_with_cols: generic matrix-rownames branch picks up dgCMatrix expression", {
    skip_if_no_mini()
    g  <- .mini_g()
    e  <- GiottoClass::getExpression(g, output = "exprObj")
    gene <- rownames(e@exprMat)[[1L]]
    hit <- .find_store_with_cols(g, gene)
    expect_identical(hit$kind, "data.table")
    expect_true(gene %in% names(hit$source))
    expect_equal(nrow(hit$source), ncol(e@exprMat))
})

test_that("view filter on expression value: surviving cell_IDs match the manual filter", {
    skip_if_no_mini()
    fx <- .mk_g_with_backed_expr()
    myc_vals <- as.numeric(fx$gene_mat["MYC", ])
    expected <- colnames(fx$gene_mat)[myc_vals > 0]
    expect_setequal(.keep_ids(fx$g, .mk_view(.vfilter(MYC > 0))), expected)
})

test_that("view filter on gene expression: in-mem dgCMatrix path matches manual subset", {
    skip_if_no_mini()
    g  <- .mini_g()
    e  <- GiottoClass::getExpression(g, output = "exprObj")
    gene <- rownames(e@exprMat)[[1L]]
    vals <- as.numeric(e@exprMat[gene, ])
    expected <- colnames(e@exprMat)[vals > 0]
    v <- .mk_view(.vfilter_lang(bquote(.(as.name(gene)) > 0)))
    expect_setequal(.keep_ids(g, v), expected)
})

test_that("view filter on gene expression intersects with cellMeta predicate", {
    skip_if_no_mini()
    fx <- .mk_g_with_backed_expr()
    cmeta <- GiottoClass::getCellMetadata(fx$g, output = "data.table")
    expressing <- colnames(fx$gene_mat)[as.numeric(fx$gene_mat["MYC", ]) > 0]
    clus1 <- cmeta$cell_ID[cmeta$leiden_clus == 1]
    v <- .mk_view(.vfilter(leiden_clus == 1), .vfilter(MYC > 0))
    expect_setequal(.keep_ids(fx$g, v), intersect(expressing, clus1))
})


# giottoMulti verification --------------------------------------------------

.mk_multi <- function() {
    g1 <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    g2 <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    GiottoClass::createGiottoMulti(list(s1 = g1, s2 = g2))
}

test_that("multi: spatIDs works as the bootstrap for the surviving set", {
    skip_if_no_mini()
    ids <- GiottoClass::spatIDs(.mk_multi())
    expect_type(ids, "character")
    expect_gt(length(ids), 0L)
})

test_that("multi: .find_store_with_cols hits joint cellMeta with list_ID column", {
    skip_if_no_mini()
    hit <- .find_store_with_cols(.mk_multi(), "leiden_clus")
    expect_identical(hit$kind, "data.table")
    expect_true(all(c("list_ID", "leiden_clus") %in% names(hit$source)))
})

test_that("multi: a filter matches the joint cellMeta filter", {
    skip_if_no_mini()
    mg <- .mk_multi()
    cm <- GiottoClass::getCellMetadata(mg, output = "data.table")
    expect_setequal(.keep_ids(mg, .mk_view(.vfilter(leiden_clus == 1))),
        unique(cm$cell_ID[cm$leiden_clus == 1]))
})

test_that("multi: a crop unions cell_IDs across per-sample spatLocs", {
    skip_if_no_mini()
    mg <- .mk_multi()
    sl_list <- GiottoClass::getSpatialLocations(mg, output = "spatLocsObj")
    sd1 <- sl_list[[1L]]@coordinates
    box <- c(min(sd1$sdimx) + 1000, min(sd1$sdimx) + 3000,
        min(sd1$sdimy) + 1000, min(sd1$sdimy) + 3000)
    expect_setequal(.keep_ids(mg, .mk_view(.vcrop(box))),
        .ns_ids(sl_list, box))
})

test_that("multi: a filter and a crop intersect", {
    skip_if_no_mini()
    mg <- .mk_multi()
    sl_list <- GiottoClass::getSpatialLocations(mg, output = "spatLocsObj")
    sd1 <- sl_list[[1L]]@coordinates
    box <- c(min(sd1$sdimx), min(sd1$sdimx) + 4000,
        min(sd1$sdimy), min(sd1$sdimy) + 4000)
    cm <- GiottoClass::getCellMetadata(mg, output = "data.table")
    clus1 <- unique(cm$cell_ID[cm$leiden_clus == 1])
    v <- .mk_view(.vfilter(leiden_clus == 1), .vcrop(box))
    expect_setequal(.keep_ids(mg, v), intersect(.ns_ids(sl_list, box), clus1))
})
