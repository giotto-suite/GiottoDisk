# Tests for analyzeData(parquetEdgeStore, labelProportionsParam).
#
# Strategy: build a small igraph + a parquetEdgeStore from the same edges,
# call analyzeData on both with the same labels DT and param, compare.
# The arrow/dplyr path should produce identical proportions to the in-mem
# igraph reference path in GiottoClass.

.lp_fixture <- function(directed = FALSE, with_weight = TRUE) {
    # 6-vertex network, two roughly disconnected components for variety.
    edges <- data.table::data.table(
        from   = c("a","a","b","b","c","d","d","e"),
        to     = c("b","c","c","d","d","e","f","f"),
        weight = c(0.9, 0.5, 0.6, 0.4, 0.8, 0.3, 0.7, 0.2)
    )
    if (!with_weight) edges[, weight := NULL]
    g <- igraph::graph_from_data_frame(edges, directed = directed)
    list(edges = edges, g = g)
}

.lp_labels <- function(g) {
    nm <- igraph::V(g)$name
    data.table::data.table(
        cell_ID  = nm,
        celltype = c("X","Y","X","Y","X","Y")[seq_along(nm)]
    )
}


test_that("disk method matches in-mem igraph reference (unweighted, alpha = 1)", {
    fx <- .lp_fixture(directed = FALSE, with_weight = TRUE)
    s  <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                     fx$g, type = "sNN")
    labs <- .lp_labels(fx$g)
    p <- GiottoClass::labelProportionsParam(
        labels = "celltype",
        group_method = "spatialnetwork",
        alpha = 1, weights = FALSE
    )

    ref <- suppressWarnings(GiottoClass::analyzeData(fx$g, p, labels = labs))
    got <- suppressWarnings(GiottoClass::analyzeData(s,    p, labels = labs))

    # Align row order (group col may be int-derived in disk path)
    ref <- ref[order(ref$group), ]
    got <- got[order(got$group), ]
    expect_equal(ref, got, tolerance = 1e-8)
})


test_that("disk method matches in-mem reference (weighted, alpha = 1)", {
    fx <- .lp_fixture(directed = FALSE, with_weight = TRUE)
    s  <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                     fx$g, type = "sNN")
    labs <- .lp_labels(fx$g)
    p <- GiottoClass::labelProportionsParam(
        labels = "celltype",
        group_method = "spatialnetwork",
        alpha = 1, weights = TRUE
    )

    ref <- GiottoClass::analyzeData(fx$g, p, labels = labs)
    got <- GiottoClass::analyzeData(s,    p, labels = labs)

    ref <- ref[order(ref$group), ]
    got <- got[order(got$group), ]
    expect_equal(ref, got, tolerance = 1e-8)
})


test_that("disk method matches reference (weighted, alpha = 0 → no self)", {
    fx <- .lp_fixture(directed = FALSE, with_weight = TRUE)
    s  <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                     fx$g, type = "sNN")
    labs <- .lp_labels(fx$g)
    p <- GiottoClass::labelProportionsParam(
        labels = "celltype",
        group_method = "spatialnetwork",
        alpha = 0, weights = TRUE
    )

    ref <- GiottoClass::analyzeData(fx$g, p, labels = labs)
    got <- GiottoClass::analyzeData(s,    p, labels = labs)

    ref <- ref[order(ref$group), ]
    got <- got[order(got$group), ]
    expect_equal(ref, got, tolerance = 1e-8)
})


test_that("disk method matches reference (weighted, alpha = 0.5)", {
    fx <- .lp_fixture(directed = FALSE, with_weight = TRUE)
    s  <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                     fx$g, type = "sNN")
    labs <- .lp_labels(fx$g)
    p <- GiottoClass::labelProportionsParam(
        labels = "celltype",
        group_method = "spatialnetwork",
        alpha = 0.5, weights = TRUE
    )

    ref <- GiottoClass::analyzeData(fx$g, p, labels = labs)
    got <- GiottoClass::analyzeData(s,    p, labels = labs)

    ref <- ref[order(ref$group), ]
    got <- got[order(got$group), ]
    expect_equal(ref, got, tolerance = 1e-8)
})


test_that("disk method errors without labels", {
    fx <- .lp_fixture()
    s  <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                     fx$g, type = "sNN")
    p <- GiottoClass::labelProportionsParam(
        labels = "celltype",
        group_method = "spatialnetwork"
    )
    expect_error(GiottoClass::analyzeData(s, p), "labels.*required")
})


test_that("parquetGeomBase stub errors with a clear not-yet-supported message", {
    # Construct a minimal parquetGeomStore handle (no on-disk files needed
    # — dispatch happens before any I/O).
    s <- parquetGeomStore(path = tempfile())
    p <- GiottoClass::labelProportionsParam(
        labels = "celltype", group_method = "polygon", spat_info = "x"
    )
    expect_error(GiottoClass::analyzeData(s, p),
                 "spatRelate.*not yet implemented")
})


test_that("disk method warns + falls back to adjacency when weights=TRUE but no weight col", {
    fx <- .lp_fixture(with_weight = FALSE)
    s  <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                     fx$g, type = "sNN")
    labs <- .lp_labels(fx$g)
    p <- GiottoClass::labelProportionsParam(
        labels = "celltype",
        group_method = "spatialnetwork",
        alpha = 1, weights = TRUE
    )
    expect_warning(GiottoClass::analyzeData(s, p, labels = labs),
                   "falling back to adjacency")
})
