# Tests for parquetEdgeStore: disk-backed network store.
# Locks down:
#   - storeCreate / round-trip via storeWrite dispatch on data.table,
#     igraph, and edge input markers
#   - three storeRead output modes (arrow / tibble / igraph)
#   - [i, j] subsetting semantics (induced subgraph, adjacency slice,
#     to-endpoint filter, named-args, negate); undirected orientation
#     symmetry and directed strict slicing
#   - node_universe preserves isolated vertices
#   - sourceAdopt fix-up of nested @nodes handle
#   - sourceWrite via gDirSource for both data.table and igraph

# ---- helpers --------------------------------------------------------------

.tiny_undirected_dt <- function() {
    # 5-vertex hand-coded undirected sNN-like edge table — canonical from <= to.
    data.table::data.table(
        from   = c("a","a","a","b","b","c","d"),
        to     = c("b","c","d","c","d","d","e"),
        weight = c(0.9, 0.7, 0.5, 0.6, 0.4, 0.8, 0.3)
    )
}

.tiny_directed_dt <- function() {
    data.table::data.table(
        from = c("a","a","b","c","c","d"),
        to   = c("b","c","c","a","d","a")
    )
}


# ---- Class basics --------------------------------------------------------

test_that("storeCreate dispatches parquetEdgeStore by type alias", {
    s1 <- storeCreate(type = "parquetEdgeStore")
    s2 <- storeCreate(type = "parquetEdge")
    expect_s4_class(s1, "parquetEdgeStore")
    expect_s4_class(s2, "parquetEdgeStore")
})

test_that("empty parquetEdgeStore is constructible without files on disk", {
    s <- parquetEdgeStore(path = tempfile())
    expect_s4_class(s, "parquetEdgeStore")
    expect_s4_class(s@nodes, "parquetStore")
    expect_false(storeExists(s))
    expect_false(storeExists(s@nodes))
})


# ---- storeWrite paths -----------------------------------------------------

test_that("storeWrite(parquetEdgeStore, data.table) writes both subdirs", {
    dt <- .tiny_undirected_dt()
    s <- storeCreate(type = "parquetEdgeStore")
    s <- storeWrite(s, dt, type = "sNN", directed = FALSE)

    expect_equal(s@n_cells, 5)
    expect_equal(s@n_edges, 7)
    expect_equal(s@type, "sNN")
    expect_false(s@directed)
    expect_true(file.exists(file.path(s@path, "edges", "edges.parquet")))
    expect_true(file.exists(file.path(s@path, "nodes", "nodes.parquet")))
})

test_that("storeWrite(parquetEdgeStore, igraph) extracts edges + vertex universe", {
    g <- igraph::make_graph(c("a","b", "b","c", "c","a", "d","e"),
                            directed = FALSE)
    igraph::E(g)$weight <- c(0.5, 0.7, 0.9, 0.3)
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), g, type = "sNN")
    expect_equal(s@n_cells, 5)        # a, b, c, d, e (including isolates none)
    expect_equal(s@n_edges, 4)
    # directed inferred from is_directed(g)
    expect_false(s@directed)
})

test_that("storeWrite picks up directed from igraph::is_directed()", {
    g <- igraph::make_graph(c("a","b", "b","c"), directed = TRUE)
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), g, type = "kNN")
    expect_true(s@directed)
})

test_that("node_universe arg preserves isolated vertices", {
    dt <- data.table::data.table(from = c("a","b"), to = c("b","c"))
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), dt,
                    type = "spatial", directed = FALSE,
                    node_universe = c("a","b","c","d","e"))
    # n_cells = full universe even though d/e have no edges
    expect_equal(s@n_cells, 5)
    nodes <- dplyr::collect(storeRead(s@nodes))
    expect_setequal(nodes$node_id, c("a","b","c","d","e"))
})

test_that("storeWrite via edgeDTInput escape hatch", {
    dt <- data.table::data.table(
        u = c("a","b","b"),
        v = c("b","c","a"),
        w = c(0.1, 0.2, 0.3)
    )
    inp <- edgeDTInput(dt, from_col = "u", to_col = "v",
                       weight_col = "w")
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), inp,
                    type = "kNN", directed = TRUE)
    expect_equal(s@n_edges, 3)
    r <- storeRead(s, output = "tibble")
    expect_setequal(names(r), c("from_id", "to_id", "weight"))
})


# ---- storeRead output modes -----------------------------------------------

test_that("storeRead arrow returns a queryable arrow dataset with raw ints", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    a <- storeRead(s, output = "arrow")
    r <- dplyr::collect(a)
    expect_true(is.integer(r$from_id) || is.numeric(r$from_id))
    expect_equal(nrow(r), 7L)
})

test_that("storeRead tibble joins char node IDs via sidecar", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    r <- storeRead(s, output = "tibble")
    expect_s3_class(r, "data.table")
    expect_true(is.character(r$from_id))
    expect_setequal(unique(c(r$from_id, r$to_id)), c("a","b","c","d","e"))
})

test_that("storeRead igraph constructs with char V()$name and int internals", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    g <- storeRead(s, output = "igraph")
    expect_s3_class(g, "igraph")
    expect_setequal(igraph::V(g)$name, c("a","b","c","d","e"))
    expect_equal(igraph::ecount(g), 7L)
    expect_false(igraph::is_directed(g))
})


# ---- Subsetting -----------------------------------------------------------

test_that("x[v_set] is induced subgraph (both endpoints in v_set)", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    sub <- s[c("a","b","c")]
    r <- storeRead(sub, output = "tibble")
    # only a-b, a-c, b-c qualify
    expect_equal(nrow(r), 3L)
    expect_true(all(r$from_id %in% c("a","b","c")))
    expect_true(all(r$to_id   %in% c("a","b","c")))
})

test_that("undirected x[v_a, v_b] ORs both orientations", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    # in canonical storage, "a-d" lives as from=a, to=d.
    # Querying x[c("d"), c("a")] should still find it via OR.
    r <- storeRead(s[c("d"), c("a")], output = "tibble")
    expect_equal(nrow(r), 1L)
})

test_that("directed x[v_a, v_b] is strict (no orientation OR)", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_directed_dt(),
                    type = "kNN", directed = TRUE)
    # Directed edges: a->b, a->c, b->c, c->a, c->d, d->a
    # x["a", "b"] should match a->b only
    r_ab <- storeRead(s["a", "b"], output = "tibble")
    expect_equal(nrow(r_ab), 1L)
    # x["b", "a"] should match nothing (no b->a edge)
    r_ba <- storeRead(s["b", "a"], output = "tibble")
    expect_equal(nrow(r_ba), 0L)
})

test_that("x[, v] is to-endpoint filter; undirected matches either side", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    # Edges incident to "a" in the tiny DT: a-b, a-c, a-d  (3 edges)
    r <- storeRead(s[, "a"], output = "tibble")
    expect_equal(nrow(r), 3L)
    # All three edges have "a" on one side
    expect_true(all("a" %in% r$from_id | "a" %in% r$to_id))
})

test_that("named-arg form is equivalent to positional", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_directed_dt(),
                    type = "kNN", directed = TRUE)
    r1 <- storeRead(s["a", "b"], output = "tibble")
    r2 <- storeRead(s[from = "a", to = "b"], output = "tibble")
    expect_equal(r1, r2)
})

test_that("negate = TRUE returns complement", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    full <- storeRead(s, output = "tibble")
    induced <- storeRead(s[c("a","b","c")], output = "tibble")
    neg     <- storeRead(s[c("a","b","c"), negate = TRUE],
                         output = "tibble")
    expect_equal(nrow(induced) + nrow(neg), nrow(full))
})


# ---- spatIDs -------------------------------------------------------------

test_that("spatIDs(parquetEdgeStore) returns full node universe when no subset", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    expect_setequal(spatIDs(s), c("a", "b", "c", "d", "e"))
})

test_that("spatIDs respects active subset (induced subgraph drops isolates)", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    # induced on {a,b,c,d} drops vertex e (only edge d-e leaves the subset)
    expect_setequal(spatIDs(s[c("a", "b", "c", "d")]),
                    c("a", "b", "c", "d"))
})

test_that("spatIDs on an empty parquetEdgeStore returns character(0)", {
    s <- parquetEdgeStore(path = tempfile())
    expect_identical(spatIDs(s), character(0L))
})


# ---- sourceAdopt ----------------------------------------------------------

test_that("sourceAdopt(parquetEdgeStore) updates both @path and @nodes@path", {
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"),
                    .tiny_undirected_dt(),
                    type = "sNN", directed = FALSE)
    gdir <- file.path(tempdir(), paste0("edge_adopt_", basename(tempfile())))
    on.exit(unlink(gdir, recursive = TRUE), add = TRUE)
    src <- gDirSource(path = gdir)
    adopted <- sourceAdopt(src, s)

    expect_true(storeExists(adopted))
    expect_true(storeExists(adopted@nodes))   # was the regression
    # Both paths must live inside the vault
    vault <- normalizePath(file.path(gdir, "artifacts"), mustWork = FALSE)
    expect_true(startsWith(normalizePath(adopted@path), vault))
    expect_true(startsWith(normalizePath(adopted@nodes@path), vault))
    # Read still works
    r <- storeRead(adopted, output = "tibble")
    expect_equal(nrow(r), 7L)
    expect_true(sourceContains(src, adopted))
})


# ---- sourceWrite (vault-resident from the start) -------------------------

test_that("sourceWrite(gDirSource, igraph) lands as a parquetEdgeStore", {
    gdir <- file.path(tempdir(), paste0("edge_srcwrite_", basename(tempfile())))
    on.exit(unlink(gdir, recursive = TRUE), add = TRUE)
    src <- gDirSource(path = gdir)
    g <- igraph::make_graph(c("a","b","b","c","c","a"), directed = FALSE)
    s <- sourceWrite(src, g)
    expect_s4_class(s, "parquetEdgeStore")
    expect_true(sourceContains(src, s))
    expect_equal(s@n_edges, 3L)
})

test_that("sourceWrite(gDirSource, data.table, store_type='parquetEdgeStore')", {
    gdir <- file.path(tempdir(), paste0("edge_srcwrite_dt_", basename(tempfile())))
    on.exit(unlink(gdir, recursive = TRUE), add = TRUE)
    src <- gDirSource(path = gdir)
    s <- sourceWrite(src, .tiny_undirected_dt(),
                     store_type = "parquetEdgeStore",
                     type = "sNN", directed = FALSE)
    expect_s4_class(s, "parquetEdgeStore")
    expect_equal(s@n_cells, 5)
})


# ---- sourceWrite(gDirSource, giotto) — gobject promotion -----------------

.tiny_gobject <- function(n_cells = 20) {
    skip_if_not_installed("GiottoClass")
    cells <- paste0("c_", sprintf("%02d", seq_len(n_cells)))
    mat <- matrix(rpois(n_cells * 50, 2), nrow = 50, ncol = n_cells,
                  dimnames = list(paste0("g_", seq_len(50)), cells))
    sl_dt <- data.table::data.table(
        cell_ID = cells, sdimx = runif(n_cells), sdimy = runif(n_cells)
    )
    sl <- GiottoClass::createSpatLocsObj(
        coordinates = sl_dt, spat_unit = "cell", provenance = "cell"
    )
    g <- GiottoClass::createGiottoObject(expression = mat)
    g <- GiottoClass::setSpatialLocations(g, sl, verbose = FALSE)

    ig <- igraph::sample_gnm(n_cells, 60, directed = FALSE)
    igraph::V(ig)$name <- cells
    nn <- methods::new("nnNetObj", network = ig, nn_type = "sNN",
        name = "sNN.test", spat_unit = "cell", feat_type = "rna",
        provenance = "cell")
    options("giotto.check_valid" = FALSE)
    GiottoClass::setNearestNetwork(g, nn, verbose = FALSE)
}

test_that("sourceWrite(gDirSource, giotto) promotes in-mem -> disk-backed", {
    skip_if_not_installed("GiottoClass")
    g <- .tiny_gobject()
    expect_null(g@source)

    gdir <- file.path(tempdir(), paste0("gw_promote_", basename(tempfile())))
    on.exit(unlink(gdir, recursive = TRUE), add = TRUE)
    src <- gDirSource(path = gdir)
    gb <- sourceWrite(src, g)

    expect_s4_class(gb@source, "gDirSource")
    nn_back <- GiottoClass::getNearestNetwork(gb, output = "nnNetObj",
        spat_unit = "cell", feat_type = "rna",
        nn_type = "sNN", name = "sNN.test")
    expect_s4_class(nn_back@network, "parquetEdgeStore")
})

test_that("sourceWrite(gDirSource, giotto) is idempotent at subobject level", {
    skip_if_not_installed("GiottoClass")
    g <- .tiny_gobject()
    gdir <- file.path(tempdir(), paste0("gw_idem_", basename(tempfile())))
    on.exit(unlink(gdir, recursive = TRUE), add = TRUE)
    src <- gDirSource(path = gdir)

    # First write: in-mem -> disk-backed.
    g1 <- sourceWrite(src, g)
    nn_class_1 <- class(GiottoClass::getNearestNetwork(g1, output = "nnNetObj",
        spat_unit = "cell", feat_type = "rna",
        nn_type = "sNN", name = "sNN.test")@network)

    # Second write: passes through setter dispatch again, but the
    # backend-aware setters' inherits(dataStore) guards skip re-writing.
    # No new manifest entries, same storage classes.
    n_before <- length(list.files(file.path(gdir, "artifacts")))
    g2 <- sourceWrite(src, g1)
    n_after <- length(list.files(file.path(gdir, "artifacts")))
    expect_equal(n_after, n_before)

    nn_class_2 <- class(GiottoClass::getNearestNetwork(g2, output = "nnNetObj",
        spat_unit = "cell", feat_type = "rna",
        nn_type = "sNN", name = "sNN.test")@network)
    expect_identical(nn_class_1, nn_class_2)

    # Different gDirSource handle pointing at same path: same behavior
    src2 <- gDirSource(path = gdir)
    n_before <- length(list.files(file.path(gdir, "artifacts")))
    g3 <- sourceWrite(src2, g1)
    n_after <- length(list.files(file.path(gdir, "artifacts")))
    expect_equal(n_after, n_before)
})

test_that("sourceWrite(gDirSource, giotto) warns on cross-source", {
    skip_if_not_installed("GiottoClass")
    g <- .tiny_gobject()
    gdir1 <- file.path(tempdir(), paste0("gw_x1_", basename(tempfile())))
    gdir2 <- file.path(tempdir(), paste0("gw_x2_", basename(tempfile())))
    on.exit(unlink(c(gdir1, gdir2), recursive = TRUE), add = TRUE)
    g_back <- sourceWrite(gDirSource(path = gdir1), g)
    expect_warning(
        sourceWrite(gDirSource(path = gdir2), g_back),
        "already backed by a different source"
    )
})
