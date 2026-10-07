# Motif enrichment on a parquetEdgeStore. Labels are named by node ID, as
# Giotto's router passes them; the method maps them to int_id through the
# node sidecar and sends only that lookup alongside the edge stream.

.motif_store <- function(n = 60L, k = 4L, seed = 3L) {
    set.seed(seed)
    ids <- sprintf("n%02d", seq_len(n))
    # a ring with k-nearest chords: connected, triangle-rich, deterministic
    from <- rep(seq_len(n), each = k)
    to <- ((from - 1L + rep(seq_len(k), n)) %% n) + 1L
    dt <- data.table::data.table(
        from = ids[pmin(from, to)], to = ids[pmax(from, to)],
        weight = 1
    )
    dt <- unique(dt)
    s <- storeWrite(storeCreate(type = "parquetEdgeStore"), dt,
        type = "sNN", directed = FALSE
    )
    lab <- stats::setNames(sample(c("A", "B", "C"), n, TRUE), ids)
    list(store = s, lab = lab)
}

.smotif_streams <- function() {
    requireNamespace("smotif", quietly = TRUE) &&
        "motif_enrichment_stream" %in% getNamespaceExports("smotif")
}

test_that("the label lookup goes through the sidecar, by name", {
    f <- .motif_store()
    lk <- .edge_label_lookup(f$store, rev(f$lab))
    nodes <- data.table::as.data.table(dplyr::collect(storeRead(f$store@nodes)))
    expect_setequal(lk$int_id, nodes$int_id)
    got <- lk$label[match(nodes$int_id, lk$int_id)]
    expect_identical(got, unname(f$lab[as.character(nodes$node_id)]))
})

test_that("the lookup leaves out unlabelled nodes and needs names", {
    f <- .motif_store()
    lk <- .edge_label_lookup(f$store, f$lab[1:10])
    expect_length(lk$int_id, 10L)
    expect_error(.edge_label_lookup(f$store, unname(f$lab)), "must be named")
})

test_that("strata and anchored_on are refused, not approximated", {
    f <- .motif_store()
    p <- Giotto::motifParam(method = "smotif", size = 3L, n_perm = 9L)
    expect_error(
        analyzeData(f$store, p, cell_type = f$lab, strata = f$lab),
        "not supported on a backed network"
    )
    expect_error(
        analyzeData(f$store, p, cell_type = f$lab, anchored_on = "n01"),
        "not supported on a backed network"
    )
})

test_that("the store gives the igraph counts and honours a pending subset", {
    skip_if_not(.smotif_streams(), "smotif has no motif_enrichment_stream()")
    f <- .motif_store()
    p <- Giotto::motifParam(method = "smotif", size = 3L, n_perm = 19L)

    for (s in list(f$store, f$store[names(f$lab)[1:25]])) {
        got <- analyzeData(s, p, cell_type = rev(f$lab))
        ig <- igraph::as.igraph(s)
        v <- igraph::V(ig)$name
        want <- analyzeData(ig, p, cell_type = unname(f$lab[v]))
        # observed counts do not depend on node order; null draws do
        j <- merge(
            data.table::as.data.table(got)[, list(motif_id, a = observed)],
            data.table::as.data.table(want)[, list(motif_id, b = observed)],
            by = "motif_id"
        )
        expect_setequal(got$motif_id, want$motif_id)
        expect_equal(j$a, j$b)
    }
})
