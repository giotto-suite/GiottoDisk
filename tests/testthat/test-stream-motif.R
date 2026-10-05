# Motif enrichment on a parquetEdgeStore. Labels are named by node ID, as
# Giotto's router passes them, and must be realigned by name.

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

test_that("the active node set matches the igraph the store reads into", {
    f <- .motif_store()
    for (s in list(f$store, f$store[names(f$lab)[1:25]])) {
        nodes <- .edge_active_nodes(s)
        expect_setequal(nodes$node_id, igraph::V(igraph::as.igraph(s))$name)
        expect_false(is.unsorted(nodes$int_id))
    }
})

test_that("labels must be named, and every node needs one", {
    f <- .motif_store()
    expect_error(.motif_realign(unname(f$lab), names(f$lab), "cell_type"),
        "must be named")
    expect_error(.motif_realign(f$lab[-1], names(f$lab), "cell_type"),
        "1 network node")
    expect_identical(.motif_realign(rev(f$lab), names(f$lab), "cell_type"),
        unname(f$lab))
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
