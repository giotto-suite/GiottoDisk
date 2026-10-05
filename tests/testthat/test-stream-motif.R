# Motif enrichment on a parquetEdgeStore: "auto" resolution, the igraph
# fallback, and the stream path. Labels are named by node ID, as Giotto's
# router passes them, and every path must realign by name.

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
    shuffled <- rev(f$lab)
    expect_identical(.motif_realign(shuffled, names(f$lab), "cell_type"),
        unname(f$lab))
})

test_that("smotifrsParam is a motifParam built through Giotto's factory", {
    p <- smotifrsParam(size = 4L, null = "conditional", n_perm = 9L)
    expect_s4_class(p, "smotifrsParam")
    expect_true(methods::is(p, "motifParam"))
    expect_identical(p$size, 4L)
    expect_identical(p$null, "conditional")
    expect_error(smotifrsParam(size = 5L))
})

test_that("unsupported options fall back to the igraph engine", {
    skip_if_not_installed("smotif")
    f <- .motif_store()
    strata <- stats::setNames(rep(c("r1", "r2"), length.out = length(f$lab)),
        names(f$lab))
    p <- smotifrsParam(size = 3L, null = "stratified", n_perm = 19L)

    got <- analyzeData(f$store, p, cell_type = rev(f$lab), strata = strata)

    ig <- igraph::as.igraph(f$store)
    v <- igraph::V(ig)$name
    want <- analyzeData(ig, methods::new("smotifParam", param = p@param),
        cell_type = unname(f$lab[v]), strata = unname(strata[v])
    )
    expect_identical(got$motif_id, want$motif_id)
    expect_equal(got$observed, want$observed)
    expect_equal(got$p_enrich, want$p_enrich)
})

test_that("auto on a store gives the igraph answer, by either route", {
    skip_if_not_installed("smotif")
    f <- .motif_store()
    p <- Giotto::motifParam(size = 3L, n_perm = 19L)

    got <- analyzeData(f$store, p, cell_type = rev(f$lab))

    ig <- igraph::as.igraph(f$store)
    v <- igraph::V(ig)$name
    want <- analyzeData(ig, p, cell_type = unname(f$lab[v]))
    # observed counts are engine-independent; the null draws are not, since
    # the stream permutes over int_id order rather than vertex order
    j <- merge(
        data.table::as.data.table(got)[, list(motif_id, a = observed)],
        data.table::as.data.table(want)[, list(motif_id, b = observed)],
        by = "motif_id"
    )
    expect_setequal(got$motif_id, want$motif_id)
    expect_equal(j$a, j$b)
})

test_that("the stream path honours a pending subset", {
    skip_if_not(.smotifrs_streams(), "smotifrs has no motif_enrichment_stream()")
    f <- .motif_store()
    sub <- f$store[names(f$lab)[1:25]]
    expect_gt(length(sub@ops), 0L)
    p <- smotifrsParam(size = 3L, n_perm = 19L)

    got <- analyzeData(sub, p, cell_type = rev(f$lab))
    ig <- igraph::as.igraph(sub)
    v <- igraph::V(ig)$name
    want <- analyzeData(ig, methods::new("smotifParam", param = p@param),
        cell_type = unname(f$lab[v])
    )
    j <- merge(
        data.table::as.data.table(got)[, list(motif_id, a = observed)],
        data.table::as.data.table(want)[, list(motif_id, b = observed)],
        by = "motif_id"
    )
    expect_setequal(got$motif_id, want$motif_id)
    expect_equal(j$a, j$b)
    expect_lt(sum(got$observed),
        sum(analyzeData(f$store, p, cell_type = f$lab)$observed))
})
