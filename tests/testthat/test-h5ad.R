# AnnData .h5ad inputs (R/h5ad.R), against the synthetic file of
# helper-h5ad-fixture.R. Every read is compared with the source R objects,
# never with a second read through the same path.

skip_if_not_installed("hdf5r")

proj <- withr::local_tempdir()
gsrc <- suppressMessages(gDirSource(proj))


# expression ####

test_that("CSR X streams into a cell-major store equal to the source", {
    fx <- make_h5ad(format = "csr")
    inp <- h5adInput(fx$path, batch_records = 40)
    expect_identical(inp@cell_ids, fx$cells)
    expect_identical(inp@feat_ids, .disambiguate_feat_ids(fx$genes))
    pe <- sourceWrite(gsrc, inp, store_type = "parquetExpr")
    expect_equal(dim(pe), c(ncol(fx$M), nrow(fx$M)))
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
    expect_true(.parts_cell_major(pe))
})

test_that("CSC X is reordered to the same cell-major store", {
    fx <- make_h5ad(format = "csc", seed = 2L)
    pe <- sourceWrite(gsrc, h5adInput(fx$path, batch_records = 40),
                      store_type = "parquetExpr")
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
    expect_true(.parts_cell_major(pe))
})

test_that("the session option pins the batch size; the argument wins", {
    fx <- make_h5ad(seed = 3L)
    withr::local_options(giottodisk.h5_batch_rows = 40)
    pe <- sourceWrite(gsrc, h5adInput(fx$path), store_type = "parquetExpr")
    expect_true(.parts_cell_major(pe))                      # several parts
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
    expect_identical(.h5ad_max_rows(NULL), 40)
    expect_identical(.h5ad_max_rows(1e3), 1e3)
})

test_that("layers and raw resolve; a dense layer streams its non-zeros", {
    fx <- make_h5ad()
    pe <- sourceWrite(gsrc, h5adInput(fx$path, layer = "counts"),
                      store_type = "parquetExpr")
    expect_equal(.store_triplets(pe), .source_triplets(fx$M * 2))
    raw <- h5adInput(fx$path, layer = "raw")
    expect_identical(raw@feat_ids, fx$genes[1:10])
    expect_identical(raw@n_genes, 10L)
    dn <- sourceWrite(gsrc, h5adInput(fx$path, layer = "dense", batch_records = 100),
                      store_type = "parquetExpr")
    expect_true(.parts_cell_major(dn))                      # several blocks
    expect_equal(.store_triplets(dn), .source_triplets(fx$M))
    expect_error(h5adInput(fx$path, layer = "nope"), "Available: X, raw, counts, dense")
})


test_that("a dense X streams, and the sparse layers are named", {
    fx <- make_h5ad(seed = 6L)
    h5 <- hdf5r::H5File$new(fx$path, mode = "r+")
    h5$link_delete("X")
    .w_array(h5, "X", as.matrix(fx$M))
    h5$close_all()
    expect_message(inp <- h5adInput(fx$path, batch_records = 100), "'counts'")
    pe <- sourceWrite(gsrc, inp, store_type = "parquetExpr")
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
})


# networks ####

test_that("the kNN streams without self-loops, weighted by connectivity", {
    fx <- make_h5ad()
    inp <- h5adEdgeInput(fx$path, connectivities_key = "connectivities",
        distances_key = "distances", batch_records = 20)
    st <- sourceWrite(gsrc, inp, store_type = "parquetEdge", type = "kNN")
    expect_true(st@directed)
    expect_identical(st@n_edges, as.numeric(length(fx$nb)))
    expect_gt(length(list.files(file.path(st@path, "edges"))), 1L)

    got <- data.table::as.data.table(storeRead(st, output = "tibble"))
    data.table::setnames(got, c("from_id", "to_id"), c("from", "to"))
    from <- rep(seq_along(fx$cells), each = ncol(fx$nb)); to <- as.vector(t(fx$nb))
    want <- data.table::data.table(from = fx$cells[from], to = fx$cells[to],
        weight = fx$conn[cbind(from, to)], distance = fx$D[cbind(from, to)])
    data.table::setorder(got, from, to); data.table::setorder(want, from, to)
    expect_equal(got[, .(from, to, weight, distance)], want)
    expect_true(all(got[, .N, by = from]$N == ncol(fx$nb)))  # uniform degree

    # same edges through the in-memory writer give the same store contents
    ref <- sourceWrite(gsrc, edgeDTInput(want, weight_col = "weight",
        distance_col = "distance"), store_type = "parquetEdge", type = "kNN")
    ref_dt <- data.table::as.data.table(storeRead(ref, output = "tibble"))
    data.table::setnames(ref_dt, c("from_id", "to_id"), c("from", "to"))
    data.table::setorder(ref_dt, from, to)
    expect_equal(got[, .(from, to, weight, distance)],
                 ref_dt[, .(from, to, weight, distance)])
})

test_that("a connectivities-only graph is undirected, upper triangle", {
    fx <- make_h5ad()
    inp <- h5adEdgeInput(fx$path, connectivities_key = "spatial_connectivities")
    st <- sourceWrite(gsrc, inp, store_type = "parquetEdge", type = "spatial",
        directed = inp@params$directed)
    expect_false(st@directed)
    up <- methods::as(Matrix::triu(fx$conn, 1), "TsparseMatrix")
    expect_identical(st@n_edges, as.numeric(length(up@x)))
    e <- data.table::setDT(dplyr::collect(arrow::open_dataset(file.path(st@path, "edges"))))
    expect_true(all(e$from_id < e$to_id))
    expect_equal(sort(e$weight), sort(up@x))
})
