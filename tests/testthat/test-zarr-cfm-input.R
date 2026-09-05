# tenxZarrInput: the exprInput over the zarr cell_feature_matrix, and its
# storeRead iterator (both modes) + the sourceWrite round trip.

skip_if_no_zarr_deps()

fx <- make_zarr_fixture()

.drain <- function(itr) {
    parts <- list()
    repeat {
        b <- itr$next_batch()
        if (is.null(b)) break
        if (nrow(b)) parts[[length(parts) + 1L]] <- b
    }
    data.table::rbindlist(parts)
}

test_that("constructor is metadata-only and exact", {
    inp <- tenxZarrInput(fx$paths$cell_feature_matrix)
    expect_s4_class(inp, "tenxZarrInput")
    expect_identical(inp@n_cells, 12L)
    expect_identical(inp@n_genes, 6L) # aggregate_gene dropped
    expect_identical(inp@cell_ids, fx$truth$cell_ids)
    expect_identical(inp@feat_ids, fx$truth$cfm_feat_keys)
    # zarr snake_case classes translated to 10x display strings
    expect_identical(
        inp@feat_types,
        c(rep("Gene Expression", 4L), "Negative Control Probe",
            "Unassigned Codeword")
    )
    # ensembl ids via feature_id_col = 1
    inp_ens <- tenxZarrInput(fx$paths$cell_feature_matrix,
        feature_id_col = 1L)
    expect_identical(inp_ens@feat_ids, fx$truth$cfm_feat_ids)
})

test_that("duplicate symbols get the 10x --N disambiguation", {
    # ad-hoc archive with a repeated feature_keys entry (metadata only —
    # the constructor never touches indices/data)
    td <- withr::local_tempdir()
    .zf_write_json(td, "cell_features/.zattrs", list(
        feature_ids = c("ENSGX1", "ENSGX2", "ENSGX3"),
        feature_keys = c("DUP", "DUP", "OK"),
        feature_types = rep("gene", 3L),
        number_cells = 2L, number_features = 3L
    ))
    .zf_write_array(td, "cell_features/indptr", c(0, 1, 2, 3), "<u8",
        chunks = 4L)
    .zf_write_array(td, "cell_features/indices", c(0L, 1L, 0L), "<u4",
        chunks = 3L)
    .zf_write_array(td, "cell_features/data", c(1L, 2L, 3L), "<u4",
        chunks = 3L)
    .zf_write_array(td, "cell_features/cell_id",
        cbind(c(1, 2), c(1, 1)), "<u4", chunks = c(2L, 2L))
    inp <- tenxZarrInput(td)
    expect_identical(
        inp@feat_ids,
        GiottoDisk:::.disambiguate_feat_ids(c("DUP", "DUP", "OK"))
    )
    expect_identical(anyDuplicated(inp@feat_ids), 0L)
})

test_that("full-mode iterator obeys the protocol and matches truth", {
    inp <- tenxZarrInput(fx$paths$cell_feature_matrix, mode = "full")
    itr <- storeRead(inp)
    expect_named(itr,
        c("next_batch", "close", "cell_ids", "feat_ids",
            "n_cells", "n_genes"))
    expect_identical(itr$cell_ids(), inp@cell_ids)
    tri <- .drain(itr)
    # ascending cells; within-cell ascending features; NULL after EOF;
    # idempotent close
    expect_false(is.unsorted(tri$row_id))
    expect_true(all(tri[, diff(col_id) > 0, by = "row_id"]$V1))
    expect_null(itr$next_batch())
    expect_no_error(itr$close())
    m <- matrix(0, nrow = 12L, ncol = 6L)
    m[cbind(tri$row_id, tri$col_id)] <- tri$value
    expect_identical(m, t(unname(fx$truth$cfm)))
})

test_that("cellblock mode emits identical triplets to full mode", {
    inp_f <- tenxZarrInput(fx$paths$cell_feature_matrix, mode = "full")
    inp_b <- tenxZarrInput(fx$paths$cell_feature_matrix,
        mode = "cellblock", cells_per_block = 5L)
    tri_f <- .drain(storeRead(inp_f))
    itr_b <- storeRead(inp_b)
    parts <- list()
    n_batches <- 0L
    repeat {
        b <- itr_b$next_batch()
        if (is.null(b)) break
        n_batches <- n_batches + 1L
        parts[[n_batches]] <- b
    }
    expect_gte(n_batches, 3L) # 12 cells / 5 per block: it really windowed
    expect_identical(tri_f, data.table::rbindlist(parts))
})

test_that("aggregate_gene can be kept when asked", {
    inp <- tenxZarrInput(fx$paths$cell_feature_matrix,
        drop_aggregate = FALSE)
    expect_identical(inp@n_genes, 7L)
    expect_identical(inp@feat_types[7L], "Aggregate Gene")
    tri <- .drain(storeRead(inp))
    agg <- tri[tri$col_id == 7L, ]
    truth_tot <- colSums(fx$truth$cfm)
    expect_identical(agg$value, unname(truth_tot[truth_tot > 0]))
})

test_that("sourceWrite streams zarr triplets into a parquetExprStore", {
    gsrc <- gDirSource(path = withr::local_tempdir())
    inp <- tenxZarrInput(fx$paths$cell_feature_matrix)
    pe <- sourceWrite(gsrc, inp, store_type = "parquetExpr",
        verbose = FALSE)
    expect_s4_class(pe, "parquetExprStore")
    expect_identical(pe@cell_ids, fx$truth$cell_ids)
    expect_identical(dim(pe), c(6, 12))
    tri <- dplyr::collect(storeRead(pe))
    m <- matrix(0, nrow = 12L, ncol = 6L)
    m[cbind(tri$row_id, tri$col_id)] <- tri$value
    expect_identical(m, t(unname(fx$truth$cfm)))
})
