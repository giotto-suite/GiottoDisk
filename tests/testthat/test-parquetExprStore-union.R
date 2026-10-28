# unionParquetExprStore: cbind2 of parquetExprStores, [-subset across
# substores, the duckdb carrier on a union, storeWrite / saveRDS round-trips.

# ---- cbind / unionParquetExprStore --------------------------------------
# These tests lock the cbind semantics. They reflect the *current behavior*
# rather than an ideal target -- the schema and key-uniqueness questions
# are unresolved (see comment on each test for what the test is asserting
# vs what an ideal contract would be).

.tiny_distinct_mat <- function(n_genes = 3L, cell_prefix = "c1_",
    cell_offset = 0L, value_offset = 0L, seed = 1L) {
    n_cells <- 4L
    set.seed(seed)
    m <- Matrix::sparseMatrix(
        i = c(1L, 2L, 3L, 1L),
        j = c(1L, 2L, 3L, 4L),
        x = as.double(c(10L, 20L, 30L, 40L) + value_offset),
        dims = c(n_genes, n_cells)
    )
    rownames(m) <- paste0("gene", seq_len(n_genes))
    colnames(m) <- paste0(cell_prefix, seq_len(n_cells) + cell_offset)
    m
}

test_that("cbind2(parquetExprStore, parquetExprStore) returns unionParquetExprStore", {
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "a", seed = 1))
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "b", cell_offset = 4L,
            value_offset = 100L, seed = 2))
    u <- cbind2(pe1, pe2)
    expect_s4_class(u, "unionParquetExprStore")
    expect_equal(u@n_genes, 3)
    expect_equal(u@n_cells, 8)
    expect_equal(u@feat_ids, c("gene1", "gene2", "gene3"))
    expect_equal(u@cell_ids, c("a1", "a2", "a3", "a4", "b5", "b6", "b7", "b8"))
})

test_that("union storeRead contains all rows from each substore", {
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "a", seed = 1))
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "b", cell_offset = 4L,
            value_offset = 100L, seed = 2))
    u <- cbind2(pe1, pe2)
    df <- as.data.frame(dplyr::collect(storeRead(u)))
    # Each substore contributes 4 non-zero entries; total = 8.
    expect_equal(nrow(df), 8L)
    # Values from both substores are present.
    expect_setequal(df$value, c(10, 20, 30, 40, 110, 120, 130, 140))
})

test_that("union storeRead carries source_id for substore disambiguation", {
    # NB: this locks the *current* behavior -- source_id is in the schema
    # because (row_id, col_id) is only locally unique per substore, so the
    # composite (source_id, row_id, col_id) is what uniquely identifies a
    # (cell, gene) entry post-union. If parquetExprStore is later changed
    # to re-number row_ids globally at union time, this test should be
    # updated together with the schema contract.
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "a", seed = 1))
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "b", cell_offset = 4L,
            value_offset = 100L, seed = 2))
    u <- cbind2(pe1, pe2)
    df <- as.data.frame(dplyr::collect(storeRead(u)))
    expect_true("source_id" %in% names(df))
    expect_setequal(unique(df$source_id), c(pe1@uid, pe2@uid))
    # Composite key (source_id, row_id, col_id) is globally unique.
    composite <- paste(df$source_id, df$row_id, df$col_id)
    expect_equal(length(unique(composite)), nrow(df))
})

test_that("cbind2 is associative across multi-arity calls", {
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "a", seed = 1))
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "b", cell_offset = 4L, seed = 2))
    pe3 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "c", cell_offset = 8L, seed = 3))
    # cbind2(cbind2(pe1, pe2), pe3) is the left-fold form base R uses
    u <- cbind2(cbind2(pe1, pe2), pe3)
    expect_s4_class(u, "unionParquetExprStore")
    expect_equal(u@n_cells, 12)
    expect_equal(length(u@stores), 3L)
    # all uids represented in the read
    df <- as.data.frame(dplyr::collect(storeRead(u)))
    expect_setequal(unique(df$source_id), c(pe1@uid, pe2@uid, pe3@uid))
})


# `[` subset + storeRead on a unionParquetExprStore ####
# Previously broken: the union wrapped each substore's read_fun with
# `dplyr::filter` for per-substore cell_idx/gene_idx, returning
# arrow_dplyr_query -- which `arrow::open_dataset(list(...))` can't
# accept. Fix: open each substore as a Dataset cleanly, then apply a
# composite source_id-aware filter expression once at the union level.

test_that("union [-subset on a single substore's cells reads correctly", {
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "a", seed = 1))
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "b", cell_offset = 4L, seed = 2))
    u <- cbind2(pe1, pe2)

    sub <- u[, "a1"]
    expect_s4_class(sub, "unionParquetExprStore")
    df <- as.data.frame(dplyr::collect(storeRead(sub)))
    # Only rows from pe1, only with row_id == 1 (the first cell)
    expect_setequal(unique(df$source_id), pe1@uid)
    expect_true(all(df$row_id == 1L))
})

test_that("union [-subset across substores reads correctly", {
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "a", seed = 1))
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_distinct_mat(cell_prefix = "b", cell_offset = 4L, seed = 2))
    u <- cbind2(pe1, pe2)

    # One cell from each substore
    sub <- u[, c("a2", "b6")]
    df <- as.data.frame(dplyr::collect(storeRead(sub)))
    expect_setequal(unique(df$source_id), c(pe1@uid, pe2@uid))
    # pe1 row_id should be 2 (a2 is its 2nd cell), pe2 row_id should be 2 (b6 is its 2nd cell)
    pe1_rows <- df[df$source_id == pe1@uid, ]
    pe2_rows <- df[df$source_id == pe2@uid, ]
    expect_true(all(pe1_rows$row_id == 2L))
    expect_true(all(pe2_rows$row_id == 2L))
})

test_that("union [-subset round-trips to a matrix matching in-mem cbind subset", {
    m1 <- Matrix::sparseMatrix(
        i = c(1, 2, 3), j = c(1, 2, 1), x = c(11, 12, 13),
        dims = c(3, 2),
        dimnames = list(c("g1", "g2", "g3"), c("c1", "c2"))
    )
    m2 <- Matrix::sparseMatrix(
        i = c(1, 3), j = c(1, 2), x = c(21, 22),
        dims = c(3, 2),
        dimnames = list(c("g1", "g2", "g3"), c("c3", "c4"))
    )
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m1)
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m2)
    u <- cbind2(pe1, pe2)

    reconstruct <- function(store) {
        df <- as.data.frame(dplyr::collect(storeRead(store)))
        uid_to_idx <- stats::setNames(
            seq_along(store@stores),
            vapply(store@stores, function(s) s@uid, character(1L))
        )
        offsets <- c(0L, cumsum(vapply(store@stores,
            function(s) s@n_cells, numeric(1L))))
        df$global_cell_idx <- df$row_id +
            offsets[uid_to_idx[df$source_id]]
        Matrix::sparseMatrix(
            i = df$col_id, j = df$global_cell_idx, x = df$value,
            dims = c(store@n_genes, store@n_cells),
            dimnames = list(store@feat_ids, store@cell_ids)
        )
    }

    ref_full <- as.matrix(cbind(m1, m2))
    expect_equal(as.matrix(reconstruct(u)), ref_full)

    ref_cross <- as.matrix(cbind(m1, m2)[, c("c1", "c3")])
    expect_equal(as.matrix(reconstruct(u[, c("c1", "c3")])), ref_cross)

    ref_sub2 <- as.matrix(cbind(m1, m2)[, c("c3", "c4")])
    expect_equal(as.matrix(reconstruct(u[, c("c3", "c4")])), ref_sub2)
})


# ---- union subset gaps -------------------------------------------------

.cbind_pair <- function() {
    m1 <- Matrix::sparseMatrix(
        i = c(1, 2, 3), j = c(1, 2, 1), x = c(11, 12, 13),
        dims = c(3, 2),
        dimnames = list(c("g1", "g2", "g3"), c("c1", "c2"))
    )
    m2 <- Matrix::sparseMatrix(
        i = c(1, 3), j = c(1, 2), x = c(21, 22),
        dims = c(3, 2),
        dimnames = list(c("g1", "g2", "g3"), c("c3", "c4"))
    )
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m1)
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m2)
    list(m1 = m1, m2 = m2, pe1 = pe1, pe2 = pe2, u = cbind2(pe1, pe2))
}

# Reconstruct a sparse Matrix from union storeRead output, resolving
# (source_id, row_id_local) -> global cell position via per-substore offsets.
.reconstruct_from_union <- function(u) {
    df <- as.data.frame(dplyr::collect(storeRead(u)))
    uids <- vapply(u@stores, function(s) s@uid, character(1L))
    offsets <- c(0L, cumsum(vapply(u@stores, function(s) s@n_cells,
        numeric(1L))))
    uid_to_idx <- stats::setNames(seq_along(uids), uids)
    df$global_row <- df$row_id + offsets[uid_to_idx[df$source_id]]
    Matrix::sparseMatrix(
        i = df$col_id, j = df$global_row, x = df$value,
        dims = c(u@n_genes, u@n_cells),
        dimnames = list(u@feat_ids, u@cell_ids)
    )
}


test_that("union [-subset: gene-only subset reads correctly", {
    s <- .cbind_pair()
    sub <- s$u["g1", ]
    expect_equal(sub@n_genes, 1L)
    expect_equal(sub@n_cells, 4L)
    expect_equal(sub@feat_ids, "g1")
    df <- as.data.frame(dplyr::collect(storeRead(sub)))
    # g1 corresponds to col_id 1 in originals; after subset that's still
    # the only col_id present
    expect_true(all(df$col_id == 1L))
    # cells c1 (val 11, src pe1), c3 (val 21, src pe2) have g1
    expect_setequal(df$value, c(11, 21))
})

test_that("union [-subset: combined gene + cell subset reads correctly", {
    s <- .cbind_pair()
    sub <- s$u["g1", c("c1", "c3")]
    expect_equal(sub@n_genes, 1L)
    expect_equal(sub@n_cells, 2L)
    df <- as.data.frame(dplyr::collect(storeRead(sub)))
    expect_equal(nrow(df), 2L)
    expect_setequal(df$value, c(11, 21))
})

test_that("union [-subset: cell selection that zeroes a substore drops it", {
    s <- .cbind_pair()
    sub <- s$u[, c("c1", "c2")]  # only pe1's cells
    expect_equal(length(sub@stores), 1L)
    expect_equal(sub@stores[[1L]]@uid, s$pe1@uid)
    expect_equal(sub@n_cells, 2L)
})

test_that("duckdb carrier: union read matches Acero, subset and not", {
    # The union's per-substore subset state is one composite source_id-aware
    # expression, used unchanged on both carriers.
    skip_if_no_duckdb()
    s <- .cbind_pair()
    .expect_dd_parity(s$u)
    .expect_dd_parity(s$u[, c("c1", "c3", "c4")])
})

test_that("duckdb carrier: a union-level multiply op matches Acero", {
    # Payload is keyed by substore uid, so this exercises the composite
    # (source_id, row_id) join across both arms of the UNION ALL.
    skip_if_no_duckdb()
    s <- .cbind_pair()
    set.seed(21)
    s$u@ops <- list(list(type = "multiply", axis = "cell",
        factors = stats::setNames(
            list(runif(as.integer(s$pe1@n_cells)),
                 runif(as.integer(s$pe2@n_cells))),
            c(s$pe1@uid, s$pe2@uid))))
    .expect_dd_parity(s$u)
})


# ---- storeWrite from a (possibly subset) union -----------------------------

test_that("storeWrite from full union round-trips to a new flat store", {
    s <- .cbind_pair()
    out <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), s$u)
    expect_equal(out@n_cells, 4L)
    expect_equal(out@cell_ids, c("c1", "c2", "c3", "c4"))
    df <- as.data.frame(dplyr::collect(storeRead(out)))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(out@n_genes, out@n_cells),
        dimnames = list(out@feat_ids, out@cell_ids)
    )
    expect_equal(as.matrix(rt), as.matrix(cbind(s$m1, s$m2)))
})

test_that("storeWrite from union subset (cross-substore) round-trips correctly", {
    s <- .cbind_pair()
    u_sub <- s$u[, c("c1", "c3")]
    out <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), u_sub)
    df <- as.data.frame(dplyr::collect(storeRead(out)))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(out@n_genes, out@n_cells),
        dimnames = list(out@feat_ids, out@cell_ids)
    )
    expect_equal(as.matrix(rt), as.matrix(cbind(s$m1, s$m2)[, c("c1", "c3")]))
})

test_that("storeWrite from union with combined gene+cell subset", {
    s <- .cbind_pair()
    u_sub <- s$u[c("g1", "g3"), c("c1", "c3", "c4")]
    out <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), u_sub)
    df <- as.data.frame(dplyr::collect(storeRead(out)))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(out@n_genes, out@n_cells),
        dimnames = list(out@feat_ids, out@cell_ids)
    )
    expect_equal(as.matrix(rt),
        as.matrix(cbind(s$m1, s$m2)[c("g1", "g3"), c("c1", "c3", "c4")]))
})


# ---- saveRDS roundtrip of a union store ----------------------------------

test_that("saveRDS roundtrip of unionParquetExprStore preserves slots + read", {
    s <- .cbind_pair()
    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(s$u, tmp)
    rt <- readRDS(tmp)
    expect_s4_class(rt, "unionParquetExprStore")
    expect_equal(rt@n_cells, s$u@n_cells)
    expect_equal(rt@cell_ids, s$u@cell_ids)
    expect_equal(rt@feat_ids, s$u@feat_ids)
    # storeRead post-roundtrip yields same matrix
    expect_equal(
        as.matrix(.reconstruct_from_union(rt)),
        as.matrix(cbind(s$m1, s$m2))
    )
})


# ---- error path: duplicate cell_ids across substores --------------------

test_that("union construction rejects substores with duplicate cell_ids", {
    mat <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(1, 2),
        dimnames = list("g1", c("c1", "c2")))
    pe1 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    pe2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    # Both have cell_ids c("c1", "c2") -- union must reject
    expect_error(cbind2(pe1, pe2), "duplicate cell_ids")
})
