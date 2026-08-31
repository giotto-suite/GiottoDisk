# Tests for parquetExprStore: storage class for streaming-friendly
# long-format expression matrices.

# helper: build a tiny labeled dgCMatrix
.tiny_mat <- function(n_genes = 5L, n_cells = 8L, density = 0.4, seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("gene", seq_len(n_genes))
    colnames(m) <- paste0("cell", seq_len(n_cells))
    m
}

# helper: write a tiny MatrixMarket triple to a directory in 10x layout
.write_mtx_triple <- function(mat, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(colnames(mat), file.path(out_dir, "barcodes.tsv"))
    feats <- data.frame(
        ensg   = paste0("ENSG", sprintf("%011d", seq_len(nrow(mat)))),
        symbol = rownames(mat),
        type   = "Gene Expression"
    )
    write.table(feats, file.path(out_dir, "features.tsv"),
        sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    Matrix::writeMM(mat, file.path(out_dir, "matrix.mtx"))
    out_dir
}


# ---- Class basics ----------------------------------------------------------

test_that("parquetExprStore constructor builds a valid object", {
    pe <- parquetExprStore(
        path     = tempfile(fileext = ".parquet"),
        cell_ids = c("c1", "c2", "c3"),
        feat_ids = c("g1", "g2")
    )
    expect_s4_class(pe, "parquetExprStore")
    expect_s4_class(pe, "fileStore")
    expect_equal(pe@n_cells, 3)
    expect_equal(pe@n_genes, 2)
    expect_equal(nrow(pe), 2)        # genes
    expect_equal(ncol(pe), 3)        # cells
    expect_equal(dim(pe), c(2, 3))
    expect_equal(rownames(pe), c("g1", "g2"))
    expect_equal(colnames(pe), c("c1", "c2", "c3"))
})

test_that("constructor validates id length vs n_cells / n_genes", {
    expect_error(
        parquetExprStore(cell_ids = c("c1", "c2"), n_cells = 5),
        "length\\(cell_ids\\)"
    )
    expect_error(
        parquetExprStore(feat_ids = c("g1"), n_genes = 4),
        "length\\(feat_ids\\)"
    )
})


# ---- storeWrite (memoryMatrix) round-trip ---------------------------------

test_that("storeWrite from dgCMatrix is lossless", {
    mat <- .tiny_mat()
    pe  <- parquetExprStore(path = tempfile(fileext = ".parquet"))
    pe  <- storeWrite(pe, mat)

    expect_equal(pe@n_genes, nrow(mat))
    expect_equal(pe@n_cells, ncol(mat))
    expect_equal(pe@feat_ids, rownames(mat))
    expect_equal(pe@cell_ids, colnames(mat))

    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    # source_id is part of the schema -- it's the substore selector that
    # makes (source_id, row_id) -> cell_ID lookup work across unions, where
    # each substore has its own @cell_ids namespace. See the cbind / union
    # tests below.
    expect_setequal(names(df), c("row_id", "col_id", "value", "source_id"))

    # Reconstruct and compare element-wise
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    rownames(rt) <- pe@feat_ids
    colnames(rt) <- pe@cell_ids
    expect_equal(as.matrix(mat), as.matrix(rt))
})

test_that("Parquet payload has row_id, col_id, value, source_id (sorted by row_id)", {
    mat <- .tiny_mat(seed = 2)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    df  <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_equal(sort(names(df)),
        c("col_id", "row_id", "source_id", "value"))
    expect_true(!is.unsorted(df$row_id))
})


# ---- storeRead returns a lazy Arrow Dataset ------------------------------

test_that("storeRead returns an Arrow Dataset (lazy, not materialized)", {
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    ds  <- storeRead(pe)
    expect_true(inherits(ds, "Dataset"))
    # Can be filtered without collecting
    sub <- ds |> dplyr::filter(row_id == 1L) |> dplyr::collect()
    expect_true(all(sub$row_id == 1L))
})


# ---- storeRead on the duckdb carrier -------------------------------------
# Same store, same modifications, different engine executing the scan. See
# helper-pestore-parity.R for why these sit beside the Acero tests.

test_that("storeRead(output = 'duckdb') returns a lazy tbl_dbi", {
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat())
    expect_s3_class(storeRead(pe, output = "duckdb"), "tbl_dbi")
})

test_that("storeRead(output = 'duckdb') creates a connection when given none", {
    # The Arrow bridge this replaces hard-errored without duckdb_params$conn.
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat())
    expect_no_error(dplyr::collect(storeRead(pe, output = "duckdb")))
})

test_that("storeRead(output = 'duckdb') uses a supplied connection as given", {
    skip_if_no_duckdb()
    conn <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat())
    tbl <- storeRead(pe, output = "duckdb", duckdb_params = list(conn = conn))
    expect_identical(dbplyr::remote_con(tbl), conn)
})

test_that("storeRead(output = 'duckdb') errors when @uid is not the on-disk partition", {
    # Minting a store from a path gives it a fresh uid, which points the scan
    # at a source_id= directory that is not there. Must not scan empty.
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat())
    pe@uid <- "definitely_not_on_disk"
    expect_error(storeRead(pe, output = "duckdb"), "must match the on-disk")
})

test_that("storeRead(output = 'duckdb') honours duckdb_params$name", {
    # Documented, and honoured by the arrow bridge this path replaced -- it
    # named the registered scanner, which was likewise what the returned tbl
    # pointed at.
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat())
    tbl <- storeRead(pe, output = "duckdb",
        duckdb_params = list(name = "my_pe_view"))
    expect_equal(dbplyr::remote_name(tbl), "my_pe_view")
    expect_gt(nrow(dplyr::collect(tbl)), 0L)
})

test_that("storeRead(output = 'duckdb') rejects conn / name passed directly", {
    # They land in `...`, which no duckdb path reads. Silently ignoring them
    # hands back a working tbl_dbi on a connection the caller did not choose.
    skip_if_no_duckdb()
    conn <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat())
    expect_error(storeRead(pe, output = "duckdb", conn = conn),
        "must be passed inside `duckdb_params`")
    expect_error(storeRead(pe, output = "duckdb", name = "x"),
        "must be passed inside `duckdb_params`")
    # ... and only on the duckdb path -- `...` is legitimate elsewhere.
    expect_no_error(storeRead(pe, output = "query"))
})

test_that("storeRead(output = 'duckdb') honours fields", {
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat())
    df <- dplyr::collect(storeRead(pe, output = "duckdb",
        fields = c("row_id", "value")))
    expect_setequal(names(df), c("row_id", "value"))
    expect_gt(nrow(df), 0L)
})

test_that("duckdb carrier: unmodified read matches Acero", {
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat(n_genes = 12L, n_cells = 20L, seed = 21))
    .expect_dd_parity(pe)
})

test_that("duckdb carrier: values match the source matrix, not just Acero", {
    # Parity alone cannot catch a bug both carriers share, since they run the
    # same dplyr. Anchor one case against the matrix the store was written from.
    skip_if_no_duckdb()
    mat <- .tiny_mat(n_genes = 6L, n_cells = 8L, density = 0.6, seed = 7)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    d <- data.table::as.data.table(
        dplyr::collect(storeRead(pe, output = "duckdb")))
    expect_equal(nrow(d), length(mat@x))

    # row_id is the cell axis (matrix column), col_id the gene axis (row).
    trip <- Matrix::summary(mat)
    data.table::setorderv(d, c("row_id", "col_id"))
    trip <- trip[order(trip$j, trip$i), ]
    expect_equal(d$row_id, as.integer(trip$j))
    expect_equal(d$col_id, as.integer(trip$i))
    expect_equal(d$value, trip$x)
})


# ---- @ops chain on the duckdb carrier ------------------------------------
# One dplyr executor serves both engines, so these assert the chain lowers
# rather than that duckdb has its own translation.

test_that("duckdb carrier: a log op matches Acero", {
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat(n_genes = 12L, n_cells = 20L, seed = 22))
    pe@ops <- list(list(type = "log", base = 2))
    .expect_dd_parity(pe)
})

test_that("duckdb carrier: multiply matches Acero on every axis", {
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat(n_genes = 12L, n_cells = 20L, seed = 23))

    pe_all <- pe
    pe_all@ops <- list(list(type = "multiply", axis = "all", factors = 2.5))
    .expect_dd_parity(pe_all)

    set.seed(11)
    pe_cell <- pe
    pe_cell@ops <- list(list(type = "multiply", axis = "cell",
        factors = stats::setNames(list(runif(as.integer(pe@n_cells))), pe@uid)))
    .expect_dd_parity(pe_cell)

    pe_feat <- pe
    pe_feat@ops <- list(list(type = "multiply", axis = "feat",
        factors = stats::setNames(list(runif(as.integer(pe@n_genes))), pe@uid)))
    .expect_dd_parity(pe_feat)
})

test_that("duckdb carrier: op order is preserved, not flattened", {
    # log-then-multiply and multiply-then-log are different values. Each must
    # match its own Acero twin -- a carrier that reordered or collapsed the
    # chain would still agree with itself but not with Acero.
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat(n_genes = 12L, n_cells = 20L, seed = 24))
    lg <- list(type = "log", base = 2)
    ml <- list(type = "multiply", axis = "all", factors = 3)

    fwd <- pe; fwd@ops <- list(lg, ml)
    rev <- pe; rev@ops <- list(ml, lg)

    a <- .expect_dd_parity(fwd)
    b <- .expect_dd_parity(rev)
    expect_false(isTRUE(all.equal(a$value, b$value)))
})

test_that("duckdb carrier: a subset and an op chain combine as on Acero", {
    skip_if_no_duckdb()
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat(n_genes = 12L, n_cells = 20L, seed = 25))
    sub <- pe[c(2L, 6L, 10L), 4:16]
    sub@ops <- list(list(type = "log", base = 2))
    .expect_dd_parity(sub)
})


# ---- Streaming mtxInput + storeWrite -------------------------------------

test_that("mtxInput + storeWrite on a synthetic 10x triple is lossless", {
    mat     <- .tiny_mat(n_genes = 6, n_cells = 10, density = 0.5, seed = 3)
    src_dir <- .write_mtx_triple(mat, file.path(tempdir(), "tiny_10x"))
    out     <- file.path(tempdir(), "tiny_10x_out")

    inp <- mtxInput(mtx_path = file.path(src_dir, "matrix.mtx"),
                    feature_id_col = 2L)
    pe  <- storeWrite(parquetExprStore(path = out), inp)

    expect_s4_class(pe, "parquetExprStore")
    expect_equal(pe@n_genes, nrow(mat))
    expect_equal(pe@n_cells, ncol(mat))
    expect_equal(pe@cell_ids, colnames(mat))
    expect_equal(pe@feat_ids, rownames(mat))

    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_equal(nrow(df), Matrix::nnzero(mat))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    expect_equal(unname(as.matrix(mat)), unname(as.matrix(rt)))

    # Cleanup
    unlink(src_dir, recursive = TRUE)
    unlink(out, recursive = TRUE)
})

test_that("mtxInput + storeWrite produces multiple chunks when batch_lines is small", {
    mat     <- .tiny_mat(n_genes = 4, n_cells = 20, density = 0.6, seed = 4)
    src_dir <- .write_mtx_triple(mat, file.path(tempdir(), "tiny_10x_b"))
    out     <- file.path(tempdir(), "tiny_10x_b_out")

    # Force multi-batch by setting a tiny batch_lines
    inp <- mtxInput(mtx_path = file.path(src_dir, "matrix.mtx"),
                    feature_id_col = 2L,
                    batch_lines    = 5L)
    pe  <- storeWrite(parquetExprStore(path = out), inp)

    expect_true(dir.exists(out))
    # storeWrite lays chunks under <out>/source_id=<uid>/part-*.parquet
    chunks <- list.files(out, pattern = "\\.parquet$",
                         full.names = TRUE, recursive = TRUE)
    expect_gt(length(chunks), 1L)

    # Aggregate read still works (Arrow handles hive partitioning)
    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_equal(nrow(df), Matrix::nnzero(mat))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    expect_equal(unname(as.matrix(mat)), unname(as.matrix(rt)))

    unlink(src_dir, recursive = TRUE)
    unlink(out,     recursive = TRUE)
})

test_that("mtxInput + storeWrite rejects mismatched barcode / mtx header counts", {
    mat     <- .tiny_mat(n_genes = 4, n_cells = 6, seed = 5)
    src_dir <- .write_mtx_triple(mat, file.path(tempdir(), "tiny_10x_mm"))
    # corrupt barcodes file: too few lines (mtx header still says 6 cells)
    writeLines(colnames(mat)[1:3], file.path(src_dir, "barcodes.tsv"))
    out <- file.path(tempdir(), "tiny_10x_mm_out")

    # mtxInput happily reads the truncated sidecar (n_cells = 3); the
    # disagreement is caught when storeRead opens the mtx header.
    inp <- mtxInput(mtx_path = file.path(src_dir, "matrix.mtx"),
                    feature_id_col = 2L)
    expect_error(
        storeWrite(parquetExprStore(path = out), inp),
        "disagrees with sidecar metadata"
    )
    unlink(src_dir, recursive = TRUE)
})


# ---- exprObj + giotto integration -----------------------------------------

test_that("parquetExprStore can be embedded in an exprObj@exprMat slot", {
    skip_if_not_installed("GiottoClass")
    mat <- .tiny_mat()
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    # Construct directly via new() — createExprObj would route through
    # .evaluate_expr_matrix which doesn't yet recognize parquetExprStore.
    eo  <- new("exprObj", name = "raw", exprMat = pe,
               spat_unit = "cell", feat_type = "rna")
    expect_s4_class(eo, "exprObj")
    expect_s4_class(slot(eo, "exprMat"), "parquetExprStore")
})

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


# ---- storeWrite from a (possibly subset) parquetExprStore / union --------

test_that("storeWrite from single-store subset round-trips to a new store", {
    mat <- Matrix::sparseMatrix(
        i = c(1, 2, 3, 1), j = c(1, 2, 3, 4), x = c(11, 12, 13, 14),
        dims = c(3, 4),
        dimnames = list(c("g1", "g2", "g3"), c("c1", "c2", "c3", "c4"))
    )
    pe <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)
    sub <- pe[, c("c2", "c3")]
    out <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), sub)
    expect_equal(out@n_cells, 2L)
    expect_equal(out@cell_ids, c("c2", "c3"))
    # Reconstruct as matrix from the new store's read; compare to in-mem subset
    df <- as.data.frame(dplyr::collect(storeRead(out)))
    rt <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(out@n_genes, out@n_cells),
        dimnames = list(out@feat_ids, out@cell_ids)
    )
    expect_equal(as.matrix(rt), as.matrix(mat[, c("c2", "c3")]))
})

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


test_that("parquetExprStore swaps into a giotto object via setExpression", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("GiottoClass")
    mat <- .tiny_mat(n_genes = 8, n_cells = 12, density = 0.5, seed = 9)
    pe  <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), mat)

    # Build a giotto skeleton with the in-memory matrix, then swap pe in.
    g  <- Giotto::createGiottoObject(expression = mat, verbose = FALSE)
    eo <- new("exprObj", name = "raw", exprMat = pe,
              spat_unit = "cell", feat_type = "rna")
    g  <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    # @exprMat is now the parquetExprStore; structural shape matches.
    em <- GiottoClass::getExpression(g)
    expect_s4_class(slot(em, "exprMat"), "parquetExprStore")
    expect_equal(nrow(em), nrow(mat))
    expect_equal(ncol(em), ncol(mat))
    expect_equal(rownames(em), rownames(mat))
    expect_equal(colnames(em), colnames(mat))
})
