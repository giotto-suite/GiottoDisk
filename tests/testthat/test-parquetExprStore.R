# Tests for parquetExprStore: storage class for streaming-friendly
# long-format expression matrices.

# helper: build a tiny labeled dgCMatrix
.tiny_mat <- function(n_genes = 5L, n_cells = 8L, density = 0.4, seed = 1L)
    .rand_counts(n_genes, n_cells, density, seed,
        gene_prefix = "gene", cell_prefix = "cell")

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

test_that("storeRead(output = 'duckdb') replaces a reused name, never reuses it", {
    # `name` labels the view; it is not a cache key. Every call compiles a
    # fresh scan and installs it with CREATE OR REPLACE, so a second call under
    # the same name repoints the tbl_dbi the first call handed out.
    skip_if_no_duckdb()
    conn <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    a <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat(n_cells = 8L, seed = 1))
    b <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .tiny_mat(n_cells = 20L, seed = 2))

    t1 <- storeRead(a, output = "duckdb",
        duckdb_params = list(conn = conn, name = "shared"))
    n_before <- nrow(dplyr::collect(t1))

    t2 <- storeRead(b, output = "duckdb",
        duckdb_params = list(conn = conn, name = "shared"))

    expect_false(n_before == nrow(dplyr::collect(t2)))
    # the first handle now reads the second store
    expect_equal(nrow(dplyr::collect(t1)), nrow(dplyr::collect(t2)))
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

# ---- storeWrite from a subset parquetExprStore ----------------------------

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
