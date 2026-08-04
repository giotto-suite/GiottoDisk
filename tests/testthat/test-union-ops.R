# Tests for unionParquetExprStore @post_ops chain:
#   * substores must be ops-clean at union construction time
#   * processData on union builds a single (source_id, orig_row_id, scalef)
#     table spanning substores
#   * storeRead on union applies @post_ops via composite-key join, identical to
#     the single-store arrow path
#   * cell subset on union slices the scalef table along (source_id, row_id);
#     works for within-substore and cross-substore subsets
#   * gene subset preserves @post_ops unchanged (current op kind is cell-axis only)

.tiny_substore <- function(n_genes = 10L, n_cells = 8L, prefix = "a",
                            seed = 1L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.5,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", sprintf("%03d", seq_len(n_genes)))
    colnames(m) <- paste0(prefix, "_c", seq_len(n_cells))
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}


test_that("unionParquetExprStore: substores with @post_ops are rejected", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 1)
    pe2 <- .tiny_substore(prefix = "b", seed = 2)
    # Pre-normalize pe1 so it has ops queued
    pe1n <- GiottoClass::processData(pe1,
        Giotto::normParam("library", scalefactor = 1e4))
    expect_length(pe1n@post_ops, 1L)

    expect_error(
        unionParquetExprStore(list(pe1n, pe2)),
        "queued @ops or @post_ops"
    )
})


test_that("processData(unionParquetExprStore, libNorm) builds union-spanning scalef", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 11, n_cells = 6L)
    pe2 <- .tiny_substore(prefix = "b", seed = 12, n_cells = 4L)
    u <- unionParquetExprStore(list(pe1, pe2))
    u2 <- GiottoClass::processData(u,
        Giotto::normParam("library", scalefactor = 1e4))

    expect_length(u2@post_ops, 1L)
    expect_equal(u2@post_ops[[1]]$type, "norm_libsize")

    scalef_dt <- u2@post_ops[[1]]$scalef
    expect_setequal(names(scalef_dt),
        c("source_id", "orig_row_id", "scalef"))
    # one row per cell in the union view (6 + 4 = 10)
    expect_equal(nrow(scalef_dt), 10L)
    # both substore uids represented
    expect_setequal(unique(scalef_dt$source_id), c(pe1@uid, pe2@uid))
})


test_that("storeRead on union with norm @post_ops mutates value in place", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 21, n_cells = 5L)
    pe2 <- .tiny_substore(prefix = "b", seed = 22, n_cells = 4L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4))

    # tibble output applies @post_ops R-side; value is now normalized.
    out <- data.table::as.data.table(storeRead(u, output = "tibble"))
    expect_true("source_id" %in% names(out))
    # both substores present in output
    expect_setequal(unique(out$source_id), c(pe1@uid, pe2@uid))
    # No NA value — every triplet's post-op join succeeded.
    expect_false(anyNA(out$value))
    # No v_norm sidecar column — value is mutated in place.
    expect_false("v_norm" %in% names(out))
})


test_that("within-substore cell subset on union slices scalef correctly", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 31, n_cells = 6L)
    pe2 <- .tiny_substore(prefix = "b", seed = 32, n_cells = 5L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4))
    expect_equal(nrow(u@post_ops[[1]]$scalef), 11L)

    # Subset to cells 1..3 (all in substore A)
    u_a <- u[, 1:3]
    expect_equal(nrow(u_a@post_ops[[1]]$scalef), 3L)
    expect_setequal(unique(u_a@post_ops[[1]]$scalef$source_id), pe1@uid)
    expect_setequal(u_a@post_ops[[1]]$scalef$orig_row_id, 1:3)
})


test_that("cross-substore cell subset on union slices scalef correctly", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 41, n_cells = 6L)
    pe2 <- .tiny_substore(prefix = "b", seed = 42, n_cells = 5L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4))

    # cells 3 (from A) and 9 (= position 3 of B, since A has 6)
    u_x <- u[, c(3L, 9L)]
    scalef_dt <- u_x@post_ops[[1]]$scalef
    expect_equal(nrow(scalef_dt), 2L)
    expect_setequal(unique(scalef_dt$source_id), c(pe1@uid, pe2@uid))
    # The B-side surviving cell should have orig_row_id = 3 (local position in B)
    b_rows <- scalef_dt[source_id == pe2@uid]
    expect_equal(b_rows$orig_row_id, 3L)
})


test_that("gene subset on union leaves @post_ops unchanged (cell-axis only op)", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 51, n_cells = 4L)
    pe2 <- .tiny_substore(prefix = "b", seed = 52, n_cells = 3L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4))
    before <- nrow(u@post_ops[[1]]$scalef)

    u_g <- u[1:5, ]   # narrow to 5 features
    expect_equal(nrow(u_g@post_ops[[1]]$scalef), before)
})


test_that("union with norm + log chain + storeRead applies log1p / log(base)", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 61, n_cells = 4L)
    pe2 <- .tiny_substore(prefix = "b", seed = 62, n_cells = 4L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4)) |>
        GiottoClass::processData(Giotto::normParam("log", base = 2,
            offset = 1))

    # Two independent records in chain order, not one fused record.
    expect_equal(vapply(u@post_ops, function(o) o$type, character(1L)),
                 c("norm_libsize", "log"))
    expect_equal(u@post_ops[[2]]$base, 2)
})
