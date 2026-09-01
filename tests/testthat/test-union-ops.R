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


test_that("unionParquetExprStore: substores with queued ops are rejected", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 1)
    pe2 <- .tiny_substore(prefix = "b", seed = 2)
    # Pre-normalize pe1 so it has ops queued
    pe1n <- GiottoClass::processData(pe1,
        Giotto::normParam("library", scalefactor = 1e4))
    expect_length(pe1n@ops, 1L)

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

    expect_length(u2@ops, 1L)
    expect_equal(u2@ops[[1]]$type, "multiply")

    # One vector per substore, each indexed by that substore's on-disk row_id.
    factors <- u2@ops[[1]]$factors
    expect_setequal(names(factors), c(pe1@uid, pe2@uid))
    expect_length(factors[[pe1@uid]], 6L)
    expect_length(factors[[pe2@uid]], 4L)
    # every cell in the union view has a factor (6 + 4 = 10)
    expect_equal(sum(vapply(factors, function(v) sum(!is.na(v)), numeric(1L))),
                 10)
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


test_that("within-substore cell subset on union leaves the payload intact", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 31, n_cells = 6L)
    pe2 <- .tiny_substore(prefix = "b", seed = 32, n_cells = 5L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4))
    before <- u@ops[[1]]$factors

    # Subset to cells 1..3 (all in substore A). Factors are keyed by on-disk
    # id, so narrowing the view cannot misalign them -- the payload is carried
    # through untouched, and the dropped substore's vector simply goes unread.
    u_a <- u[, 1:3]
    expect_identical(u_a@ops[[1]]$factors[[pe1@uid]], before[[pe1@uid]])
})


test_that("cross-substore cell subset on union keeps both vectors", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 41, n_cells = 6L)
    pe2 <- .tiny_substore(prefix = "b", seed = 42, n_cells = 5L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4))

    # cells 3 (from A) and 9 (= position 3 of B, since A has 6)
    u_x <- u[, c(3L, 9L)]
    factors <- u_x@ops[[1]]$factors
    expect_setequal(names(factors), c(pe1@uid, pe2@uid))
    # The surviving B cell reads its factor at on-disk id 3, its position in B
    expect_equal(factors[[pe2@uid]][3L],
                 u@ops[[1]]$factors[[pe2@uid]][3L])
})


test_that("gene subset on union leaves the op chain unchanged (cell-axis op)", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 51, n_cells = 4L)
    pe2 <- .tiny_substore(prefix = "b", seed = 52, n_cells = 3L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4))
    before <- u@ops[[1]]$factors

    u_g <- u[1:5, ]   # narrow to 5 features
    # cell-axis payload, so a gene subset cannot touch it
    expect_identical(u_g@ops[[1]]$factors, before)
})


test_that("union with norm + log chain records two ops in chain order", {
    skip_if_not_installed("Giotto")
    pe1 <- .tiny_substore(prefix = "a", seed = 61, n_cells = 4L)
    pe2 <- .tiny_substore(prefix = "b", seed = 62, n_cells = 4L)
    u <- unionParquetExprStore(list(pe1, pe2)) |>
        GiottoClass::processData(Giotto::normParam("library",
            scalefactor = 1e4)) |>
        GiottoClass::processData(Giotto::normParam("log", base = 2,
            offset = 1))

    # Two independent records in chain order, not one fused record.
    expect_equal(vapply(u@ops, function(o) o$type, character(1L)),
                 c("multiply", "log"))
    expect_equal(u@ops[[2]]$base, 2)
})
