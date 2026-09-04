# Shared harness for the `output = "query"` vs `output = "duckdb"` equivalence
# checks, used by test-parquetExprStore.R and test-parquetExprStore-subset.R.
#
# The two outputs are not separate implementations: `.pestore_to_duckdb()`
# rebuilds the base scan from duckdb's `read_parquet` and then applies the same
# `.pe_apply_axis_pred()` / `.pe_apply_ops()` the Acero path applies, because
# both are dplyr. So these live beside the Acero tests rather than in a suite of
# their own -- a divergence means a modification reached one carrier and not the
# other, which is a fact about the shared code, not about duckdb.

skip_if_no_duckdb <- function() {
    testthat::skip_if_not_installed("duckdb")
    testthat::skip_if_not_installed("dbplyr")
}

# Collect both carriers and compare after a common sort. Neither engine
# promises row order -- duckdb's parallel scans and hash joins reorder freely.
#
# `min_rows` is load-bearing rather than defensive: a @uid that does not match
# the on-disk source_id partition yields an empty scan, and every comparison
# below would then pass against two empty frames.
.expect_dd_parity <- function(pe, min_rows = 1L, tolerance = 1e-9) {
    q <- data.table::as.data.table(
        dplyr::collect(storeRead(pe, output = "query")))
    d <- data.table::as.data.table(
        dplyr::collect(storeRead(pe, output = "duckdb")))

    testthat::expect_gte(nrow(q), min_rows)
    testthat::expect_equal(nrow(d), nrow(q))
    testthat::expect_setequal(names(d), names(q))

    data.table::setcolorder(d, names(q))
    k <- intersect(c("source_id", "row_id", "col_id"), names(q))
    data.table::setorderv(q, k)
    data.table::setorderv(d, k)
    testthat::expect_equal(as.data.frame(d), as.data.frame(q),
        tolerance = tolerance)
    invisible(q)
}
