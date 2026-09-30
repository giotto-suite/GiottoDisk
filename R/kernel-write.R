# Kernel: parquet writes -- codec options and the per-window write task.
#
# R/kernel-*.R is the part of GiottoDisk that runs where GiottoDisk is not
# loaded. `.kernel_bundle()` (R/utils-isolate.R) copies these functions out of
# the namespace for a worker that has only arrow, dplyr and data.table, so each
# one takes plain data (op records, axis plans, lookup tables, paths) plus a
# lazy query or a data.table, and:
#   * calls other packages only as `pkg::fn`;
#   * calls no GiottoDisk function outside R/kernel-*.R;
#   * reaches no S4 generic or class -- including a base name GiottoDisk masks
#     with a generic, such as `unique`, which is why the kernel writes
#     `base::unique`;
#   * reads options only with a default (the worker gets a giottodisk.*
#     snapshot, not the parent's session).
# test-kernel-isolation.R enforces the list. See AGENTS.md, "Kernel".

# Global arrow writer compression settings -------------------------------
#
# Option `giottodisk.parquet_compression`: codec passed to arrow's parquet
# writer. One of "snappy" (default), "zstd", "gzip", "lz4", "brotli", or
# "uncompressed". Snappy keeps reads fast on local NVMe; zstd buys ~30-40%
# smaller files at a small CPU cost — better when disk space matters or
# when shipping artifacts.
# Option `giottodisk.parquet_compression_level`: numeric. Codec-specific
# level. For zstd, 1-22 (default 3 if unset); for gzip, 0-9. Ignored for
# snappy / uncompressed.
# Both options registered in zzz.R via init_option().
.parquet_compression <- function() {
    codec <- getOption("giottodisk.parquet_compression", default = "snappy")
    checkmate::assert_string(codec)
    codec
}

.parquet_compression_level <- function() {
    lvl <- getOption("giottodisk.parquet_compression_level", default = NULL)
    if (is.null(lvl)) return(NULL)
    checkmate::assert_number(lvl)
    as.integer(lvl)
}

# arrow keeps a data.frame's R attributes in the parquet metadata and restores
# them on collect. data.table's `sorted` and `index` describe the object in
# memory, not the file: restored onto a read that spans several files, or
# after a filter, they are stale, and data.table then answers keyed and
# indexed subsets from them silently wrong. Rebuild the table over the same
# columns (no data copied) so neither is written.
.drop_dt_state <- function(x) {
    if (data.table::is.data.table(x)) x <- data.table::setDT(as.list(x))
    x
}

# Internal wrapper: arrow::write_parquet that honors the GiottoDisk-global
# compression options. Direct callers should use this instead of
# arrow::write_parquet so the codec is consistent across the package.
.write_parquet_file <- function(x, sink, ...) {
    x <- .drop_dt_state(x)
    args <- list(...)
    if (is.null(args$compression)) {
        args$compression <- .parquet_compression()
    }
    lvl <- .parquet_compression_level()
    if (is.null(args$compression_level) && !is.null(lvl)) {
        args$compression_level <- lvl
    }
    do.call(arrow::write_parquet, c(list(x = x, sink = sink), args))
}

# One window of the lazy-chain write from lowered data rather than a store:
# the same scan, remap, shift and sort as `.pestore_write_windows_serial()`.
# Runs in a worker without GiottoDisk, so it returns only the row count.
.pestore_write_window_task <- function(task) {
    row_id <- NULL  # NSE
    sc <- task$scan
    q <- .pe_compose_scan(.pe_read_dataset(sc$path), sc$ci_plan, sc$gi_plan,
        sc$ops)
    off <- task$off
    tb <- .pestore_apply_remap(q, task$luts) |>
        dplyr::mutate(row_id = row_id + off) |>
        dplyr::compute()
    n <- tb$num_rows
    if (n > 0L) .write_parquet_file(tb, task$file, chunk_size = 1048576L)
    # same reason as the serial loop: R does not see Arrow's allocation
    rm(tb)
    invisible(gc(verbose = FALSE))
    n
}
