# Run on library loading


.onAttach <- function(libname, pkgname) {
    # check_ver <- getOption("giotto.check_version", TRUE)
    # if (isTRUE(check_ver)) {
    #     GiottoUtils::check_github_suite_ver("GiottoDisk")
    #     options("giotto.check_version" = FALSE)
    # }

    # initialize options
    # -- file testing
    init_option("giottodisk.uri_test_exists", FALSE)
    # -- default formats
    init_option("giotto.gdsrc_sparsematrix_format", "parquetExpr")
    init_option("giotto.gdsrc_densematrix_format",  "h5")
    init_option("giotto.gdsrc_spatvector_format",   "parquetGeom")
    init_option("giotto.gdsrc_dataframe_format",    "parquet")
    # -- random IDs
    init_option("giottodisk.uid_include_node", FALSE)
    init_option("giottodisk.uid_include_pid", TRUE)
    # -- gDirSource
    init_option("giottodisk.use_locking", TRUE)
    # -- arrow schema
    init_option("giottodisk.arrow_timestamp_locale", "us")
    # -- parquet writer compression. See .parquet_compression() in
    #    utils-arrow.R. Snappy default keeps reads fast; switch to
    #    "zstd" for ~30-40% smaller files at a small CPU cost.
    init_option("giottodisk.parquet_compression", "snappy")
    # -- dump management
    init_option("giottodisk.artifact_dump", file.path(tempdir(), "gdisk_dump"))
    init_option("giottodisk.adopt_external", FALSE)
    # -- duckdb
    init_option("giottodisk.duckdb_memory_limit", "8GB")
    # Above this many ids, an axis membership predicate is registered and
    # joined rather than inlined by dbplyr as a literal IN list. Distinct from
    # the sedona threshold below: same number, different engine, and a
    # different rewrite (semi/anti join vs a VALUES subquery).
    init_option("giottodisk.duckdb_in_subquery_threshold", 1000L)
    # -- plot
    init_option("giottodisk.plot_sample_max", 1e5)
    # -- sedona SQL translation
    init_option("giottodisk.sedona_in_subquery_threshold", 1000L)
}
