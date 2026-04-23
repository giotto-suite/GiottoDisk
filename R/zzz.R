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
    init_option("giotto.gdsrc_sparsematrix_format", "bpcells")
    init_option("giotto.gdsrc_densematrix_format", "h5")
    init_option("giotto.gdsrc_spatvector_format", "parquetGeom")
    init_option("giotto.gdsrc_dataframe_format", "parquet")
    # -- random IDs
    init_option("giottodisk.uid_include_node", FALSE)
    init_option("giottodisk.uid_include_pid", TRUE)
    # -- gDirSource
    init_option("giottodisk.use_locking", TRUE)
    # -- arrow schema
    init_option("giottodisk.arrow_timestamp_locale", "us")
    # -- dump management
    init_option("giottodisk.artifact_dump", file.path(tempdir(), "gdisk_dump"))
    # -- duckdb
    init_option("giottodisk.duckdb_memory_limit", "8GB")
    # -- plot
    init_option("giottodisk.plot_sample_max", 1e5)
    # -- sedona SQL translation
    init_option("giottodisk.sedona_in_subquery_threshold", 1000L)
}
