# Run on library loading


.onAttach <- function(libname, pkgname) {
    # check_ver <- getOption("giotto.check_version", TRUE)
    # if (isTRUE(check_ver)) {
    #     GiottoUtils::check_github_suite_ver("GiottoDisk")
    #     options("giotto.check_version" = FALSE)
    # }

    # initialize options
    # -- default formats
    init_option("giotto.gdsrc_matrix_format", "h5")
    init_option("giotto.gdsrc_spatvector_format", "parquetGeom")
    init_option("giotto.gdsrc_dataframe_format", "parquet")
    # -- random IDs
    init_option("giottodisk.uid_include_node", FALSE)
    init_option("giottodisk.uid_include_pid", TRUE)
    # -- gDirSource
    init_option("giottodisk.use_locking", TRUE)
}
