#' @include class-parquetExprStore.R
NULL

# Optional compiled kernels (GiottoKernels). A caller asks .gd_kernel(name)
# and runs its own R path when the answer is FALSE, so the package works, with
# the same results up to rounding, whether or not GiottoKernels is installed.
# `giottodisk.use_kernels = FALSE` forces the R paths.
.gd_kernel <- function(name) {
    isTRUE(getOption("giottodisk.use_kernels", TRUE)) &&
        .gd_kernels_installed() &&
        GiottoKernels::has_kernel(name)
}

# A missing or outdated GiottoKernels is mentioned once per session, since the
# R path is only slower. `pkg` is a seam for the tests.
.gd_kernels_installed <- function(pkg = "GiottoKernels") {
    GiottoUtils::package_check(pkg,
        repository = sprintf("r-universe:giotto-suite/%s>=0.3.0", pkg),
        on_missing = "message",
        once = TRUE,
        custom_msg = sprintf(paste0(
            "%s is not installed, so this step runs in R: same results, ",
            "slower. To use the compiled kernels:\n",
            "install.packages(\"%s\",\n",
            "    repos = c(\"https://giotto-suite.r-universe.dev\", ",
            "\"https://cloud.r-project.org\"))"
        ), pkg, pkg)
    )
}

# The (row_id, col_id, value) triplets a store's view produces, as the Arrow
# record batch reader the stream kernels take.
.gd_kernel_reader <- function(pe) {
    row_id <- col_id <- value <- NULL  # NSE
    arrow::as_record_batch_reader(
        dplyr::select(storeRead(pe, output = "query"), row_id, col_id, value))
}

# the fork workers' knob sets the thread count too; left at 1, the kernels'
# own option (gkernels.n_threads) decides
.gd_kernel_threads <- function() {
    n <- .par_workers()
    if (n > 1L) n else NULL
}
