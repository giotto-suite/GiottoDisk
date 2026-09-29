# The kernel contract (R/kernel-*.R, AGENTS.md "Kernel"): kernel code runs in
# a worker that has not loaded GiottoDisk, so everything it reaches must be
# kernel code, and the fan-out runtime must keep the package off the worker.

# Top-level function names defined in each R/ file. Needs the package source,
# which an installed check does not have.
.defs_by_file <- function() {
    rdir <- testthat::test_path("..", "..", "R")
    testthat::skip_if_not(dir.exists(rdir), "package R/ sources not available")
    out <- character()
    for (fp in list.files(rdir, pattern = "[.]R$", full.names = TRUE)) {
        ln <- readLines(fp, warn = FALSE)
        hit <- regmatches(ln, regexpr("^[.A-Za-z_][.A-Za-z0-9_]* <- function", ln))
        nm <- sub(" <- function", "", hit)
        out[nm] <- basename(fp)
    }
    out
}

test_that("everything a kernel function reaches is kernel code", {
    skip_if_not_installed("codetools")
    defs <- .defs_by_file()
    kernel <- names(defs)[grepl("^kernel-", defs)]
    expect_gt(length(kernel), 10L)
    ns <- asNamespace("GiottoDisk")
    for (f in kernel) {
        b <- .kernel_bundle(get(f, envir = ns))
        reached <- ls(b$env, all.names = TRUE)
        outside <- reached[!grepl("^kernel-", defs[reached])]
        expect_identical(outside, character(0), info = f)
    }
})

test_that("a bundle refuses code that reaches a GiottoDisk generic", {
    skip_if_not_installed("codetools")
    # `unique` is an S4 generic here (the parquetBase distinct op); a worker
    # would resolve the same call to base::unique instead
    f <- function(x) unique(x)
    environment(f) <- asNamespace("GiottoDisk")
    expect_error(.kernel_bundle(f), "S4 generic")
})

test_that("an isolated map runs without GiottoDisk on the workers", {
    skip_on_cran()
    skip_if_not_installed("mirai")
    skip_if_not_installed("codetools")
    probe <- function(x) c(x, "GiottoDisk" %in% loadedNamespaces())
    environment(probe) <- asNamespace("GiottoDisk")
    res <- .isolated_map(list(1, 2, 3), probe, n_workers = 2L, site = "test")
    expect_identical(vapply(res, `[[`, 0, 1L), c(1, 2, 3))
    expect_false(any(vapply(res, `[[`, 0, 2L) == 1))
})

test_that("a failing task names its site and leaves no daemons", {
    skip_on_cran()
    skip_on_os("windows")
    skip_if_not_installed("mirai")
    skip_if_not_installed("codetools")
    n_daemons <- function() {
        length(grep("mirai::daemon\\(", system2("ps", c("-A", "-o", "command="),
            stdout = TRUE), value = TRUE))
    }
    before <- n_daemons()
    boom <- function(x) if (x == 2) stop("kaboom") else x
    environment(boom) <- asNamespace("GiottoDisk")
    expect_error(
        .isolated_map(list(1, 2), boom, n_workers = 2L, site = "test"),
        "[test] fan-out task 2 of 2 failed", fixed = TRUE)
    # shutdown is asynchronous; give the daemons a moment to exit
    for (i in 1:50) {
        if (n_daemons() <= before) break
        Sys.sleep(0.1)
    }
    expect_lte(n_daemons(), before)
})
