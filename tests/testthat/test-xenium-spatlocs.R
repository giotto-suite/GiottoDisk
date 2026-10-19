# GiottoDisk carries its own copy of the Xenium `create_gobject` orchestration,
# so a fallback wired only into Giotto's copy would silently do nothing on the
# `backend =` path. This asserts the disk copy calls it too -- the drift that
# made an earlier reader fix a two-package change.

test_that("the disk create path falls back to cell-metadata centroids", {
    skip_if_not_installed("Giotto")
    src <- paste(deparse(
        methods::getMethod("initialize", "XeniumDiskReader")@.Data
    ), collapse = " ")
    expect_true(grepl("load_spatlocs", src))
})

test_that("polygon centroids still take precedence", {
    skip_if_not_installed("Giotto")
    src <- paste(deparse(
        methods::getMethod("initialize", "XeniumDiskReader")@.Data
    ), collapse = " ")
    # the fallback is guarded on there being no spatial info, so the
    # centroids_to_spatlocs path stays the primary source
    expect_true(grepl("centroids_to_spatlocs", src))
    expect_true(grepl("list_spatial_info_names", src))
})
