# AteraDiskReader is a reader of its own, not a Xenium subclass, so the two
# technologies can diverge. It shares lower-level loaders with Xenium only
# where the data type and its reading are the same. Data-free checks; the
# end-to-end import is in test-zarr-atera-v6.R.

test_that("AteraDiskReader is standalone", {
    expect_false(extends("AteraDiskReader", "XeniumDiskReader"))
    expect_false(extends("AteraDiskReader", "XeniumReader"))
})

test_that("importAteraDisk requires a backend", {
    expect_error(importAteraDisk(), "backend` is required")
})
