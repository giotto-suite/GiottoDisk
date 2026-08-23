# The Atera output layout matches Xenium's today, so `AteraDiskReader`
# subclasses `XeniumDiskReader` and overrides nothing but the platform label.
# These assertions are the contract that keeps that true: if someone later
# copies the Xenium reader instead of extending it, the inheritance check fails.
# They are deliberately data-free so they run in CI.

test_that("AteraDiskReader extends the Xenium disk reader", {
    expect_true(extends("AteraDiskReader", "XeniumDiskReader"))
    expect_true(extends("AteraDiskReader", "XeniumReader"))
})

test_that("AteraDiskReader reports Atera in path-detection messages", {
    # the label is carried by a prototype slot, so no instance is needed
    expect_identical(getClass("AteraDiskReader")@prototype@platform, "Atera")
    # and the parent must be unchanged
    expect_identical(getClass("XeniumDiskReader")@prototype@platform, "Xenium")
})

test_that("importAteraDisk requires a backend", {
    expect_error(importAteraDisk(), "backend` is required")
})
