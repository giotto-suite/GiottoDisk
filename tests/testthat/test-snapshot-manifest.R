# snapshotSave writes a description of the saved object beside it, so a
# snapshot's contents can be read without loading the snapshot.

skip_if_not_installed("GiottoClass")
skip_if_not_installed("GiottoData")

.mk_backed_g <- function(td) {
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    g@source <- gDirSource(td)
    g
}

test_that("snapshotSave writes manifest and history sidecars", {
    td <- tempfile("manifest_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)

    gsdir <- file.path(td, "giottosave")
    expect_true(file.exists(file.path(gsdir, "snap1.manifest.json")))
    expect_true(file.exists(file.path(gsdir, "snap1.history.ndjson")))
    # the snapshot itself is still discoverable by name
    expect_true("snap1" %in% .gdsrc_detect_gsavename(gsdir))
})

test_that("snapshotManifest reads a snapshot without loading it", {
    td <- tempfile("manifest_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)

    m <- snapshotManifest(g@source, "snap1")
    expect_s3_class(m, "gmanifest")
    expect_identical(m$schema_version, "0.1.0")
    expect_identical(m$object$uid, GiottoClass:::.gobject_uid(g))
    expect_identical(m$slots$expression$cell$rna$raw$shape[[2]], ncol(g))
    # accepts a project path as well as a source object
    expect_identical(snapshotManifest(td, "snap1")$object$uid, m$object$uid)
})

test_that("snapshotHistory returns one record per operation", {
    td <- tempfile("manifest_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)

    h <- snapshotHistory(g@source, "snap1")
    expect_length(h, length(GiottoClass::objHistory(g)))
    expect_true(all(c("step_id", "fn", "status") %in% names(h[[1]])))
})

test_that("the most recent snapshot is used when no name is given", {
    td <- tempfile("manifest_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "older", verbose = FALSE)
    Sys.sleep(1.1) # mtime resolution
    g2 <- g
    g2@expression$cell$rna$scaled <- NULL
    snapshotSave(g2@source, g2, name = "newer", verbose = FALSE)

    m <- snapshotManifest(g@source)
    expect_false("scaled" %in% names(m$slots$expression$cell$rna))
    expect_true(
        "scaled" %in% names(
            snapshotManifest(g@source, "older")$slots$expression$cell$rna
        )
    )
})

test_that("snapshotDelete removes the sidecars with the snapshot", {
    td <- tempfile("manifest_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)
    snapshotDelete(g@source, "snap1")

    expect_length(
        list.files(file.path(td, "giottosave"), pattern = "^snap1\\."), 0L
    )
})

test_that("a snapshot written without a manifest reports that clearly", {
    td <- tempfile("manifest_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)
    unlink(file.path(td, "giottosave", "snap1.manifest.json"))

    expect_error(snapshotManifest(g@source, "snap1"), "no manifest")
})
