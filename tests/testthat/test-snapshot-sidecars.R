# snapshotSave writes a description of the saved object beside it, so a
# snapshot's contents can be read without loading the snapshot.

skip_if_not_installed("GiottoData")

.mk_backed_g <- function(td) {
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    g@source <- gDirSource(td)
    g
}

.expr_names <- function(m) {
    vapply(m$slots$expression, function(r) r$keys$name, character(1L))
}

test_that("snapshotSave writes manifest and history sidecars", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)

    snapdir <- dirname(.snap_file(td, "snap1"))
    expect_true(file.exists(file.path(snapdir, "manifest.json")))
    expect_true(file.exists(file.path(snapdir, "history.ndjson")))
    # the sidecars are not mistaken for the snapshot object
    expect_identical(basename(.snap_file(td, "snap1")), "gobject.RDS")
})

test_that("objManifest reads a snapshot without loading it", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)

    m <- objManifest(g@source, "snap1")
    expect_s3_class(m, "gmanifest")
    expect_identical(m$object$uid, GiottoClass:::.gobject_uid(g))
    # read back in the live form, so it diffs against the object
    expect_false(GiottoClass::manifestDiff(m, objManifest(g))$changed)
    expect_type(objManifest(g@source, "snap1", output = "json"), "character")
})

test_that("objHistory reads a snapshot's history without loading it", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)

    expect_identical(objHistory(g@source, "snap1"), objHistory(g))
})

test_that("a project directory path reads through its source", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)

    expect_identical(objManifest(td, name = "snap1"),
        objManifest(g@source, "snap1"))
    expect_identical(objHistory(td), objHistory(g@source))
})

test_that("the most recent snapshot is used when no name is given", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "older", verbose = FALSE)
    Sys.setFileTime(.snap_file(td, "older"), Sys.time() - 60L)
    g2 <- g
    g2@expression$cell$rna$scaled <- NULL
    snapshotSave(g2@source, g2, name = "newer", verbose = FALSE)

    expect_false("scaled" %in% .expr_names(objManifest(g@source)))
    expect_true("scaled" %in% .expr_names(objManifest(g@source, "older")))
})

test_that("snapshotDelete removes the sidecars with the snapshot", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)
    snapdir <- dirname(.snap_file(td, "snap1"))
    snapshotDelete(g@source, "snap1")

    expect_false(dir.exists(snapdir))
})

test_that("a snapshot written without a manifest reports that clearly", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)
    unlink(file.path(dirname(.snap_file(td, "snap1")), "manifest.json"))

    expect_error(objManifest(g@source, "snap1"), "no manifest.json")
    expect_error(objManifest(g@source, "nope"), "no snapshot named 'nope'")
})

test_that("an overwrite replaces the sidecars with the snapshot", {
    td <- tempfile("sidecars_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    snapshotSave(g@source, g, name = "snap1", verbose = FALSE)
    g2 <- g
    g2@expression$cell$rna$scaled <- NULL
    snapshotSave(g2@source, g2, name = "snap1", overwrite = TRUE,
        verbose = FALSE)

    expect_false("scaled" %in% .expr_names(objManifest(g@source, "snap1")))
})
