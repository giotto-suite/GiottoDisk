# snapshotLoad: error paths (missing dir, unknown name w/ fuzzy match,
# fallback listing) + success paths (valid name, NULL picks most recent).

# Helper: gDirSource with optional set of snapshot names on disk.
.mk_snapshot_src <- function(names = character()) {
    td <- tempfile("snapLoad_"); dir.create(td)
    src <- gDirSource(td)
    dir.create(.snap_gsdir(td))
    for (nm in names) .write_snapshot_fixture(td, nm, list(name = nm))
    list(src = src, td = td)
}


# Error paths ####

test_that("snapshotLoad: errors clearly when no snapshots exist", {
    x <- .mk_snapshot_src()
    on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    expect_error(
        snapshotLoad(x$src, name = "anything"),
        "no snapshots found",
        fixed = TRUE
    )
})

test_that("snapshotLoad: typo close to an available snapshot suggests it", {
    x <- .mk_snapshot_src(c("snap_v1", "other_name"))
    on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    # adist("snap_v2", "snap_v1") = 1; threshold = max(2, ceil(7 * 0.25)) = 2
    err <- expect_error(snapshotLoad(x$src, name = "snap_v2"))
    expect_match(conditionMessage(err), "Did you mean", fixed = TRUE)
    expect_match(conditionMessage(err), "snap_v1",       fixed = TRUE)
})

test_that("snapshotLoad: no close typo falls back to listing available", {
    x <- .mk_snapshot_src(c("snap_v1", "other_name"))
    on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    err <- expect_error(snapshotLoad(x$src, name = "totally_unrelated_zzz"))
    expect_match(conditionMessage(err), "Available",   fixed = TRUE)
    expect_match(conditionMessage(err), "snap_v1",     fixed = TRUE)
    expect_match(conditionMessage(err), "other_name",  fixed = TRUE)
})


# Success path ####

test_that("snapshotLoad: valid name returns the serialized object", {
    x <- .mk_snapshot_src(c("snap_v1"))
    on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    obj <- snapshotLoad(x$src, name = "snap_v1")
    expect_equal(obj, list(name = "snap_v1"))
})

test_that("snapshotLoad: NULL name picks the most recent by mtime", {
    x <- .mk_snapshot_src()
    on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    # Write two snapshots with staggered mtimes -- older, then newer.
    older <- .write_snapshot_fixture(x$td, "snap_old", list(tag = "old"))
    Sys.setFileTime(older, Sys.time() - 60L)
    .write_snapshot_fixture(x$td, "snap_new", list(tag = "new"))
    obj <- snapshotLoad(x$src, name = NULL)
    expect_equal(obj$tag, "new")
})


# qs2 ####

test_that("snapshots save and load with qs2", {
    skip_if_not_installed("qs2")
    skip_if_not_installed("GiottoData")
    withr::local_options(giotto.use_conda = FALSE)
    td <- withr::local_tempdir()
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    g@source <- sourceCreate(td, type = "gDirSource")

    snapshotSave(g@source, g, name = "snap_qs2", method = "qs2",
                 verbose = FALSE)
    snaps <- file.path(td, "giottosave")
    expect_identical(basename(.snap_file(td, "snap_qs2")), "gobject.qs2")
    g2 <- snapshotLoad(g@source, name = "snap_qs2")
    expect_identical(GiottoClass::spatIDs(g2), GiottoClass::spatIDs(g))

    # an older caller hands on its unresolved default vector; it resolves to
    # the first element, as match.arg() resolves a default
    snapshotSave(g@source, g, name = "snap_default",
                 method = c("RDS", "qs"), verbose = FALSE)
    expect_identical(basename(.snap_file(td, "snap_default")), "gobject.RDS")
    expect_setequal(.gdsrc_snapshot_names(snaps),
                    c("snap_qs2", "snap_default"))
})
