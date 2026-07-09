# Tests for snapshotLoad error paths.
#
# The load routine wraps `list.files()` over the source's `giottosave/`
# subdir. Prior to the fix, a typo in `name` returned NA and errored
# unhelpfully inside `.load_serialized`. Now it errors up front with a
# fuzzy-matched suggestion when possible, or lists available snapshots.


# Helper: build a gDirSource with an optional set of snapshot names.
# Each name is saveRDS'd with a trivial object so `.load_serialized` can
# read them back on the success path.
.mk_snapshot_src <- function(names = character()) {
    td <- tempfile("snapLoad_"); dir.create(td)
    snaps_dir <- file.path(td, "giottosave")
    dir.create(snaps_dir)
    for (nm in names) {
        saveRDS(list(name = nm), file.path(snaps_dir, paste0(nm, ".rds")))
    }
    list(src = gDirSource(td), td = td)
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
    snaps_dir <- file.path(x$td, "giottosave")
    # Write two snapshots with staggered mtimes -- older, then newer.
    older <- file.path(snaps_dir, "snap_old.rds")
    newer <- file.path(snaps_dir, "snap_new.rds")
    saveRDS(list(tag = "old"), older)
    Sys.setFileTime(older, Sys.time() - 60L)
    saveRDS(list(tag = "new"), newer)
    obj <- snapshotLoad(x$src, name = NULL)
    expect_equal(obj$tag, "new")
})
