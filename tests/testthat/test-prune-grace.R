# sourcePrune() keeps what a snapshot deleted within the grace period reached
# (adr/0025), and reports what it removes or would remove.

# callers remove `td` with on.exit()
.grace_src <- function() {
    td <- tempfile("grace_")
    list(td = td, src = gDirSource(td))
}

# a snapshot named `name` and a recorded artifact tagged by it
.grace_snapshot <- function(x, name) {
    tag <- basename(dirname(.write_snapshot_fixture(x$td, name, list(a = 1))))
    u <- sourceWrite(x$src, data.frame(a = 1:3))@uid
    .gdsrc_json_consolidate(x$td)
    .gdsrc_json_add_tags(x$td, u, tag)
    list(tag = tag, uid = u)
}

# pretend the snapshot `tag` was deleted `secs` seconds ago
.grace_backdate <- function(td, tag, secs) {
    dir <- file.path(.snap_gsdir(td), tag)
    info <- .gdsrc_stored_info(dir)
    info$deleted_at <- format(Sys.time() - secs, .GDSRC_INFO_TIME, tz = "UTC")
    .ss_gdsrc_write_info(info, dir)
}

.grace_prune <- function(src, ...) suppressMessages(sourcePrune(src, ...))
.grace_exists <- function(x, uid) dir.exists(.gdsrc_artifact_dir(x$td, uid))

describe("a deleted snapshot", {
    it("records when it was deleted", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        .grace_snapshot(x, "qc")
        .grace_snapshot(x, "live")
        snapshotDelete(x$src, "qc")
        info <- snapshotInfo(x$src, deleted = TRUE)
        expect_s3_class(info$deleted_at, "POSIXct")
        at <- info$deleted_at[info$name == "qc"]
        expect_true(abs(as.numeric(Sys.time()) - as.numeric(at)) < 60)
        expect_true(is.na(info$deleted_at[info$name == "live"]))
    })

    it("keeps what it reached through a prune within the grace period", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- .grace_snapshot(x, "qc")
        snapshotDelete(x$src, "qc")
        tab <- .grace_prune(x$src)
        expect_true(.grace_exists(x, s$uid))
        expect_identical(tab$kept_by[tab$uid == s$uid], s$tag)
        # its tag stays on the record, so the next prune keeps it too
        expect_identical(x$src[s$uid, "giottosave"], s$tag)
        .grace_prune(x$src)
        expect_true(.grace_exists(x, s$uid))
    })

    it("keeps what its artifacts depend on", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- .grace_snapshot(x, "qc")
        dep <- sourceWrite(x$src, data.frame(b = 1))@uid
        .gdsrc_json_edit(x$td, uid = s$uid, add = list(depends = dep))
        snapshotDelete(x$src, "qc")
        tab <- .grace_prune(x$src)
        expect_true(.grace_exists(x, dep))
        expect_identical(tab$kept_by[tab$uid == dep], s$tag)
    })

    it("releases them once the grace period has passed", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- .grace_snapshot(x, "qc")
        snapshotDelete(x$src, "qc")
        .grace_backdate(x$td, s$tag, 2 * 86400)
        tab <- .grace_prune(x$src)
        expect_false(.grace_exists(x, s$uid))
        expect_true(is.na(tab$kept_by[tab$uid == s$uid]))
        expect_false(s$uid %in% names(as.list(x$src)))
    })

    it("releases them at once with grace = FALSE", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- .grace_snapshot(x, "qc")
        snapshotDelete(x$src, "qc")
        .grace_prune(x$src, grace = FALSE)
        expect_false(.grace_exists(x, s$uid))
    })

    it("without a recorded time gets no grace period", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- .grace_snapshot(x, "qc")
        # deleted the way GiottoDisk did before it recorded the time
        unlink(.gdsrc_snapshot_objfiles(file.path(.snap_gsdir(x$td), s$tag)))
        .grace_prune(x$src)
        expect_false(.grace_exists(x, s$uid))
    })

    it("replaced by a save of the same name is kept for the grace period", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        old <- .grace_snapshot(x, "qc")
        new <- basename(dirname(.write_snapshot_fixture(x$td, "qc", list(a = 2))))
        .ss_gdsrc_drop_older(.snap_gsdir(x$td), name = "qc", keep = new)
        expect_identical(snapshotInfo(x$src)$tag, new)
        tab <- .grace_prune(x$src)
        expect_true(.grace_exists(x, old$uid))
        expect_identical(tab$kept_by, old$tag)
    })
})

describe("the grace period", {
    it("comes from giottodisk.grace, or from the call", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- .grace_snapshot(x, "qc")
        snapshotDelete(x$src, "qc")
        .grace_backdate(x$td, s$tag, 120)
        .grace_prune(x$src, grace = as.difftime(1, units = "hours"))
        expect_true(.grace_exists(x, s$uid))
        GiottoUtils::gwith_options(list(giottodisk.grace = 60),
            .grace_prune(x$src))
        expect_false(.grace_exists(x, s$uid))
    })

    it("does not cover artifacts no snapshot reached", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        u <- sourceWrite(x$src, data.frame(a = 1))@uid
        .gdsrc_json_consolidate(x$td)
        tab <- .grace_prune(x$src)
        expect_false(.grace_exists(x, u))
        expect_identical(tab$uid, u)
    })

    it("is a non-negative number of seconds, a difftime or a flag", {
        expect_identical(.gdsrc_grace_secs(TRUE), 86400)
        expect_identical(.gdsrc_grace_secs(FALSE), 0)
        expect_identical(.gdsrc_grace_secs(as.difftime(2, units = "mins")), 120)
        expect_error(.gdsrc_grace_secs(-1), "grace")
        expect_error(.gdsrc_grace_secs("1 day"), "grace")
    })
})

describe("sourcePrune(dry_run = TRUE)", {
    it("lists what would go and changes nothing", {
        x <- .grace_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        kept <- .grace_snapshot(x, "qc")
        graced <- .grace_snapshot(x, "old")
        snapshotDelete(x$src, "old")
        u <- sourceWrite(x$src, data.frame(a = 1:10))@uid
        top <- max(.gdsrc_registry_versions(x$td))
        expect_message(tab <- sourcePrune(x$src, dry_run = TRUE),
            "Would prune 1 unprotected artifact .*1 .* kept for snapshots deleted in the last 1 day")
        expect_identical(tab$uid, c(u, graced$uid))
        expect_identical(tab$kept_by, c(NA, graced$tag))
        expect_identical(tab$store, c("parquet", "parquet"))
        expect_true(all(tab$bytes > 0))
        expect_s3_class(tab$created, "POSIXct")
        # the pending edit is read, not folded in
        expect_identical(max(.gdsrc_registry_versions(x$td)), top)
        expect_true(.grace_exists(x, u))
        expect_false(kept$uid %in% tab$uid)
    })
})
