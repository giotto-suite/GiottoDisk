# What a project uses and what a prune would free: snapshotInfo(sizes = TRUE)
# and sourceInfo() (GiottoDisk#159).

# callers remove `x$td` with on.exit()
.info_src <- function() {
    td <- tempfile("info_")
    list(td = td, src = gDirSource(td))
}

# a snapshot called `name`; returns its tag
.info_snapshot <- function(x, name) {
    basename(dirname(.write_snapshot_fixture(x$td, name, list(a = 1))))
}

# a recorded artifact of `n` rows, tagged by `tags`; returns its uid
.info_artifact <- function(x, n, tags = character()) {
    u <- sourceWrite(x$src, data.frame(a = seq_len(n)))@uid
    .gdsrc_json_consolidate(x$td)
    for (t in tags) .gdsrc_json_add_tags(x$td, u, t)
    u
}

.info_bytes <- function(x, uids) {
    sum(vapply(uids, function(u) .gdsrc_dir_bytes(.gdsrc_artifact_dir(x$td, u)), 0))
}

describe("snapshotInfo(sizes = TRUE)", {
    it("splits what each snapshot reaches into its own and shared bytes", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        a <- .info_snapshot(x, "a")
        b <- .info_snapshot(x, "b")
        only_a <- .info_artifact(x, 10, a)
        both <- .info_artifact(x, 1000, c(a, b))
        only_b <- .info_artifact(x, 100, b)
        # reached by `a` only through a dependency of its own artifact
        dep <- .info_artifact(x, 5)
        .gdsrc_json_edit(x$td, uid = only_a, add = list(depends = dep))

        info <- snapshotInfo(x$src, sizes = TRUE)
        ra <- info[info$name == "a", ]
        rb <- info[info$name == "b", ]
        expect_identical(ra$own_bytes, .info_bytes(x, c(only_a, dep)))
        expect_identical(ra$shared_bytes, .info_bytes(x, both))
        expect_identical(rb$own_bytes, .info_bytes(x, only_b))
        expect_identical(rb$shared_bytes, .info_bytes(x, both))
        expect_identical(ra$object_bytes, .gdsrc_dir_bytes(file.path(.snap_gsdir(x$td), a)))
        expect_true(ra$object_bytes > 0)
        # without sizes, no size columns
        expect_false("own_bytes" %in% names(snapshotInfo(x$src)))
    })

    it("gives a snapshot's bytes to the other once it is deleted", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        a <- .info_snapshot(x, "a")
        b <- .info_snapshot(x, "b")
        both <- .info_artifact(x, 100, c(a, b))
        snapshotDelete(x$src, "a")
        info <- snapshotInfo(x$src, deleted = TRUE, sizes = TRUE)
        expect_identical(info$own_bytes[info$name == "b"], .info_bytes(x, both))
        # the deleted one still carries its tag, shared with `b`
        expect_identical(info$own_bytes[info$name == "a"], 0)
        expect_identical(info$shared_bytes[info$name == "a"], .info_bytes(x, both))
        suppressMessages(sourcePrune(x$src, grace = FALSE))
        info <- snapshotInfo(x$src, deleted = TRUE, sizes = TRUE)
        expect_identical(info$shared_bytes[info$name == "a"], 0)
    })

    it("has the size columns on a project with no snapshots", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        info <- snapshotInfo(x$src, sizes = TRUE)
        expect_identical(nrow(info), 0L)
        expect_true(all(c("object_bytes", "own_bytes", "shared_bytes") %in% names(info)))
    })
})

describe("sourceInfo()", {
    it("sums the vault, counts snapshots and says what a prune would free", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        keep <- .info_snapshot(x, "keep")
        temp <- .info_snapshot(x, "temp")
        snapshotMetadata(x$src, "temp", keep = FALSE)
        gone <- .info_snapshot(x, "gone")
        held <- .info_artifact(x, 50, keep)
        graced <- .info_artifact(x, 20, gone)
        loose <- .info_artifact(x, 30)
        snapshotDelete(x$src, "gone")

        s <- sourceInfo(x$src)
        expect_s3_class(s, "gdsrcInfo")
        expect_identical(s$bytes[["vault"]], .info_bytes(x, c(held, graced, loose)))
        expect_identical(s$bytes[["images"]], 0)
        expect_true(s$bytes[["snapshots"]] > 0)
        expect_identical(s$snapshots,
            c(live = 2L, kept = 1L, unkept = 1L, deleted = 1L))
        expect_identical(s$artifacts, 3L)
        expect_identical(s$prunable$n, 1L)
        expect_identical(s$prunable$bytes, .info_bytes(x, loose))
        expect_s3_class(s$prunable$newest, "POSIXct")
        expect_identical(s$graced$n, 1L)
        expect_identical(s$graced$bytes, .info_bytes(x, graced))
        # it agrees with the prune preview
        tab <- suppressMessages(sourcePrune(x$src, dry_run = TRUE))
        expect_identical(sum(is.na(tab$kept_by)), s$prunable$n)

        s0 <- sourceInfo(x$src, grace = FALSE)
        expect_identical(s0$prunable$n, 2L)
        expect_identical(s0$graced$n, 0L)
        expect_output(print(s), "1 deleted")
        expect_identical(sourceInfo(x$td)$artifacts, 3L)
    })

    it("reports an empty project", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- sourceInfo(x$src)
        expect_identical(s$artifacts, 0L)
        expect_true(is.na(s$prunable$newest))
        expect_output(print(s), "prunable +: none")
    })
})
