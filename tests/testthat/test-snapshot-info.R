# What a save records about itself (adr/0024): snapshot.json beside the
# object, listed by snapshotInfo(), changed by snapshotMetadata(), and kept
# when the snapshot is deleted.

# callers remove `td` with on.exit()
.info_src <- function() {
    td <- tempfile("info_")
    list(td = td, src = gDirSource(td))
}

# a snapshot as snapshotSave() leaves it, object and snapshot.json
.info_snapshot <- function(td, name, ...) {
    f <- .write_snapshot_fixture(td, name, list(a = 1))
    tag <- basename(dirname(f))
    .ss_gdsrc_write_info(.ss_gdsrc_info(name = name, tag = tag, ...), dirname(f))
    tag
}

describe("snapshotInfo", {
    it("lists what each save recorded, without loading it", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        t1 <- .info_snapshot(x$td, "raw")
        t2 <- .info_snapshot(x$td, "work", keep = FALSE, creator = "GiottoLens",
            session = "s1", parents = t1, labels = c("qc", "draft"),
            extra = list(GiottoLens = list(step = "cluster")))
        info <- snapshotInfo(x$src)
        expect_identical(info$tag, c(t1, t2))
        w <- info[info$name == "work", ]
        expect_identical(w$creator, "GiottoLens")
        expect_identical(w$session, "s1")
        expect_identical(w$parents[[1]], t1)
        expect_false(w$keep)
        expect_identical(w$labels[[1]], c("qc", "draft"))
        expect_identical(w$extra[[1]]$GiottoLens$step, "cluster")
        expect_s3_class(info$created, "POSIXct")
        expect_identical(info$keep, c(TRUE, FALSE))
        expect_identical(nrow(snapshotInfo(x$td, name = "raw")), 1L)
    })

    it("reads a save from before snapshot.json as kept, creator unknown", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        f <- .write_snapshot_fixture(x$td, "old", list(a = 1))
        info <- snapshotInfo(x$src)
        expect_identical(info$name, "old")
        expect_true(info$keep)
        expect_true(is.na(info$creator))
    })

    it("lists deleted snapshots only when asked", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        t1 <- .info_snapshot(x$td, "raw")
        t2 <- .info_snapshot(x$td, "work", parents = t1)
        snapshotDelete(x$src, "raw")
        expect_identical(snapshotInfo(x$src)$tag, t2)
        all <- snapshotInfo(x$src, deleted = TRUE)
        expect_identical(all$deleted, c(TRUE, FALSE))
        # the child still names a folder that exists
        expect_true(dir.exists(file.path(.snap_gsdir(x$td), all$parents[[2]])))
    })
})

describe("snapshotMetadata", {
    it("changes keep, description, labels and extra, nothing else", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        tag <- .info_snapshot(x$td, "work", keep = FALSE, creator = "GiottoLens",
            extra = list(GiottoLens = list(step = "cluster")))
        before <- snapshotInfo(x$src)
        row <- snapshotMetadata(x$src, "work", keep = TRUE,
            description = "leiden at 0.5", labels = "final",
            extra = list(note = "x"))
        expect_true(row$keep)
        expect_identical(row$description, "leiden at 0.5")
        expect_identical(row$labels[[1]], "final")
        expect_identical(names(row$extra[[1]]), c("GiottoLens", "note"))
        expect_identical(row[c("name", "tag", "created", "creator", "session",
            "parents")], before[c("name", "tag", "created", "creator", "session",
            "parents")])
        # left NULL: unchanged
        row <- snapshotMetadata(x$td, "work", labels = c("final", "v2"))
        expect_identical(row$description, "leiden at 0.5")
        expect_identical(row$labels[[1]], c("final", "v2"))
    })

    it("gives a save from before snapshot.json a description", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        f <- .write_snapshot_fixture(x$td, "old", list(a = 1))
        row <- snapshotMetadata(x$src, "old", keep = FALSE)
        expect_false(row$keep)
        expect_true(file.exists(file.path(dirname(f), "snapshot.json")))
    })

    it("refuses a snapshot that does not exist", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        expect_error(snapshotMetadata(x$src, "nope", keep = TRUE), "no snapshot")
    })
})

describe("a deleted snapshot", {
    it("releases its artifacts to sourcePrune, and its name", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        tag <- .info_snapshot(x$td, "qc")
        u <- sourceWrite(x$src, data.frame(a = 1))@uid
        .gdsrc_json_add_tags(x$td, u, tag)
        snapshotDelete(x$src, "qc")
        suppressMessages(sourcePrune(x$src))
        expect_false(dir.exists(.gdsrc_artifact_dir(x$td, u)))
        tag2 <- .info_snapshot(x$td, "qc")
        expect_false(identical(tag2, tag))
        expect_identical(snapshotInfo(x$src)$tag, tag2)
    })
})

describe("the snapshot an object last passed", {
    g0 <- GiottoUtils::gwith_options(list(giotto.use_conda = FALSE),
        GiottoData::loadGiottoMini("visium", verbose = FALSE)
    )
    .save <- function(src, g, name, ...) {
        suppressMessages(snapshotSave(src, g, name = name, verbose = FALSE,
            ...))
    }
    .row <- function(src, name) snapshotInfo(src, name = name)

    it("is set by a save, recorded as the next save's parent, and by a load", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        expect_null(.gdsrc_obj_snapshot(g))

        g <- .save(x$src, g, "first")
        first <- .row(x$src, "first")
        expect_identical(.gdsrc_obj_snapshot(g), first$tag)
        expect_identical(first$object, GiottoClass::objUID(g0))
        expect_length(first$parents[[1]], 0L)

        g <- .save(x$src, g, "second")
        expect_identical(.row(x$src, "second")$parents[[1]], first$tag)

        loaded <- suppressMessages(snapshotLoad(x$src, name = "first"))
        expect_identical(.gdsrc_obj_snapshot(loaded), first$tag)
    })

    it("records every snapshot a combined object came from, then only its own", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        .gdsrc_obj_snapshot(g) <- c("a.aaaaaaaa", "b.bbbbbbbb")
        g <- .save(x$src, g, "joined")
        expect_identical(.row(x$src, "joined")$parents[[1]],
            c("a.aaaaaaaa", "b.bbbbbbbb"))
        expect_identical(.gdsrc_obj_snapshot(g), .row(x$src, "joined")$tag)
        g <- .save(x$src, g, "next")
        expect_identical(.row(x$src, "next")$parents[[1]],
            .row(x$src, "joined")$tag)
    })

    it("gives way to parents the caller names", {
        x <- .info_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        g <- .save(x$src, g, "first")
        .save(x$src, g, "second", parents = c("a.abcdefgh", "b.abcdefgh"))
        expect_identical(.row(x$src, "second")$parents[[1]],
            c("a.abcdefgh", "b.abcdefgh"))
    })
})
