# Snapshot folders (giottosave/<name>.<id>/gobject.<ext>), the tags they
# protect artifacts by, and image tagging on saves after the first.

# callers remove `td` with on.exit()
.layout_src <- function() {
    td <- tempfile("layout_")
    list(td = td, src = gDirSource(td))
}

# a vault artifact tagged with `tags`, as a save would leave it
.layout_artifact <- function(src, tags) {
    st <- sourceWrite(src, data.frame(a = 1:3))
    .gdsrc_json_consolidate(src@path)
    src[st@uid, "giottosave"] <- tags
    st@uid
}

.layout_pruned <- function(src, uid) {
    suppressMessages(sourcePrune(src))
    !dir.exists(.gdsrc_artifact_dir(src@path, uid))
}

describe("snapshot folders", {
    it("are named <name>.<id> and listed by name", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        f <- .write_snapshot_fixture(x$td, "qc", list(a = 1))
        tag <- basename(dirname(f))
        expect_match(tag, "^qc\\.[a-z0-9]{8}$")
        expect_identical(.gdsrc_snapshot_names(.snap_gsdir(x$td)), "qc")
        expect_identical(.gdsrc_snapshot_tags(.snap_gsdir(x$td)), tag)
    })

    it("keep a name that itself contains dots", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        .write_snapshot_fixture(x$td, "run.v2", list(a = 1))
        expect_identical(.gdsrc_snapshot_names(.snap_gsdir(x$td)), "run.v2")
    })

    it("are not listed while staging or without an object file", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        gsdir <- .snap_gsdir(x$td)
        dir.create(file.path(gsdir, ".staging_qc.abcdefgh"), recursive = TRUE)
        saveRDS(1, file.path(gsdir, ".staging_qc.abcdefgh", "gobject.RDS"))
        dir.create(file.path(gsdir, "half.abcdefgh"))
        expect_length(.gdsrc_snapshot_tags(gsdir), 0L)
    })

    it("reject a name that is not one path component", {
        expect_error(.gdsrc_assert_snapshot_name("a/b"), "path separator")
        expect_error(.gdsrc_assert_snapshot_name(".hidden"), "start with")
        expect_silent(.gdsrc_assert_snapshot_name("ok_name.1"))
    })
})

describe("snapshot tags", {
    it("protect an artifact while their folder exists", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        f <- .write_snapshot_fixture(x$td, "qc", list(a = 1))
        uid <- .layout_artifact(x$src, basename(dirname(f)))
        expect_false(.layout_pruned(x$src, uid))
    })

    it("stop protecting once the snapshot is deleted", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        f <- .write_snapshot_fixture(x$td, "qc", list(a = 1))
        uid <- .layout_artifact(x$src, basename(dirname(f)))
        snapshotDelete(x$src, "qc")
        expect_length(.gdsrc_snapshot_names(.snap_gsdir(x$td)), 0L)
        expect_true(.layout_pruned(x$src, uid))
    })

    it("are not revived by a new snapshot of the same name", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        f <- .write_snapshot_fixture(x$td, "qc", list(a = 1))
        uid <- .layout_artifact(x$src, basename(dirname(f)))
        snapshotDelete(x$src, "qc")
        .write_snapshot_fixture(x$td, "qc", list(a = 2))
        expect_true(.layout_pruned(x$src, uid))
    })
})

describe("snapshotSave", {
    skip_if_not_installed("GiottoData")
    g0 <- GiottoUtils::gwith_options(list(giotto.use_conda = FALSE),
        GiottoData::loadGiottoMini("visium", verbose = FALSE)
    )

    it("writes a folder holding gobject.RDS and tags with its name", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        snapshotSave(x$src, g, name = "snap", verbose = FALSE)
        f <- .snap_file(x$td, "snap")
        expect_identical(basename(f), "gobject.RDS")
        expect_identical(dirname(dirname(f)), .snap_gsdir(x$td))
    })

    it("leaves one folder per name on overwrite, under a new tag", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        snapshotSave(x$src, g, name = "snap", verbose = FALSE)
        tag1 <- .snap_tag(x$td, "snap")
        snapshotSave(x$src, g, name = "snap", overwrite = TRUE,
            verbose = FALSE)
        tags <- .gdsrc_snapshot_tags(.snap_gsdir(x$td))
        expect_length(tags, 1L)
        expect_false(identical(tags, tag1))
    })

    it("tags an image already in the vault on every save", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        orig <- terra::sources(g[["images"]][[1]]@raster_object)
        orig_md5 <- tools::md5sum(orig)

        g <- snapshotSave(x$src, g, name = "v1", verbose = FALSE)
        img <- g[["images"]][[1]]
        expect_true(sourceContains(x$src, img@raster_object))
        uid <- .ss_gdsrc_detect_uid_images(g)[[1]]

        n_vault <- length(list.dirs(.gdsrc_vault_dir(x$td), recursive = FALSE))
        g <- snapshotSave(x$src, g, name = "v2", verbose = FALSE)
        tags <- unlist(x$src[uid, "giottosave"])
        expect_true(.snap_tag(x$td, "v2") %in% tags)
        # the second save adopts nothing new
        expect_length(list.dirs(.gdsrc_vault_dir(x$td), recursive = FALSE),
            n_vault)

        # the v1 / v2 / delete v1 / prune sequence keeps the image v2 needs
        snapshotDelete(x$src, "v1")
        suppressMessages(sourcePrune(x$src))
        expect_true(dir.exists(.gdsrc_artifact_dir(x$td, uid)))
        # the dataset's own file was copied, never moved
        expect_identical(tools::md5sum(orig), orig_md5)
    })

    it("does not tag an image outside the vault", {
        x <- .layout_src()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        g <- g0
        g@source <- x$src
        g <- snapshotSave(x$src, g, name = "v1", export_image = FALSE,
            verbose = FALSE)
        expect_length(.ss_gdsrc_detect_uid_images(g), 0L)
    })
})
