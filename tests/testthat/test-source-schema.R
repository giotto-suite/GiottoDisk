# Project directory schema: version stamping, the refusal outside the
# supported range, and the v1 -> v2 upgrade.

# a project as GiottoDisk with schema v1 left it: flat snapshot files and
# artifacts tagged with the bare snapshot name. Callers remove `td`.
.v1_project <- function() {
    td <- tempfile("schema_v1_")
    src <- gDirSource(td)
    a <- sourceWrite(src, data.frame(a = 1:3))@uid
    b <- sourceWrite(src, data.frame(b = 1:3))@uid
    gsdir <- .snap_gsdir(td)
    dir.create(gsdir)
    saveRDS(list(snap = "qc"), file.path(gsdir, "qc.rds"))
    saveRDS(list(snap = "old"), file.path(gsdir, "old.qs")) # unreadable now
    manifest <- .gdsrc_json_read(td, consolidate = TRUE)$content
    manifest[[a]]$giottosave <- "qc"
    manifest[[b]]$giottosave <- "old"
    .gdsrc_json_write(td, manifest, version = 1L)
    list(td = td, src = src, a = a, b = b)
}

describe("sourceSchema", {
    it("stamps a new project with the current version", {
        td <- tempfile("schema_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        s <- sourceSchema(gDirSource(td))
        expect_identical(s$version, .GDSRC_SCHEMA_VERSION)
        expect_identical(s$current, .GDSRC_SCHEMA_VERSION)
        expect_true(s$supported)
    })

    it("reports an old project as unsupported", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        s <- sourceSchema(x$td)
        expect_identical(s$version, 1L)
        expect_false(s$supported)
    })
})

describe("an outdated project", {
    it("is refused by the snapshot verbs, pruning and manifest writes", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        expect_error(snapshotLoad(x$src), "sourceUpdate")
        expect_error(snapshotDelete(x$src, "qc"), "sourceUpdate")
        expect_error(sourcePrune(x$src), "sourceUpdate")
        expect_error(x$src[x$a, "giottosave"] <- "x", "sourceUpdate")
        n_vault <- length(list.dirs(.gdsrc_vault_dir(x$td), recursive = FALSE))
        expect_error(sourceWrite(x$src, data.frame(c = 1)), "sourceUpdate")
        expect_error(setArtifactDumpDir(x$src), "sourceUpdate")
        # refused before writing, so nothing reached the vault
        expect_length(list.dirs(.gdsrc_vault_dir(x$td), recursive = FALSE),
            n_vault)
    })

    it("is let through with a warning by the escape hatch, keeping v1", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        GiottoUtils::gwith_options(list(giottodisk.schema_check = FALSE), {
            expect_warning(x$src[x$a, "meta"] <- "x", "schema_check")
            expect_warning(.gdsrc_json_read(x$td, consolidate = TRUE),
                "schema_check")
        })
        expect_identical(sourceSchema(x$td)$version, 1L)
    })

    it("written by a newer GiottoDisk is refused too", {
        td <- tempfile("schema_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        gDirSource(td)
        .gdsrc_json_write(td, list(), version = .GDSRC_SCHEMA_VERSION + 1L)
        expect_error(sourcePrune(gDirSource(td)), "newer GiottoDisk")
        expect_error(sourceUpdate(td), "newer GiottoDisk")
    })
})

describe("sourceUpdate", {
    it("moves flat snapshots into folders and retags their artifacts", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        sourceUpdate(x$td)
        expect_identical(sourceSchema(x$td)$version, 2L)

        gsdir <- .snap_gsdir(x$td)
        expect_setequal(.gdsrc_snapshot_names(gsdir), c("qc", "old"))
        expect_length(list.files(gsdir, pattern = "\\.(rds|qs)$"), 0L)
        expect_identical(basename(.snap_file(x$td, "old")), "gobject.qs")
        expect_identical(snapshotLoad(x$src, name = "qc"), list(snap = "qc"))

        expect_identical(x$src[x$a, "giottosave"], .snap_tag(x$td, "qc"))
        expect_identical(x$src[x$b, "giottosave"], .snap_tag(x$td, "old"))
    })

    it("keeps the artifacts of a snapshot it cannot load", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        sourceUpdate(x$td)
        suppressMessages(sourcePrune(x$src))
        expect_true(dir.exists(.gdsrc_artifact_dir(x$td, x$b)))
        expect_true(dir.exists(.gdsrc_artifact_dir(x$td, x$a)))
    })

    it("finishes an upgrade that was interrupted after moving the files", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        gsdir <- .snap_gsdir(x$td)
        d <- file.path(gsdir, .gdsrc_snapshot_tag("qc"))
        dir.create(d)
        file.rename(file.path(gsdir, "qc.rds"), file.path(d, "gobject.RDS"))

        sourceUpdate(x$td)
        expect_identical(x$src[x$a, "giottosave"], basename(d))
    })

    it("leaves a current project alone", {
        td <- tempfile("schema_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        src <- gDirSource(td)
        before <- tools::md5sum(.gdsrc_json_path(td))
        sourceUpdate(src)
        expect_identical(tools::md5sum(.gdsrc_json_path(td)), before)
    })
})
