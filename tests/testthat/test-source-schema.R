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
    .as_legacy_layout(td, manifest, version = 1L)
    list(td = td, src = src, a = a, b = b)
}

# turn a project into the layout before schema v3: the registry as
# giottodir.json itself, no versions, no v3 edit files
.as_legacy_layout <- function(td, manifest, version) {
    .gdsrc_legacy_write(td, manifest, version = version)
    unlink(file.path(td, list.files(td, pattern = .GDSRC_REGISTRY_RE)))
    unlink(list.files(.gdsrc_json_pending_dir(td), pattern = .GDSRC_EDIT_RE,
        full.names = TRUE))
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
        expect_error(.gdsrc_json_add_tags(x$td, x$a, "x"), "sourceUpdate")
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
            expect_warning(.gdsrc_json_read(x$td, consolidate = TRUE),
                "schema_check")
            # an older layout takes no writes, escape hatch or not
            expect_error(suppressWarnings(.gdsrc_json_add_tags(x$td, x$a, "x")),
                "older layout")
        })
        expect_identical(sourceSchema(x$td)$version, 1L)
    })

    it("written by a newer GiottoDisk is refused too", {
        td <- tempfile("schema_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        gDirSource(td)
        .gdsrc_legacy_write(td, list(), version = .GDSRC_SCHEMA_VERSION + 1L)
        expect_error(sourcePrune(gDirSource(td)), "newer GiottoDisk")
        expect_error(sourceUpdate(td), "newer GiottoDisk")
    })
})

describe("sourceUpdate", {
    it("moves flat snapshots into folders and retags their artifacts", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        sourceUpdate(x$td)
        expect_identical(sourceSchema(x$td)$version, .GDSRC_SCHEMA_VERSION)

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

    it("moves a v2 registry into versioned files, with its pending edits", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        .gdsrc_upgrade_1_to_2(x$td)
        # an edit as GiottoDisk < v3 queued it: one file per uid, merged
        d <- .gdsrc_json_pending_dir(x$td)
        dir.create(d, showWarnings = FALSE)
        jsonlite::write_json(list(content = list(giottosave = c(
            .snap_tag(x$td, "qc"), "late.abcdefgh"))), file.path(d, paste0(x$a, ".json")))

        sourceUpdate(x$td)
        expect_identical(sourceSchema(x$td)$version, 3L)
        expect_length(.gdsrc_registry_versions(x$td), 1L)
        expect_setequal(x$src[x$a, "giottosave"],
            c(.snap_tag(x$td, "qc"), "late.abcdefgh"))
        expect_length(list.files(d, pattern = "^[^e].*\\.json$"), 0L)
        # giottodir.json keeps only the stamp that turns older versions away
        stub <- .json_read(.gdsrc_json_path(x$td))
        expect_identical(stub$header$versions$manifest, 3L)
        expect_length(stub$content, 0L)
    })

    it("finishes a v2 -> v3 upgrade interrupted before the first version", {
        x <- .v1_project()
        on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
        .gdsrc_upgrade_1_to_2(x$td)
        manifest <- .gdsrc_legacy_read(x$td)$content
        .gdsrc_legacy_write(x$td, manifest, version = 3L) # stamped, not moved
        expect_identical(sourceSchema(x$td)$version, 2L)
        expect_error(sourcePrune(x$src), "sourceUpdate")

        sourceUpdate(x$td)
        expect_identical(sourceSchema(x$td)$version, 3L)
        expect_identical(x$src[x$b, "giottosave"], .snap_tag(x$td, "old"))
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
