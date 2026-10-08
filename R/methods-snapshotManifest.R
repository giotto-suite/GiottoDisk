#' @include AllGenerics.R
#' @include class-gsource.R
NULL

# Snapshot sidecars: a description of a saved gobject, kept beside it as
# text rather than inside it. Payload untouched, provenance greppable.
#
#   giottosave/<name>.rds             the object
#   giottosave/<name>.manifest.json   what it is    (GiottoClass::objManifest)
#   giottosave/<name>.history.ndjson  why it is     (one operation per line)
#
# The point is reading a snapshot's contents without loading it. Both files
# are derived from the object, so failing to write one loses a convenience,
# never data.

.gdsrc_manifest_path <- function(p, name) {
    file.path(.gdsrc_giottosave_dir(p), paste0(name, ".manifest.json"))
}

.gdsrc_history_path <- function(p, name) {
    file.path(.gdsrc_giottosave_dir(p), paste0(name, ".history.ndjson"))
}

# Written by snapshotSave() once the snapshot itself has landed. Uses the
# same temp-then-rename write as the artifacts manifest, so a reader never
# sees a half-written sidecar.
.ss_gdsrc_write_sidecars <- function(x, p, name, verbose = NULL) {
    tryCatch(
        {
            .atomic_write(
                writer = function(f) {
                    GiottoClass::objManifest_json(x, file = f, level = "full")
                },
                file = .gdsrc_manifest_path(p, name)
            )
            .atomic_write(
                writer = function(f) {
                    GiottoClass::objHistory_ndjson(x, file = f)
                },
                file = .gdsrc_history_path(p, name)
            )
            invisible(TRUE)
        },
        error = function(e) {
            warning("[snapshotSave] sidecars not written: ",
                conditionMessage(e),
                call. = FALSE
            )
            invisible(FALSE)
        }
    )
}

# most recently modified snapshot, used when no name is given
.gdsrc_latest_gsavename <- function(p) {
    gsdir <- .gdsrc_giottosave_dir(p)
    names <- .gdsrc_detect_gsavename(gsdir)
    if (length(names) == 0L) {
        stop("[GiottoDisk] no snapshots in this project directory\n",
            call. = FALSE
        )
    }
    files <- list.files(gsdir,
        pattern = "\\.rds$|\\.qs$", full.names = TRUE
    )
    names[[which.max(file.info(files)$mtime)]]
}

#' @name snapshotManifest
#' @title Read a Giotto Snapshot Manifest
#' @description
#' Read the manifest written beside a snapshot: a machine-readable inventory
#' of what the saved `giotto` object contains. Reads the sidecar text, so a
#' multi-gigabyte snapshot can be inspected without loading it.
#'
#' See [GiottoClass::objManifest()] for the structure, and
#' [snapshotHistory()] for the operation log written alongside it.
#' @param src `gsource` object or `character` filepath to project dir if
#'   `gDirSource` controlled.
#' @param name `character` (optional) name of the snapshot. If NULL, the most
#'   recently written one is used.
#' @param ... additional params to pass (none implemented)
#' @returns list of class `gmanifest`
NULL

#' @rdname snapshotManifest
#' @export
setMethod("snapshotManifest", signature("character"), function(
        src, name = NULL, ...) {
    snapshotManifest(sourceCreate(src, type = "gDirSource"), name = name, ...)
})

#' @rdname snapshotManifest
#' @export
setMethod("snapshotManifest", signature("gDirSource"), function(
        src, name = NULL, ...) {
    p <- src@path
    name <- name %null% .gdsrc_latest_gsavename(p)
    f <- .gdsrc_manifest_path(p, name)
    if (!file.exists(f)) {
        stop(sprintf(
            "[snapshotManifest] no manifest for snapshot '%s'.\n%s",
            name,
            "Snapshots written before manifests existed have none; re-save to create one."
        ), call. = FALSE)
    }
    out <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    class(out) <- c("gmanifest", "list")
    out
})

#' @name snapshotHistory
#' @title Read a Giotto Snapshot History
#' @description
#' Read the operation log written beside a snapshot: one JSON record per
#' operation, in the order they ran. Provenance only - object state is read
#' from [snapshotManifest()], never reconstructed by replaying this.
#' @inheritParams snapshotManifest
#' @returns list of records
NULL

#' @rdname snapshotHistory
#' @export
setMethod("snapshotHistory", signature("character"), function(
        src, name = NULL, ...) {
    snapshotHistory(sourceCreate(src, type = "gDirSource"), name = name, ...)
})

#' @rdname snapshotHistory
#' @export
setMethod("snapshotHistory", signature("gDirSource"), function(
        src, name = NULL, ...) {
    p <- src@path
    name <- name %null% .gdsrc_latest_gsavename(p)
    f <- .gdsrc_history_path(p, name)
    if (!file.exists(f)) {
        stop(sprintf(
            "[snapshotHistory] no history for snapshot '%s'\n", name
        ), call. = FALSE)
    }
    lines <- readLines(f, warn = FALSE)
    lapply(lines[nzchar(lines)], jsonlite::fromJSON, simplifyVector = FALSE)
})
