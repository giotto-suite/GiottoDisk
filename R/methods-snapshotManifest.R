#' @include AllGenerics.R
#' @include class-gsource.R
NULL

# Snapshot sidecars describe a saved gobject beside it, as text:
#
#   giottosave/<name>.<id>/gobject.<ext>    the object
#   giottosave/<name>.<id>/manifest.json    what it is   (GiottoClass::objManifest)
#   giottosave/<name>.<id>/history.ndjson   why it is    (one operation per line)
#
# The point is reading a snapshot's contents without loading it.

.GDSRC_SNAPSHOT_MANIFEST <- "manifest.json"
.GDSRC_SNAPSHOT_HISTORY <- "history.ndjson"

# Written into the staging folder, so they land with the object in one rename.
# Derived from the object, so a failed write costs a convenience, never data.
.ss_gdsrc_write_sidecars <- function(x, dir) {
    if (!inherits(x, "giotto")) return(invisible(FALSE))
    tryCatch(
        {
            GiottoClass::objManifest_json(x,
                file = file.path(dir, .GDSRC_SNAPSHOT_MANIFEST), level = "full"
            )
            GiottoClass::objHistory_ndjson(x,
                file = file.path(dir, .GDSRC_SNAPSHOT_HISTORY)
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

# path of a sidecar of the snapshot `name` (the most recent if NULL)
.gdsrc_sidecar_path <- function(p, name, file, site) {
    .gdsrc_check_schema(p, site = site)
    snap <- .gdsrc_find_snapshot(.gdsrc_giottosave_dir(p), name = name)
    if (is.null(snap)) {
        stop(site, if (is.null(name)) " no snapshots found\n" else
            sprintf(" no snapshot named '%s'\n", name), call. = FALSE)
    }
    f <- file.path(dirname(snap$file), file)
    if (!file.exists(f)) {
        stop(sprintf(
            "%s no %s for snapshot '%s'.\n%s", site, file, snap$name,
            "Re-save the snapshot to create one."
        ), call. = FALSE)
    }
    f
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
    f <- .gdsrc_sidecar_path(src@path, name, .GDSRC_SNAPSHOT_MANIFEST,
        site = "[snapshotManifest]"
    )
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
    f <- .gdsrc_sidecar_path(src@path, name, .GDSRC_SNAPSHOT_HISTORY,
        site = "[snapshotHistory]"
    )
    lines <- readLines(f, warn = FALSE)
    lapply(lines[nzchar(lines)], jsonlite::fromJSON, simplifyVector = FALSE)
})
