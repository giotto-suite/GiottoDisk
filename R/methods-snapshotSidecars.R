#' @include AllGenerics.R
#' @include class-gsource.R
NULL

# Snapshot sidecars describe a saved gobject beside it, as text:
#
#   giottosave/<name>.<id>/gobject.<ext>    the object
#   giottosave/<name>.<id>/manifest.json    what it is   (objManifest)
#   giottosave/<name>.<id>/history.ndjson   why it is    (objHistory)
#
# The point is reading a snapshot's contents without loading it. GiottoClass
# owns both formats and reads them back from a folder; this file only writes
# them at save time and finds the folder for a snapshot name.

.GDSRC_SNAPSHOT_MANIFEST <- "manifest.json"
.GDSRC_SNAPSHOT_HISTORY <- "history.ndjson"

# Written into the staging folder, so they land with the object in one rename.
# Derived from the object, so a failed write costs a convenience, never data.
# giottoMulti has no manifest yet.
.ss_gdsrc_write_sidecars <- function(x, dir) {
    if (!inherits(x, "giotto")) return(invisible(FALSE))
    tryCatch(
        {
            objManifest(x, file = file.path(dir, .GDSRC_SNAPSHOT_MANIFEST))
            objHistory(x, file = file.path(dir, .GDSRC_SNAPSHOT_HISTORY))
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

#' @name objManifest-gDirSource
#' @title Read a snapshot's manifest or history
#' @description
#' [GiottoClass::objManifest()] and [GiottoClass::objHistory()] for a project
#' directory: read what a `giotto` snapshot holds, and how it was produced,
#' from the files [snapshotSave()] writes beside it, without loading the
#' snapshot. The result is the same as calling them on the object when it
#' was saved, so a snapshot can go to [GiottoClass::manifestDiff()].
#'
#' A project directory path works too: GiottoClass detects it with
#' [resolveSource()] and calls these methods, so `objManifest("my_project")`
#' is the same as `objManifest(gDirSource("my_project"))`.
#'
#' `giottoMulti` snapshots have no manifest or history yet.
#' @param x `gDirSource`
#' @param name `character` (optional) name of the snapshot. If NULL, the most
#'   recently written one is used.
#' @param output,file,pretty,summarized see [GiottoClass::objManifest()] and
#'   [GiottoClass::objHistory()]
#' @param ... additional params to pass (none implemented)
#' @returns as [GiottoClass::objManifest()] / [GiottoClass::objHistory()]
#' @examples
#' \dontrun{
#' src <- gDirSource("my_project")
#' objManifest(src)            # most recent snapshot
#' objHistory(src, "clustered")
#' objManifest("my_project", name = "clustered")
#' manifestDiff(objManifest(src, "normalized"), objManifest(src, "clustered"))
#' }
NULL

#' @rdname objManifest-gDirSource
#' @importFrom GiottoClass objManifest
#' @export
setMethod("objManifest", signature("gDirSource"), function(
        x, name = NULL, output = c("gmanifest", "json"), file = NULL,
        pretty = TRUE, ...) {
    f <- .gdsrc_sidecar_path(x@path, name, .GDSRC_SNAPSHOT_MANIFEST,
        site = "[objManifest]"
    )
    objManifest(f, output = output, file = file, pretty = pretty)
})

#' @rdname objManifest-gDirSource
#' @importFrom GiottoClass objHistory
#' @export
setMethod("objHistory", signature("gDirSource"), function(
        x, name = NULL, summarized = FALSE, output = c("ghistory", "ndjson"),
        file = NULL, ...) {
    f <- .gdsrc_sidecar_path(x@path, name, .GDSRC_SNAPSHOT_HISTORY,
        site = "[objHistory]"
    )
    objHistory(f, summarized = summarized, output = output, file = file)
})
