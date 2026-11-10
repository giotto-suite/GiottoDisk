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

# save description ####

# giottosave/<tag>/snapshot.json describes the save rather than the object:
# who made it, from which save, whether to keep it (adr/0024). Fixed when the
# save is written: schema, name, tag, created, creator, session, parent.
# Changeable with snapshotMetadata(): keep, description, labels, extra.

.GDSRC_SNAPSHOT_INFO <- "snapshot.json"
.GDSRC_INFO_SCHEMA <- 1L
.GDSRC_INFO_TIME <- "%Y-%m-%dT%H:%M:%OS3Z"

# A new save's description, as snapshot.json holds it
.ss_gdsrc_info <- function(name, tag, keep = TRUE, creator = "user",
    session = NULL, parent = NULL, description = NULL, labels = NULL,
    extra = NULL) {
    checkmate::assert_flag(keep)
    checkmate::assert_string(creator, min.chars = 1L)
    checkmate::assert_string(session, null.ok = TRUE)
    checkmate::assert_string(parent, null.ok = TRUE)
    checkmate::assert_string(description, null.ok = TRUE)
    checkmate::assert_character(labels, any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_list(extra, names = "unique", null.ok = TRUE)
    .gdsrc_info_fill(list(name = name, creator = creator, session = session,
        parent = parent, keep = keep, description = description,
        labels = labels, extra = extra
    ), tag = tag, created = format(Sys.time(), .GDSRC_INFO_TIME, tz = "UTC"))
}

# A save's description with every field present, in the form snapshot.json
# holds; a field `info` lacks takes its default. `created` is used only when
# `info` has none.
.gdsrc_info_fill <- function(info, tag, created) {
    scalar <- function(x, default) {
        x <- unlist(x)
        if (length(x) == 0L || is.na(x[[1L]])) default else x[[1L]]
    }
    extra <- info$extra
    list(
        schema = .GDSRC_INFO_SCHEMA,
        name = scalar(info$name, sub(.GDSRC_SNAPSHOT_DIR, "\\1", tag)),
        tag = tag,
        created = scalar(info$created, created),
        creator = scalar(info$creator, NA_character_),
        session = scalar(info$session, NA_character_),
        parent = scalar(info$parent, NA_character_),
        keep = as.logical(scalar(info$keep, TRUE)),
        description = scalar(info$description, ""),
        labels = as.character(unlist(info$labels)),
        extra = if (is.list(extra) && length(extra)) extra else
            structure(list(), names = character())
    )
}

.ss_gdsrc_write_info <- function(info, dir) {
    .json_atomic_write(info, file.path(dir, .GDSRC_SNAPSHOT_INFO))
}

# The snapshot.json of the save in folder `dir` (tag = folder name), filled.
# A save written before snapshot.json existed reads as kept, by an unknown
# creator, created when its object was written.
.gdsrc_stored_info <- function(dir) {
    f <- file.path(dir, .GDSRC_SNAPSHOT_INFO)
    info <- if (file.exists(f)) .json_read(f) else list()
    .gdsrc_info_fill(info, tag = basename(dir), created = format(
        file.info(c(.gdsrc_snapshot_objfiles(dir), dir)[[1L]])$mtime,
        .GDSRC_INFO_TIME, tz = "UTC"
    ))
}

# The description of the save in folder `dir`, as snapshotInfo() lists it
.gdsrc_read_info <- function(dir) {
    info <- .gdsrc_stored_info(dir)
    info$schema <- NULL
    info$created <- as.POSIXct(info$created, format = "%Y-%m-%dT%H:%M:%OSZ",
        tz = "UTC")
    info$deleted <- !length(.gdsrc_snapshot_objfiles(dir))
    info
}

#' @name snapshotInfo
#' @title Describe the Snapshots of a Project
#' @description
#' List a project's snapshots with what each save recorded about itself,
#' without loading any of them: who made it ("creator", such as `"user"` or
#' the package that saved it), the session it came from, the snapshot it was
#' made from (`parent`), whether it is to be kept, a description and labels.
#' Set those when saving (see [snapshotSave()]) or later with
#' [snapshotMetadata()].
#'
#' A deleted snapshot keeps its folder and description: only its object is
#' removed, which is what releases its artifacts to [sourcePrune()]. List
#' those too with `deleted = TRUE`.
#' @param src `gsource` object or `character` filepath to a project
#'   directory.
#' @param name `character` (optional) snapshot name to describe; all of them
#'   if `NULL`.
#' @param deleted `logical`. Whether to include deleted snapshots.
#' @param ... additional params to pass (none implemented)
#' @returns `data.frame`, one row per snapshot, oldest first: `name`, `tag`,
#'   `created`, `creator`, `session`, `parent` (a tag), `keep`,
#'   `description`, `labels` (list column), `extra` (list column),
#'   `deleted`.
#' @examples
#' src <- gDirSource(file.path(tempdir(), "info_example"))
#' snapshotInfo(src)
#' @export
setMethod("snapshotInfo", signature("character"), function(src, ...) {
    if (!dir.exists(src)) {
        stop("[snapshotInfo] not an existing project directory\n",
            call. = FALSE)
    }
    snapshotInfo(sourceCreate(src, type = "gDirSource"), ...)
})

#' @rdname snapshotInfo
#' @export
setMethod("snapshotInfo", signature("gDirSource"), function(src, name = NULL,
    deleted = FALSE, ...) {
    checkmate::assert_string(name, null.ok = TRUE)
    checkmate::assert_flag(deleted)
    gsdir <- .gdsrc_giottosave_dir(src@path)
    tags <- .gdsrc_snapshot_dirs(gsdir)
    infos <- lapply(file.path(gsdir, tags), .gdsrc_read_info)
    keep <- vapply(infos, function(i) {
        (deleted || !i$deleted) && (is.null(name) || identical(i$name, name))
    }, NA)
    infos <- infos[keep]
    out <- data.frame(
        name = vapply(infos, `[[`, "", "name"),
        tag = vapply(infos, `[[`, "", "tag"),
        created = .POSIXct(vapply(infos, function(i) as.numeric(i$created), 0),
            tz = "UTC"),
        creator = vapply(infos, `[[`, "", "creator"),
        session = vapply(infos, `[[`, "", "session"),
        parent = vapply(infos, `[[`, "", "parent"),
        keep = vapply(infos, `[[`, NA, "keep"),
        description = vapply(infos, `[[`, "", "description"),
        deleted = vapply(infos, `[[`, NA, "deleted"),
        stringsAsFactors = FALSE
    )
    out$labels <- lapply(infos, `[[`, "labels")
    out$extra <- lapply(infos, `[[`, "extra")
    out <- out[order(out$created), c("name", "tag", "created", "creator",
        "session", "parent", "keep", "description", "labels", "extra",
        "deleted")]
    rownames(out) <- NULL
    out
})

#' @name snapshotMetadata
#' @title Change What a Snapshot Records About Itself
#' @description
#' Change a snapshot's `keep` flag, description, labels or `extra` fields,
#' as listed by [snapshotInfo()]. Arguments left `NULL` are unchanged. Who
#' made the snapshot, when, in which session and from which snapshot are
#' fixed when it is saved. The snapshot's object and the artifacts it
#' protects are not touched.
#' @inheritParams snapshotInfo
#' @param name `character` name of the snapshot.
#' @param keep `logical`. Whether the snapshot is to be kept long-term.
#' @param description `character` free text.
#' @param labels `character` labels; replace the current ones.
#' @param extra named `list` of fields kept for the package that saved the
#'   snapshot, merged into the current ones (a `NULL` element removes one).
#' @returns the snapshot's row of [snapshotInfo()], invisibly
#' @export
setMethod("snapshotMetadata", signature("character", "character"),
    function(src, name, ...) {
    if (!dir.exists(src)) {
        stop("[snapshotMetadata] not an existing project directory\n",
            call. = FALSE)
    }
    snapshotMetadata(sourceCreate(src, type = "gDirSource"), name, ...)
})

#' @rdname snapshotMetadata
#' @export
setMethod("snapshotMetadata", signature("gDirSource", "character"),
    function(src, name, keep = NULL, description = NULL, labels = NULL,
    extra = NULL, ...) {
    checkmate::assert_flag(keep, null.ok = TRUE)
    checkmate::assert_string(description, null.ok = TRUE)
    checkmate::assert_character(labels, any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_list(extra, names = "unique", null.ok = TRUE)
    p <- src@path
    .gdsrc_check_schema(p, site = "[snapshotMetadata]")
    snap <- .gdsrc_find_snapshot(.gdsrc_giottosave_dir(p), name = name)
    if (is.null(snap)) {
        stop(sprintf("[snapshotMetadata] no snapshot named '%s'\n", name),
            call. = FALSE)
    }
    dir <- dirname(snap$file)
    info <- .gdsrc_stored_info(dir)
    if (!is.null(keep)) info$keep <- keep
    if (!is.null(description)) info$description <- description
    if (!is.null(labels)) info$labels <- labels
    if (!is.null(extra)) {
        merged <- utils::modifyList(info$extra, extra)
        info$extra <- if (length(merged)) merged else
            structure(list(), names = character())
    }
    .ss_gdsrc_write_info(info, dir)
    out <- snapshotInfo(src, name = name)
    invisible(out[out$tag == snap$tag, , drop = FALSE])
})
