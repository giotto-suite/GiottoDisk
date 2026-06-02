#' @name snapshotDelete
#' @title Delete a Giotto Snapshot
#' @description
#' Delete a Giotto snapshot. By extension, also removes
#' the protection of the associated giottosave tag for
#' artifact pruning.
#' @param src `gsource` object or `character` filepath to project dir if
#'   `gDirSource` controlled.
#' @param name `character` (optional) name of specific snapshot to delete.
#'   if NULL, the most recent is selected
#' @param ... additional params to pass (none implemented)
#' @returns `TRUE` invisibly if succeeds
NULL

#' @rdname snapshotDelete
#' @export
setMethod("snapshotDelete", signature("character", "character"), function(src, name, ...) {
    if (!dir.exists(src)) {
        stop("[snapshotDelete] not an existing project directory\n", call. = FALSE)
    }
    gsrc <- sourceCreate(src, type = "gDirSource")
    snapshotDelete(gsrc, name, ...)
})

#' @rdname snapshotDelete
#' @export
setMethod("snapshotDelete", signature("gDirSource", "character"), function(src, name, ...) {
    p <- src@path
    gsdir <- .gdsrc_giottosave_dir(p)
    save_path <- list.files(gsdir,
        pattern = paste0("^", name, "\\."),
        full.names = TRUE,
        recursive = FALSE
    )
    if (length(save_path) == 0L) {
        stop(
            sprintf("[snapshotDelete] no snapshot found with name '%s'\n", name),
            call. = FALSE
        )
    }

    # Cascade child deletes for giottoMulti snapshots. The multi save
    # produces a snapshot per child under "<name>_<child_name>" in each
    # child's own vault — without cascade, deleting the multi here would
    # orphan all those child .rds files on disk and leave their
    # protection tags pinning artifacts indefinitely.
    #
    # Loading a backed snapshot is cheap (gobjects hold file handles,
    # not data). Best-effort: a child whose source is unreachable or
    # whose snapshot was already deleted produces a warning, not a
    # hard error — the multi delete always proceeds.
    snap <- tryCatch(
        .load_serialized(save_path[[1L]]),
        error = function(e) NULL
    )
    if (inherits(snap, "giottoMulti")) {
        for (child_name in names(snap@objects)) {
            child <- snap@objects[[child_name]]
            if (is.null(child@source)) next
            child_snap_name <- paste0(name, "_", child_name)
            tryCatch(
                snapshotDelete(child@source, child_snap_name),
                error = function(e) warning(sprintf(
                    "[snapshotDelete] child snapshot '%s' delete failed: %s",
                    child_snap_name, conditionMessage(e)
                ), call. = FALSE)
            )
        }
    }

    unlink(save_path, force = TRUE)
    invisible(TRUE)
})