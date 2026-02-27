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
    unlink(save_path, force = TRUE)
    invisible(TRUE)
})