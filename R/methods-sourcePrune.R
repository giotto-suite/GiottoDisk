#' @name sourcePrune
#' @title Prune a project source
#' @description
#' On-disk processing can generate many intermediate
#' artifacts. Manually call this function in order to
#' drop artifacts that are not tagged with a
#' giottosave that depends on them (the giottosave .rds
#' file must also exist).
#'
#' What a deleted snapshot reached is kept for a grace period after it was
#' deleted, 24 hours by default, so that another session still using it does
#' not lose its data to a prune. This covers a save that replaces a snapshot
#' of the same name, which deletes the older one. Artifacts no snapshot ever
#' reached get no grace period: save work as a snapshot to keep it.
#' @param src `gsource` object. Used to provide default save locations
#' and save formats for a managed Giotto backend.
#' @param grace `TRUE` (default) to keep what snapshots deleted within the
#'   grace period reached; the period is `getOption("giottodisk.grace")`, in
#'   seconds or a `difftime`, 24 hours if unset. `FALSE` or `0` for none, or
#'   a number of seconds or a `difftime` for this call.
#' @param dry_run `logical`. If `TRUE`, report what would be pruned and change
#'   nothing.
#' @param ... additional params to pass (none implemented)
#' @returns `data.frame`, invisibly unless `dry_run`, one row per artifact
#'   no live snapshot protects: `uid`, `store` (its store type), `created`,
#'   `bytes`, and `kept_by`, the tag of the deleted snapshot whose grace
#'   period kept it, or `NA` for an artifact that was (or would be) removed.
#'   Artifacts removed first.
#' @examples
#' src <- gDirSource(file.path(tempdir(), "prune_example"))
#' sourcePrune(src, dry_run = TRUE)
#' @export
setMethod("sourcePrune", signature("gDirSource"),
    function(src, grace = TRUE, dry_run = FALSE, ...) {
        checkmate::assert_flag(dry_run)
        .gdsrc_artifact_prune(src@path, grace = .gdsrc_grace_secs(grace),
            dry_run = dry_run)
    })

#' @rdname sourcePrune
#' @export
setMethod("sourcePrune", signature("giotto"), function(src, ...) {
    if (is.null(src@source)) return(invisible())
    sourcePrune(src@source, ...)
})
