#' @include AllGenerics.R
#' @include class-gsource.R
NULL

#' @name sourceInfo
#' @title Summarize a Project's Disk Use
#' @description
#' How much space a project uses and where, how many snapshots it has, and
#' what [sourcePrune()] would free, without loading anything. Sizes are
#' computed from the files when asked, so this reads every folder of the
#' project. Per snapshot, see [snapshotInfo()] with `sizes = TRUE`.
#' @param src `gDirSource` object or `character` filepath to a project
#'   directory.
#' @param grace as in [sourcePrune()]: which deleted snapshots still keep
#'   what they reached.
#' @param ... additional params to pass (none implemented)
#' @returns A `gdsrcInfo` list:
#' - `path`
#' - `bytes`: `vault` (all artifacts), `images` (image artifacts, part of
#'   `vault`) and `snapshots` (snapshot folders: object files and sidecars)
#' - `snapshots`: counts of `live`, `kept` and `unkept` (live ones by their
#'   `keep` flag) and `deleted` snapshots
#' - `artifacts`: the number of artifacts in the vault
#' - `prunable`: the artifacts [sourcePrune()] would remove now (`n`,
#'   `bytes`), and the newest one's modification time (`newest`), so recent
#'   work no snapshot holds stands out
#' - `graced`: the artifacts kept only by snapshots deleted within the grace
#'   period (`n`, `bytes`)
#' @examples
#' src <- gDirSource(file.path(tempdir(), "info_example"))
#' sourceInfo(src)
#' @export
setMethod("sourceInfo", signature("character"), function(src, ...) {
    if (!dir.exists(src)) {
        stop("[sourceInfo] not an existing project directory\n",
            call. = FALSE)
    }
    sourceInfo(sourceCreate(src, type = "gDirSource"), ...)
})

#' @rdname sourceInfo
#' @export
setMethod("sourceInfo", signature("gDirSource"), function(src, grace = TRUE,
    ...) {
    p <- src@path
    .gdsrc_check_schema(p, site = "[sourceInfo]")
    content <- .gdsrc_current_content(p)
    bytes <- .gdsrc_vault_bytes(p)
    store <- vapply(names(bytes), function(u) {
        s <- unlist(content[[u]]$store)
        if (length(s)) as.character(s[[1L]]) else NA_character_
    }, "", USE.NAMES = FALSE)
    gsdir <- .gdsrc_giottosave_dir(p)
    snaps <- snapshotInfo(src, deleted = TRUE)
    live <- !snaps$deleted

    plan <- .gdsrc_prune_plan(p, content, .gdsrc_grace_secs(grace))
    gone <- names(plan$kept_by)[is.na(plan$kept_by)]
    kept <- names(plan$kept_by)[!is.na(plan$kept_by)]
    newest <- if (length(gone)) {
        max(file.info(.gdsrc_artifact_dir(p, gone))$mtime)
    } else {
        .POSIXct(NA_real_)
    }
    structure(list(
        path = p,
        bytes = c(
            vault = sum(bytes),
            images = sum(bytes[store %in% "IMAGE"]),
            snapshots = sum(vapply(file.path(gsdir, snaps$tag),
                .gdsrc_dir_bytes, 0))
        ),
        snapshots = c(live = sum(live), kept = sum(live & snaps$keep),
            unkept = sum(live & !snaps$keep), deleted = sum(!live)),
        artifacts = length(bytes),
        prunable = list(n = length(gone), bytes = sum(bytes[gone]),
            newest = newest),
        graced = list(n = length(kept), bytes = sum(bytes[kept]))
    ), class = "gdsrcInfo")
})

#' @export
print.gdsrcInfo <- function(x, ...) {
    size <- function(b) format(structure(b, class = "object_size"),
        units = "auto", standard = "SI")
    s <- x$snapshots
    cat("GiottoDisk project\n")
    print_list(list(
        path = x$path,
        vault = sprintf("%s in %d artifacts (images %s)", size(x$bytes[["vault"]]),
            x$artifacts, size(x$bytes[["images"]])),
        snapshots = sprintf("%d live (%d kept, %d unkept), %d deleted; %s",
            s[["live"]], s[["kept"]], s[["unkept"]], s[["deleted"]],
            size(x$bytes[["snapshots"]])),
        prunable = if (x$prunable$n) {
            sprintf("%d artifacts, %s; newest %s", x$prunable$n,
                size(x$prunable$bytes), format(x$prunable$newest, "%Y-%m-%d %H:%M"))
        } else {
            "none"
        },
        graced = if (x$graced$n) {
            sprintf("%d artifacts, %s, kept for recently deleted snapshots",
                x$graced$n, size(x$graced$bytes))
        } else {
            "none"
        }
    ), pre = "  ")
    invisible(x)
}
