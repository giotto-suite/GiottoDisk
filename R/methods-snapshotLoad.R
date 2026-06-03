#' @name snapshotLoad
#' @title Load Giotto Snapshot
#' @description
#' Load a gobject snapshot. If no name is provided, the most recent
#' snapshot is automatically loaded.
#' @param src `gsource` object or `character` filepath to project dir if
#'   `gDirSource` controlled.
#' @param name `character` (optional) name of specific snapshot to load.
#'   if NULL, the most recent is selected
#' @param load_params `list` additional parameters (if any) for loading
#'   or reading giotto object
#' @param verbose verbosity
#' @param ... additional params to pass (none implemented)
#' @returns `giotto` object (possibly with additional downstream load
#'   steps needed)
NULL

#' @rdname snapshotLoad
#' @export
setMethod("snapshotLoad", signature("character"), function(src, ...) {
    if (!dir.exists(src)) {
        stop("[snapshotLoad] not an existing project directory\n", call. = FALSE)
    }
    gsrc <- sourceCreate(src, type = "gDirSource")
    snapshotLoad(gsrc, ...)
})

#' @rdname snapshotLoad
#' @export
setMethod("snapshotLoad", signature("gDirSource"), function(src,
    name = NULL,
    load_params = list(),
    verbose = NULL,
    ...) {
    checkmate::assert_character(name, null.ok = TRUE)
    checkmate::assert_list(load_params)
    p <- src@path
    snaps_dir <- .gdsrc_giottosave_dir(p)

    existing_snaps <- list.files(snaps_dir,
        pattern = "\\.(rds|qs)$",
        full.names = TRUE,
        recursive = FALSE
    )
    if (length(existing_snaps) == 0L) {
        stop("[snapshotLoad] no snapshots found in '", snaps_dir, "'",
            call. = FALSE)
    }

    if (is.null(name)) { # find most recent snapshot if name = NULL
        modtimes <- file.info(existing_snaps)$mtime
        snap_path <- existing_snaps[which.max(modtimes)]
    } else {
        hits <- existing_snaps[grepl(paste0("^", name, "\\."),
            basename(existing_snaps))]
        if (length(hits) == 0L) {
            stop(.snapshot_name_not_found_msg(name, existing_snaps),
                call. = FALSE)
        }
        snap_path <- hits[1L]
    }

    gobject <- .load_serialized(snap_path, load_params = load_params)
    return(gobject)
    
    # to be completed by GiottoClass::loadGiotto
    # - image reconnection
    # - python path handling
    # - data.table overallocation
    # - gobject initialization
})

# internals ####

# Build a "no snapshot named X" error message, with fuzzy suggestions when
# `name` is close to one of the available snapshot names. Uses base::adist
# (edit distance) so no new dep.
.snapshot_name_not_found_msg <- function(name, existing_snaps) {
    available <- tools::file_path_sans_ext(basename(existing_snaps))
    head_msg <- paste0("[snapshotLoad] no snapshot named '", name, "'.")
    if (length(available) == 0L) return(head_msg)

    # Suggest by: (1) substring containment in either direction, then
    # (2) close edit distance (<= 2, or <= 25% of name length).
    contains <- available[grepl(name, available, fixed = TRUE) |
        vapply(available, grepl, logical(1L), x = name, fixed = TRUE)]
    d <- as.integer(utils::adist(name, available))
    threshold <- max(2L, ceiling(nchar(name) * 0.25))
    close <- available[!is.na(d) & d <= threshold]
    sugg <- unique(c(contains, close))

    if (length(sugg) > 0L) {
        sugg <- head(sugg, 5L)
        return(paste0(head_msg, " Did you mean: ",
            paste(shQuote(sugg), collapse = ", "), "?"))
    }
    paste0(head_msg, " Available: ",
        paste(shQuote(head(available, 10L)), collapse = ", "),
        if (length(available) > 10L) ", ..." else "")
}

.load_serialized <- function(path, load_params = list()) {
    fext <- tail(file_extension(path), 1L)
    fext <- match.arg(tolower(fext), choices = c("rds", "qs"))
  
    read_fun <- switch(fext,
        "rds" = readRDS,
        "qs" = {
            package_check("qs", repository = "CRAN")
            read_fun <- get("qread", asNamespace("qs"))
        }
    )
    
    do.call(read_fun, args = c(list(file = path), load_params))
}

