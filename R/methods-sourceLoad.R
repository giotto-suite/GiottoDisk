#' @name sourceLoad
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

#' @rdname sourceLoad
#' @export
setMethod("sourceLoad", signature("character"), function(src, ...) {
    if (!dir.exists(src)) {
        stop("[GiottoDisk] not an existing project directory\n", call. = FALSE)
    }
    gsrc <- sourceCreate(src, type = "gDirSource")
    sourceLoad(gsrc, ...)
})

#' @rdname sourceLoad
#' @export
setMethod("sourceLoad", signature("gDirSource"), function(src,
    name = NULL,
    load_params = list(),
    verbose = NULL,
    ...) {
    checkmate::assert_character(name, null.ok = TRUE)
    checkmate::assert_list(load_params)
    p <- src@path
    snaps_dir <- .gdsrc_giottosave_dir(p)
    
    # get path to gobject snapshot
    if (is.null(name)) { # find most recent snapshot if name = NULL
        existing_snaps <- list.files(snaps_dir, 
            pattern = "\\.rds|\\.qs", 
            full.names = TRUE,
            recursive = FALSE
        )
        modtimes <- file.info(existing_snaps)$mtime
        snap_path <- existing_snaps[which(modtimes == max(modtimes))][1L]
    } else {
        snap_path <- list.files(snaps_dir,
            pattern = paste0(name, "\\."),
            full.names = TRUE,
            recursive = FALSE
        )[1L]
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

