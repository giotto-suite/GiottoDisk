
# classes ####

#' @name gsource
#' @title Giotto Sources
#' @slot path character. Filepath to use
#' @description
#' Object defining a data backend for a Giotto object. The `gsource`
#' class is VIRTUAL. Specific backend types should extend from this class.
#' @export
setClass("gsource",
    contains = "VIRTUAL",
    slots = list(
        path = "character"
    ),
    prototype = list(
        path = NA_character_
    )
)

#' @name gDirSource-class
#' @title Giotto Directory Source
#' @slot path character. Directory to use
#' @description
#' Extends [gsource]. Used to designate and manage a directory as a
#' giotto project directory containing on-disk data artifacts.
#' Contained artifacts are tracked via a manifest with editable
#' metadata.
#' 
#' Edits to the `gDirSource` manifest are atomic. Additionally, a WAL
#' pattern is used to ensure concurrency safety.
#' Edits are only made as loose files specific to single artifacts 
#' (and by extension usually single-process).
#' Loose edits are later consolidated into the main manifest on-read
#' (preferably only from the main process). {filelock} integration
#' is available for extra concurrency safety to avoid edge cases
#' where two processes attempt consolidation at the same time but see
#' different `_pending` states.
#' 
#' Data can be written to the backing directory using `sourceWrite()`.
#' This serializes the data to disk using the preferred store types
#' for the `gsource` class.
#' 
#' # Directory Structure
#' 
#' * `giottodir.json` - json file manifest of the contained artifacts.
#' * `giottodir.json.lock` - (optional) lockfile for {filelock} to control
#' writes to `giottodir.json` for extra concurrency safety
#' * `artifacts` - vault directory containing the actual data artifacts
#' which are within subdirectories named by their `uid`.
#' * `_pending` - directory containing pending edits to the manifest.
#' * `giottosave` - directory containing .RDS of giotto projects that
#' reference assets controlled by the `gDirSource`.
#' 
#' @export
setClass("gDirSource",
    contains = "gsource"
)

# constructor ####

#' @name gDirSource
#' @title Create a Giotto Directory Source
#' @description
#' Create a `gDirSource` object. See [gDirSource-class] for more information.
#' @param path character. Directory to use as Giotto project backend.
#' @param ... additional params to pass (none implemented)
#' @export
gDirSource <- function(path, ...) {
    new("gDirSource", path = path, ...)
}

#' @name sourceCreate
#' @title Create a Giotto Source
#' @description
#' Hub function for creating Giotto source objects.
#' @param path character. Directory to use
#' @param ... additional params to pass
#' @seealso [gDirSource()]
#' @export
sourceCreate <- function(type = "gsource", path, ...) {
    type <- match.arg(type, choices = c("gsource"))
    switch(type,
        "gsource" = gDirSource(path, ...)
    )
}

# initialize ####

setMethod("initialize", signature("gDirSource"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    p <- .Object@path <- normalizePath(.Object@path, mustWork = FALSE)
    if (is.na(p)) {
        stop("[gDirSource] 'path' should be a directory path for the giotto project.\n",
            call. = FALSE
        )
    }
    if (!file.exists(.gdsrc_json_path(p))) {
        .gdsrc_json_write(p, list()) # initialize with empty list
    }
    .Object
})

# methods ####

setMethod("show", signature("gDirSource"), function(object) {
    cat(sprintf("<%s>\n", class(object)))

    gsave_ids <- list.files(.gdsrc_giottosave_dir(object@path),
        full.names = FALSE, recursive = FALSE, pattern = "\\.rds"
    )
    cat("giottosaves:", length(gsave_ids), "\n")
    cat("artifacts:", length(object), "\n")
})

setMethod("length", signature("gDirSource"), function(x) {
    vd <- .gdsrc_vault_dir(x@path)
    if (!dir.exists(vd)) {
        return(0L)
    }
    uids <- list.dirs(vd, full.names = FALSE, recursive = FALSE)
    length(uids)
})

setMethod("[", signature("gDirSource", i = "character", j = "missing"), function(x, i, j) {
    data <- .gdsrc_json_read(x@path, consolidate = TRUE)
    data$content[[i]]
})

setMethod("[", signature("gDirSource", i = "character", j = "character"), function(x, i, j) {
    data <- .gdsrc_json_read(x@path, consolidate = TRUE)
    data$content[[i]][[j]]
})

setMethod("[", signature("gDirSource", i = "missing", j = "character"), function(x, i, j) {
    data <- .gdsrc_json_read(x@path, consolidate = TRUE)
    lapply(data$content, function(art) art[[j]])
})

setMethod("[<-", signature("gDirSource", i = "character", j = "missing", value = "ANY"), function(x, i, j, ..., value) {
    checkmate::assert_list(value)
    .gdsrc_json_edit(x@path, uid = i, x = value)
    x
})

setMethod("[<-", signature("gDirSource", i = "character", j = "character", value = "ANY"), function(x, i, j, ..., value) {
    checkmate::assert_atomic(value)
    entry <- list()
    entry[[j]] <- value
    .gdsrc_json_edit(x@path, uid = i, x = entry)
    x
})

# json internal tools ####

#' @name giotto_json
#' @title Giotto Directory JSON
#' @description
#' Internal utilities for editing the [gDirSource] managed
#' `giottodir.json`.
#' @keywords internal
NULL

## path setting ####

# p is giotto directory
.gdsrc_json_path <- function(p) {
    file.path(p, "giottodir.json")
}

.gdsrc_json_pending_dir <- function(p) {
    file.path(p, "_pending")
}

.gdsrc_vault_dir <- function(p) {
    file.path(p, "artifacts")
}

.gdsrc_artifact_dir <- function(p, uid) {
    vd <- .gdsrc_vault_dir(p)
    file.path(vd, uid)
}

.gdsrc_giottosave_dir <- function(p) {
    file.path(p, "giottosave")
}

## tools ####

#' @describeIn giotto_json Safe atomic writes of information to json
#' files so that there is no instant at which the data does not exist
#' on disk.
.json_atomic_write <- function(x, file, 
    temp_file = tempfile(tmpdir = dirname(file), fileext = ".tmp"), 
    cleanup = TRUE) {
    jsonlite::write_json(x, temp_file)

    on.exit({
        if (cleanup) {
            unlink(temp_file, force = TRUE)
        }
    }, add = TRUE)

    success <- file.rename(from = temp_file, to = file)
    if (!success) {
        stop("[GiottoDisk] Failed to rename ", temp_file, " to ", file)
    }
    cleanup <- FALSE # disable cleanup on success
    invisible(TRUE)
}

#' @describeIn giotto_json json reading with specific params to work
#' with GiottoDisk's usecase where values are generally atomic vectors
.json_read <- function(file) {
    jsonlite::fromJSON(file, simplifyVector = TRUE,
        simplifyMatrix = FALSE
    )
}

#' @describeIn giotto_json general header to apply at top of json
#' files from GiottoDisk. Stamps with information for:
#' 
#' * GiottoDisk version
#' * manifest version
.gdsrc_json_header <- function() {
    list(
        "description" = "[GiottoDisk] Artifacts Manifest",
        "versions" = list(
            "package" = as.character(packageVersion("GiottoDisk")),
            "manifest" = 1L
        )
    )
}

#' @describeIn giotto_json Write function for the main
#' `giottodir.json`.
.gdsrc_json_write <- function(p, x) { # write a giotto dir json
    checkmate::assert_list(x)
    if (!dir.exists(p)) { # first creation
        message("Setting up Giotto project directory at: \n", p)
        dir.create(path = p)
    }
  
    content <- list(
        "header" = .gdsrc_json_header(),
        "content" = x
    )
  
    .json_atomic_write(x = content, file = .gdsrc_json_path(p))
    invisible(TRUE)
}

#' @describeIn giotto_json Read function for the main
#' `giottodir.json`. This function does not check or consolidate
#' pending edits by default.
.gdsrc_json_read <- function(p, consolidate = FALSE) {
    if (isTRUE(consolidate)) {
        .gdsrc_json_consolidate(p)
    }
  
    etag <- "[gDirSource]" # tag source of error
    if (!dir.exists(p)) {
        stop(etag,
            " Path to Giotto project directory does not exist\n",
            call. = FALSE
        )
    }
    json_path <- .gdsrc_json_path(p)
    if (!file.exists(json_path)) {
        stop(sprintf("%s giottodir.json not found in:\n%s\n", etag, p),
                call. = FALSE)
    }
    .json_read(json_path)
}

#' @describeIn giotto_json Queue an edit to the manifest by
#' generating a pending edit that is written to disk as a
#' loose json. This can later be consolidated into the main
#' `giottodir.json`.
.gdsrc_json_edit <- function(p, uid, x) {
    checkmate::assert_character(uid)
    checkmate::assert_list(x)
    pending_dir <- .gdsrc_json_pending_dir(p)
    file <- file.path(pending_dir, sprintf("%s.json", uid))
  
    if (!dir.exists(pending_dir)) { # first creation
        dir.create(path = pending_dir)
    }
    
    if (file.exists(file)) {
        # warning to discourage multiple worker edits on pending entries
        warning("[GiottoDisk] modifying artifact metadata before ",
            "consolidation is not well supported")
        data <- .json_read(file)
      
        # !!! future version catching !!!
      
        data <- data$content
        x <- modifyList(data, x)
    }
  
    # package edit entry
    entry <- list(
        "header" = .gdsrc_json_header(),
        "content" = x
    )
  
    .json_atomic_write(entry, file)
    invisible(TRUE)
}

#' @describeIn giotto_json Add a tracked artifact to the directory. This adds
#' the tracking metadata for the object and places it in a location specified
#' by `.gdsrc_json_pending_dir()` to await consolidation into the main
#' `giottodir.json`.
#' @param store_type character. Type of file format storage (e.g. \"h5\")
#' @param uid character. Unique ID for artifact tracking
#' @param hash character. Hash of the in-memory store object to track changes
#' @param meta list (optional). Additional list of atomic object(s) that can
#' be attached as further metadata to the particular uid
#' @keywords internal
.gdsrc_json_add_artifact <- function(p, store_type, uid, hash, meta = NULL) {
    checkmate::assert_character(store_type)
    checkmate::assert_character(uid)
    checkmate::assert_character(hash)
    checkmate::assert_list(meta, null.ok = TRUE)

    content <- list(
        "time" = .timestamp(),
        "store" = store_type,
        "giottosave" = NA_character_, # tagged giotto save(s)
        "hash" = hash
    )
    content <- c(content, meta)
    .gdsrc_json_edit(p = p, uid = uid, x = content)
}

.gdsrc_json_tag_giottosave <- function(p, uid, giottosave) {
    checkmate::assert_character(p)
    checkmate::assert_character(uid)
    checkmate::assert_character(giottosave)
    content <- list(
        "giottosave" = giottosave
    )
    .gdsrc_json_edit(p = p, uid = uid, x = content)
}

#' @describeIn giotto_json Consolidate pending edits into central
#' `giottodir.json` manifest. Scans for pending edits then applies
#' the changes in a for loop on the manifest content before
#' writing back out.
.gdsrc_json_consolidate <- function(p) {
    pending_dir <- .gdsrc_json_pending_dir(p)
    scan_pending <- function() {
        list.files(pending_dir, full.names = TRUE, pattern = "\\.json")
    }
    edits <- scan_pending()
    if (!dir.exists(pending_dir) || 
        length(edits) == 0L) {
        #  no pending edits to giottodir.json
        return(invisible(TRUE))
    }
    
    # optional: locking for safer consolidate
    use_locking <- getOption("giottodisk.use_locking", TRUE)
    if (use_locking && requireNamespace("filelock", quietly = TRUE)) {
        lock <- filelock::lock(paste0(.gdsrc_json_path(p), ".lock"), 
            timeout = 10000)
        on.exit(filelock::unlock(lock), add = TRUE)
        # rescan after getting lock (in case of changes during timeout)
        edits <- scan_pending()
        if (length(edits) == 0L) {
            # early return if other process already consolidated
            return(invisible(TRUE))
        }
    } else if (use_locking) {
        # Optional: warn once per session
        if (!getOption("giottodisk.lock_warning_shown", FALSE)) {
            message("[GiottoDisk] Install 'filelock' for safer concurrent access")
            options(giottodisk.lock_warning_shown = TRUE)
        }
    }
    
    uids <- gsub("\\.json$", "", basename(edits)) # artifact uids
    names(edits) <- uids
  
    # append to giottodir.json
    json_data <- .gdsrc_json_read(p, consolidate = FALSE)
    manifest <- json_data$content
    
    for (artifact in uids) {
        edit <- .json_read(edits[[artifact]])
        
        # !!! future version catching !!!
        
        edit <- edit$content
      
        manifest[[artifact]] <-
            modifyList(manifest[[artifact]] %||% list(), edit)
    }
  
    data <- list(
        "header" = .gdsrc_json_header(),
        "content" = manifest
    )
  
    .json_atomic_write(data, .gdsrc_json_path(p))
    for (efile in edits) { # remove known edit files
        unlink(efile, recursive = TRUE, force = TRUE)
    }
    invisible(TRUE)
}

#' @title Coerce gDirSource
#' @name coerce_gdsrc
#' @description
#' Coerce an object to data.frame
#' @param x object to coerce
#' @export
as.list.gDirSource <- function(x) {
    data <- .gdsrc_json_read(x@path, consolidate = TRUE)
    data$content
}

#' @rdname coerce_gdsrc
#' @param ... additional params to pass
#' @export
as.data.frame.gDirSource <- function(x, ...) {
    data <- as.list.gDirSource(x)
    data.table::rbindlist(data, 
        fill = TRUE, idcol = "uid"
    )
}

## pruning ####

# protecting tags:
# - giottosave
# - depends
.gdsrc_artifact_prune <- function(p) {
    gsave_dir <- .gdsrc_giottosave_dir(p)
    gsave_ids <- list.files(gsave_dir,
        pattern = "\\.rds", recursive = FALSE, full.names = FALSE
    )
    gsave_ids <- gsub("\\.rds$", "", gsave_ids) # strip extension
  
    manifest <- .gdsrc_json_read(p, consolidate = TRUE)$content
    art_ids <- names(manifest)
  
    # find directly protected artifacts
    # (tagged with existing giottosave)
    is_direct_protect <- vapply(manifest, function(art) {
        any(art$giottosave %in% gsave_ids)
    }, FUN.VALUE = logical(1L))
    protected <- art_ids[is_direct_protect]
  
    # recursively protect dependencies
    # BFS: walk `depends` tags to find all transitively protected
    # artifacts
    queue <- protected
    visited <- character(0L)
    while (length(queue) > 0L) {
        current <- queue[[1L]]
        queue <- queue[-1L]
        if (current %in% visited) next
        visited <- c(visited, current)
      
        deps <- manifest[[current]]$depends
        if (!is.null(deps) && length(deps) > 0L) {
            # only follow deps that exist in the manifest
            new_deps <- intersect(deps, art_ids)
            new_deps <- setdiff(new_deps, visited)
            queue <- c(queue, new_deps)
        }
    }
    protected <- visited
  
    # prune unprotected artifacts
    to_prune <- setdiff(art_ids, protected)
    
    for (uid in to_prune) {
        art_dir <- .gdsrc_artifact_dir(p, uid)
        if (dir.exists(art_dir)) {
            unlink(art_dir, recursive = TRUE, force = TRUE)
        }
    }
  
    # update manifest
    manifest[to_prune] <- NULL
    .gdsrc_json_write(p, manifest)
  
    message("[GiottoDisk] Pruned ", length(to_prune), 
        " unprotected artifacts")
    invisible(to_prune)
}


# .gdsrc_dir_prune <- function(gsrc) {
#     gsrc <- gsrc@read()
#     artifacts <- .gdsrc_json_artifacts(gsrc)
#     unversioned <- artifacts[is.na(version)]
#     for (f in unversioned$fullpaths) {
#         if (file.exists(f)) {
#             file.remove(f)
#         } else {
#             vmsg("File not found (already deleted?):", f)
#         }
#     }
#     vmsg("Removed ", nrow(unversioned), " unversioned artifacts")

#     # update json
#     keep_ids <- artifacts[!is.na(version), id]
#     gsrc@catalog$artifacts <- gsrc@catalog$artifacts[keep_ids]
#     gsrc@write()

#     TRUE
# }
