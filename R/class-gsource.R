
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
#' @slot catalog list. Structured list recording what assets are available
#' within the project directory.
#' @slot read function. Read from the `giottodir.json` to populate the
#' `catalog` slot. Returns a new `gDirSource` object.
#' @slot write function. Write the `catalog` slot contents to `giottodir.json`
#' @description
#' Extends [gsource]. Used to designate and manage a directory as a
#' giotto project directory containing on-disk data artifacts. Artifacts are
#' organized by subdirectories into by specific filetype. The toplevel of a
#' managed directory also contains a `giottodir.json` which is a catalog of the
#' contained artifacts.
#'
#' Data can be written to the backing directory using the `gsource` as a store
#' with `storeWrite()`.
#'
#' Note that methods of editing `giottodir.json` contents are currently
#' internal only. See [giotto_json]
#' @export
setClass("gDirSource",
    contains = "gsource",
    slots = list(
        catalog = "list",
        read = "function",
        write = "function"
    )
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
    etag <- "[gDirSource]" # tag source of error
    json_path <- file.path(p, "giottodir.json")
    if (is.na(p)) {
        stop(etag,
            " 'path' should be a directory path for the giotto project.\n",
            call. = FALSE
        )
    }

    # default catalog subdirectories
    .Object@catalog <- .gdsrc_catalog_defaults(.Object@catalog)

    write_fun <- function() { # write a giotto dir json
        if (!dir.exists(p)) {
            message("Setting up Giotto project directory at:", p)
            dir.create(path = p)
        }
        jsonlite::write_json(.Object@catalog,
            path = json_path
        )
        invisible(TRUE)
    }
    .Object@write <- write_fun

    read_fun <- function() { # read a giotto dir json as gDirSource
        if (!dir.exists(p)) {
            stop(etag,
                " Path to Giotto project directory does not exist\n",
                call. = FALSE
            )
        }
        if (!file.exists(json_path)) {
            stop(sprintf("%s giottodir.json not found in:\n%s\n", etag, p),
                 call. = FALSE)
        }
        .Object@catalog <- jsonlite::fromJSON(json_path)
        .Object
    }
    .Object@read <- read_fun

    if (!dir.exists(p) || !file.exists(json_path)) {
        .Object@write()
    }

    .Object
})


# default catalog ####

# update the catalog list of a gDirSource object with the default subdirectory
# names for specific formats if no alternatives have already been provided.
.gdsrc_catalog_defaults <- function(clist = list()) {
    checkmate::assert_list(clist)
    clist$stores$h5 <- clist$stores$h5 %null% "h5"
    clist$stores$spatvector <- clist$stores$spatvector %null% "spatvector"
    clist$stores$parquet <- clist$stores$parquet %null% "parquet"
    clist
}


# json ####

#' @name giotto_json
#' @title Giotto Directory JSON
#' @description
#' Internal utilities for editing the [gDirSource] managed
#' `giottodir.json`.
#' These functions accept a `giotto` object and return a `giotto` object.
#' @param gsrc `gDirSource`
#' @keywords internal
NULL

#' @describeIn giotto_json Hub function for editing the json. Business logic for
#' ensuring that edits follow these steps:
#'
#' 1. reading what exists on-disk into the `gDirSource` catalog
#' 2. performing the edit on the catalog using `fun`
#' 3. writing the catalog as a new version of the json
#' @param fun function. Function to edit catalog contents before writing.
#' @keywords internal
.gdsrc_json_edit_content <- function(gsrc, fun) {
    gsrc <- gsrc@read()
    gsrc@catalog <- fun(gsrc@catalog)
    gsrc <- initialize(gsrc) # needed to update the write function
    gsrc@write()
    gsrc
}

#' @describeIn giotto_json Add a tracked asset to the directory.
#' @param store_type character. Type of file format storage (e.g. \"h5\")
#' @param uid character. Unique ID for asset tracking
#' @param hash character. Hash of the in-memory store object to track changes
#' @param meta list (optional). Additional list of atomic object(s) that can
#' be attached as further metadata to the particular uid
#' @keywords internal
.gdsrc_json_add_artifact <- function(gsrc, store_type, uid, hash, meta = NULL) {
    checkmate::assert_character(store_type)
    checkmate::assert_character(uid)
    checkmate::assert_character(hash)
    checkmate::assert_list(meta, null.ok = TRUE)

    entry <- list(
        time = .timestamp(),
        store = store_type,
        version = NA_character_,
        hash = hash
    )

    .gdsrc_json_edit_content(gsrc, function(x) {
        x$artifacts[[uid]] <- entry
        if (!is.null(meta)) {
            x$meta[[uid]] <- meta
        }
        x
    })
}

#' @describeIn giotto_json Add a file store_type to the directory.
#' @keywords internal
.gdsrc_json_add_store <- function(gsrc, store_type, path) {
    checkmate::assert_character(store_type)
    checkmate::assert_character(path)
    gsrc <- gsrc@read()
    # return early if store already present
    if (store_type %in% names(gsrc@catalog$stores)) return(gsrc)

    basepath <- gsrc@path
    fullpath <- file.path(basepath, path)
    if (!dir.exists(fullpath)) {
        vmsg("creating", fullpath, "...")
        dir.create(fullpath)
    }

    entry <- list(
        path = path
    )

    .gdsrc_json_edit_content(gsrc, function(x) {
        x$stores[[store_type]] <- entry
        x
    })
}

#' @describeIn giotto_json Establish a giotto project saved version ID.
#' @keywords internal
.gdsrc_json_add_project_version <- function(gsrc,
    uid = paste0("giottosave_", .make_uid())) {
    checkmate::assert_character(uid)

    entry <- list(
        time = .timestamp()
    )

    .gdsrc_json_edit_content(gsrc, function(x) {
        x$versions[[uid]] <- entry
        x
    })
}

#' @describeIn giotto_json Tag specific artifacts (based on uid) as belonging to
#' a particular saved version of the giotto project
#' @keywords internal
.gdsrc_json_artifact_tag_version <- function(gsrc, artifacts, version) {
    checkmate::assert_character(artifacts)
    checkmate::assert_character(version)

    .gdsrc_json_edit_content(gsrc, function(x) {
        for (id in artifacts) {
            x$artifacts[[id]]$version <- version
        }
        x
    })
}

#' @describeIn giotto_json Get table of artifacts information
.gdsrc_json_artifacts <- function(gsrc) {
    gsrc <- gsrc@read()
    version <- store <- NULL # NSE var
    artifacts <- gsrc@catalog$artifacts
    artifacts <- data.table::rbindlist(artifacts, idcol = "id")
    stores <- gsrc@catalog$stores
    stores <- data.frame(
        store = names(stores),
        path = unlist(stores),
        stringsAsFactors = FALSE
    )
    artifacts <- merge(artifacts, stores, by = "store")
    artifacts[, "fullpaths" := file.path(gsrc@path, path, id)]
    artifacts
}

# pruning ####

.gdsrc_dir_prune <- function(gsrc) {
    gsrc <- gsrc@read()
    artifacts <- .gdsrc_json_artifacts(gsrc)
    unversioned <- artifacts[is.na(version)]
    for (f in unversioned$fullpaths) {
        if (file.exists(f)) {
            file.remove(f)
        } else {
            vmsg("File not found (already deleted?):", f)
        }
    }
    vmsg("Removed ", nrow(unversioned), " unversioned artifacts")

    # update json
    keep_ids <- artifacts[!is.na(version), id]
    gsrc@catalog$artifacts <- gsrc@catalog$artifacts[keep_ids]
    gsrc@write()

    TRUE
}
