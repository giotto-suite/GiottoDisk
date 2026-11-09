
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
#' Contained artifacts are tracked via a manifest (the registry) with a
#' record per artifact.
#' 
#' Writers queue edits as files in `_pending/`, one file per edit, and a
#' read folds them into the registry as its next version. A version is
#' published by creating its file only if no other writer created it first,
#' so processes sharing a project never overwrite each other's commits, and
#' no lock is needed. A record's fields are written once, when it is
#' created; its snapshot tags (`giottosave`) and dependencies (`depends`)
#' are only ever added to, and dead tags are dropped by [sourcePrune()].
#' 
#' # Artifact tracking
#' 
#' Data can be written to the backing directory using `sourceWrite()`.
#' This serializes the data to disk using the preferred store types
#' for the `gsource` class. The artifact is also assigned an
#' [artifact_uid] at this point.
#' 
#' # Directory Structure
#' 
#' * `giottodir.<n>.json` - the registry, one file per version; the highest
#' `n` is current and the one below it is kept.
#' * `giottodir.json` - stub carrying the project's schema version.
#' * `artifacts` - vault directory containing the actual data artifacts
#' which are within subdirectories named by their `uid`.
#' * `_pending` - directory containing pending edits to the manifest.
#' * `giottosave` - one folder per snapshot, `giottosave/<name>.<id>/`,
#' holding the serialized giotto object (`gobject.RDS` or `gobject.qs2`)
#' that references assets controlled by the `gDirSource`.
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
sourceCreate <- function(path, type = "gDirSource", ...) {
    type <- match.arg(tolower(type), choices = c("gdirsource"))
    switch(type,
        "gdirsource" = gDirSource(path, ...)
    )
}

# initialize ####

setMethod("initialize", signature("gDirSource"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    p <- .Object@path
    if (is.na(p)) {
        stop("[gDirSource] 'path' should be a directory path for the giotto project.\n",
            call. = FALSE
        )
    }
    if (!file.exists(.gdsrc_json_path(p))) {
        .gdsrc_registry_init(p)
    }
    # normalize only once the directory exists: normalizePath() returns a
    # missing path unchanged, which would leave a new project relative to
    # whatever getwd() is when its artifacts are later read (knitr, workers)
    .Object@path <- normalizePath(p, mustWork = TRUE)
    .Object
})

# methods ####

setMethod("show", signature("gDirSource"), function(object) {
    cat(sprintf("<%s>\n", class(object)))

    gsave_ids <- .gdsrc_snapshot_names(
        .gdsrc_giottosave_dir(object@path)
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

# A record's fields are fixed when it is created; only snapshot tags can be
# added afterwards (adr/0023)
.gdsrc_fixed_fields <- function() {
    stop("[gDirSource] registry fields are fixed when an artifact is ",
        "recorded; only snapshot tags can be added, with ",
        "src[uid, \"giottosave\"] <- tag\n", call. = FALSE)
}

setMethod("[<-", signature("gDirSource", i = "character", j = "missing", value = "ANY"), function(x, i, j, ..., value) {
    .gdsrc_fixed_fields()
})

setMethod("[<-", signature("gDirSource", i = "character", j = "character", value = "ANY"), function(x, i, j, ..., value) {
    checkmate::assert_atomic(value)
    if (!identical(j, "giottosave")) .gdsrc_fixed_fields()
    .gdsrc_json_add_tags(x@path, uid = i, tag = value)
    x
})

setMethod("[<-", signature("gDirSource", i = "missing", j = "character", value = "ANY"), function(x, i, j, ..., value) {
    checkmate::assert_atomic(value)
    if (!identical(j, "giottosave")) .gdsrc_fixed_fields()
    uids <- names(.gdsrc_json_read(x@path, consolidate = TRUE)$content)
    .gdsrc_json_add_tags(x@path, uid = uids, tag = value)
    x
})

# source check ####

#' @name resolveSource
#' @title Source Detection and Regeneration
#' @description
#' Detect if an input path is to a known [gsource] type.
#' If so, regenerate the source object. Otherwise, return
#' FALSE.
#' @param path `character` filepath to managed directory
#' @returns `gsource` inheriting object if recognized source type,
#' `FALSE` if not.
#' @export
resolveSource <- function(path) {
    if (!dir.exists(path)) return(FALSE)
    toplevel <- list.files(path, recursive = FALSE, full.names = FALSE)
    if ("giottodir.json" %in% toplevel) {
        return(sourceCreate(path, type = "gDirSource"))
    }
    FALSE
}

# json internal tools ####

#' @name giotto_json
#' @title Giotto Directory JSON
#' @description
#' Internal utilities for reading and editing the [gDirSource] registry
#' (`giottodir.<n>.json`).
#' @keywords internal
NULL

## path setting ####

# p is giotto directory. `giottodir.json` is a stub carrying only the schema
# version, so a GiottoDisk that predates schema v3 refuses the project; the
# registry itself is giottodir.<n>.json (adr/0023).
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
#' * schema version of the project directory (see [sourceSchema()])
#' 
#' A registry version also carries its number (`seq`), a hash of its
#' content, and the pending edit files it applied (`applied`).
.gdsrc_json_header <- function(version = .GDSRC_SCHEMA_VERSION, seq = NULL,
    hash = NULL, applied = NULL) {
    header <- list(
        "description" = "[GiottoDisk] Artifacts Manifest",
        "versions" = list(
            "package" = as.character(packageVersion("GiottoDisk")),
            "manifest" = as.integer(version)
        )
    )
    if (!is.null(seq)) {
        header$seq <- as.integer(seq)
        header$hash <- hash
        header$applied <- as.character(applied)
    }
    header
}

# test seam: tests interleave another writer at a named point
.gdsrc_hook <- function(point) {
    h <- getOption("giottodisk.test_hook")
    if (is.function(h)) h(point)
    invisible()
}

## registry versions ####

# The registry is a series of whole files, giottodir.<n>.json. A commit
# writes version n + 1 by hard-linking a complete temp file to its name,
# which fails if that name exists: of writers that read version n, one
# commits and the others re-read and retry. A reader takes the highest number
# whose header agrees with its name. Two versions are kept (adr/0023).

.GDSRC_REGISTRY_RE <- "^giottodir\\.([0-9]{8,})\\.json$"
.GDSRC_COMMIT_ATTEMPTS <- 50L

.gdsrc_registry_file <- function(p, n) {
    file.path(p, sprintf("giottodir.%08d.json", as.integer(n)))
}

# version numbers on disk, ascending
.gdsrc_registry_versions <- function(p) {
    f <- list.files(p, pattern = .GDSRC_REGISTRY_RE)
    sort(as.integer(sub(.GDSRC_REGISTRY_RE, "\\1", f)))
}

# The newest readable version as list(n, top, header, content). A version
# that does not parse, or whose header disagrees with its name (a write cut
# short by a crash), gives way to the one below. Only a read error counts:
# a warning raised mid-read (R closing an unused connection, say) must not,
# since basing a commit on the version below would drop the newest's edits.
# `top` is the highest number on disk, readable or not, which the next commit
# goes above. A version listed and then removed by another writer's cleanup
# is read again from a fresh listing.
.gdsrc_registry_load <- function(p) {
    for (attempt in 1:3) {
        ns <- .gdsrc_registry_versions(p)
        for (n in rev(ns)) {
            reg <- tryCatch(.json_read(.gdsrc_registry_file(p, n)),
                error = function(e) NULL
            )
            if (identical(as.integer(reg$header$seq), n)) {
                return(list(n = n, top = max(ns), header = reg$header,
                    content = reg$content %||% list()))
            }
        }
    }
    stop("[gDirSource] no readable registry (giottodir.<n>.json) in:\n  ", p,
        call. = FALSE)
}

# edit files the versions on disk record as applied
.gdsrc_registry_applied <- function(p) {
    unlist(lapply(.gdsrc_registry_versions(p), function(n) {
        reg <- tryCatch(.json_read(.gdsrc_registry_file(p, n)),
            error = function(e) NULL
        )
        unlist(reg$header$applied)
    }))
}

.gdsrc_content_hash <- function(content) {
    digest::digest(jsonlite::toJSON(content), algo = "xxhash64",
        serialize = FALSE
    )
}

.gdsrc_readonly <- function(p) {
    structure(class = c("gdsrc_readonly", "error", "condition"), list(
        message = paste0("[gDirSource] project directory is not writable:\n  ",
            p),
        call = NULL
    ))
}

# Write `content` as version n + 1, recording the edit files it applied; `n`
# is the highest version on disk when the content was read. TRUE when
# committed, FALSE when another writer took n + 1 first. Signals a
# `gdsrc_readonly` condition when the directory cannot be written.
.gdsrc_registry_commit <- function(p, n, content, applied = character()) {
    if (file.access(p, mode = 2L) != 0L) stop(.gdsrc_readonly(p))
    target <- .gdsrc_registry_file(p, n + 1L)
    if (file.exists(target)) return(FALSE)
    data <- list(
        header = .gdsrc_json_header(seq = n + 1L,
            hash = .gdsrc_content_hash(content), applied = applied
        ),
        content = content
    )
    temp <- tempfile("giottodir_", tmpdir = p, fileext = ".tmp")
    on.exit(unlink(temp, force = TRUE), add = TRUE)
    jsonlite::write_json(data, temp)
    .gdsrc_hook("commit")
    if (!.gdsrc_publish(temp, target)) return(FALSE)
    .gdsrc_registry_cleanup(p, n + 1L)
    TRUE
}

# Hard-link `from` to `to`; FALSE if `to` exists. On a filesystem without
# hard links (FAT, some network shares and synced folders), option
# giottodisk.commit_by_rename publishes by rename, which is safe only while a
# single process writes the project.
.gdsrc_publish <- function(from, to) {
    if (isTRUE(getOption("giottodisk.commit_by_rename", FALSE))) {
        if (file.exists(to)) return(FALSE)
        return(file.rename(from, to))
    }
    ok <- .gdsrc_link(from, to)
    if (ok || file.exists(to)) return(ok)
    stop("[gDirSource] could not commit the registry in:\n  ", dirname(to),
        "\nThe filesystem may not support hard links, which commits rely on.",
        "\noptions(giottodisk.commit_by_rename = TRUE) commits by rename",
        " instead, which is safe only while one process writes the project.",
        call. = FALSE
    )
}

# file.link(), apart so tests can stand in a filesystem without hard links
.gdsrc_link <- function(from, to) suppressWarnings(file.link(from, to))

# After version m is in place: the edit files applied by the versions below
# it are in m too, so they go, then those versions, keeping m - 1 for a
# reader that listed it a moment ago. Edit files go before their version, so
# an interruption leaves them still listed as applied.
.gdsrc_registry_cleanup <- function(p, m) {
    pending <- .gdsrc_json_pending_dir(p)
    ns <- .gdsrc_registry_versions(p)
    for (n in ns[ns < m]) {
        reg <- tryCatch(.json_read(.gdsrc_registry_file(p, n)),
            error = function(e) NULL
        )
        applied <- unlist(reg$header$applied)
        if (length(applied)) {
            unlink(file.path(pending, applied), force = TRUE)
        }
        if (n < m - 1L) unlink(.gdsrc_registry_file(p, n), force = TRUE)
    }
    invisible()
}

# A new project: version 1, empty, then the stub. An interrupted setup is
# finished by the next one.
.gdsrc_registry_init <- function(p) {
    if (!dir.exists(p)) {
        message("Setting up Giotto project directory at: \n", p)
        dir.create(path = p)
    }
    if (!length(.gdsrc_registry_versions(p))) {
        .gdsrc_registry_commit(p, 0L, structure(list(), names = character()))
    }
    .gdsrc_write_stub(p)
}

# giottodir.json as a stub: the schema version and nothing else
.gdsrc_write_stub <- function(p) {
    .json_atomic_write(list(
        header = .gdsrc_json_header(),
        content = structure(list(), names = character()),
        note = "The registry is giottodir.<n>.json, the highest n."
    ), .gdsrc_json_path(p))
}

## reading ####

#' @describeIn giotto_json Read the registry as list(header, content), with
#' pending edits applied when `consolidate = TRUE`. A project in an older
#' layout is read only with the schema check turned off.
.gdsrc_json_read <- function(p, consolidate = FALSE) {
    if (!dir.exists(p)) {
        stop("[gDirSource] Path to Giotto project directory does not exist\n",
            call. = FALSE
        )
    }
    if (!length(.gdsrc_registry_versions(p))) {
        .gdsrc_check_schema(p)
        return(.gdsrc_legacy_read(p, consolidate = consolidate))
    }
    if (consolidate) return(.gdsrc_json_consolidate(p))
    .gdsrc_registry_load(p)
}

## edits ####

# Fields holding sets; edits only add to them. Every other field is written
# once, when the record is created.
.GDSRC_SET_FIELDS <- c("giottosave", "depends")
.GDSRC_EDIT_RE <- "^edit_.+\\.json$"

#' @describeIn giotto_json Queue an edit as its own file in `_pending/`,
#' written whole and renamed in, so no two writers share one. `create`: the
#' fields of a new record; on an existing record only the fields it lacks are
#' filled. `add`: values joined into the set fields (`giottosave`,
#' `depends`), for one or more `uid`s. The result depends neither on the
#' order edits are applied in nor on applying one twice (adr/0023).
.gdsrc_json_edit <- function(p, uid, create = NULL, add = NULL) {
    checkmate::assert_character(uid, min.len = 1L, any.missing = FALSE)
    checkmate::assert_list(create, null.ok = TRUE)
    checkmate::assert_list(add, null.ok = TRUE)
    if (!is.null(create) && length(uid) != 1L) {
        stop("[gDirSource] a record is created for one uid at a time\n",
            call. = FALSE)
    }
    bad <- setdiff(names(add), .GDSRC_SET_FIELDS)
    if (length(bad)) {
        stop("[gDirSource] only ", paste(.GDSRC_SET_FIELDS, collapse = ", "),
            " can be added to; not: ", paste(bad, collapse = ", "), "\n",
            call. = FALSE)
    }
    .gdsrc_check_schema(p)
    if (!length(.gdsrc_registry_versions(p))) {
        stop("[gDirSource] project directory is in an older layout; ",
            "run sourceUpdate(\"", p, "\") before writing to it\n",
            call. = FALSE)
    }
    d <- .gdsrc_json_pending_dir(p)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    host <- substr(gsub("[^A-Za-z0-9-]", "-", Sys.info()[["nodename"]]), 1, 16)
    file <- file.path(d, sprintf("edit_%s_%s_%s.json", host, Sys.getpid(),
        .make_uid(n = 10L, include_pid = FALSE)
    ))
    entry <- list(header = .gdsrc_json_header(), uid = uid, create = create,
        add = add)
    .json_atomic_write(entry[!vapply(entry, is.null, NA)], file)
    invisible(TRUE)
}

#' @describeIn giotto_json Add a tracked artifact to the directory. This adds
#' the tracking metadata for the object and places it in a location specified
#' by `.gdsrc_json_pending_dir()` to await consolidation into the registry.
#' @param store_type character. Type of file format storage (e.g. \"h5\")
#' @param uid character. Unique ID for artifact tracking
#' @param hash character. Hash of the in-memory store object to track changes
#' @param meta list (optional). Additional list of atomic object(s) that can
#' @param giottosave `character` tagged giotto save(s) if any.
#' be attached as further metadata to the particular uid
#' @keywords internal
.gdsrc_json_add_artifact <- function(p, store_type, uid, hash, meta = NULL, giottosave = NA_character_, depends = NULL) {
    checkmate::assert_character(store_type)
    checkmate::assert_character(uid)
    checkmate::assert_character(hash)
    checkmate::assert_list(meta, null.ok = TRUE)
    checkmate::assert_character(giottosave)
    checkmate::assert_character(depends, null.ok = TRUE)

    create <- c(list(
        "time" = .timestamp(),
        "store" = store_type,
        "hash" = hash
    ), meta)
    # giottosave is always present on a record, NA while untagged
    add <- list("giottosave" = giottosave)
    if (length(depends) > 0L) add[["depends"]] <- depends
    .gdsrc_json_edit(p = p, uid = uid, create = create, add = add)
}

#' @describeIn giotto_json Tag artifacts with a snapshot tag
.gdsrc_json_add_tags <- function(p, uid, tag) {
    if (!length(uid)) return(invisible(FALSE))
    .gdsrc_json_edit(p = p, uid = unique(uid), add = list(giottosave = tag))
}

# one edit applied to the registry content
.gdsrc_apply_edit <- function(content, edit) {
    for (uid in unlist(edit$uid)) {
        rec <- content[[uid]] %||% list()
        for (f in names(edit$create)) {
            if (is.null(rec[[f]])) rec[[f]] <- edit$create[[f]]
        }
        for (f in names(edit$add)) {
            v <- unique(as.character(c(unlist(rec[[f]]), unlist(edit$add[[f]]))))
            v <- v[!is.na(v)]
            rec[[f]] <- if (length(v)) v else if (f == "giottosave") NA
        }
        content[[uid]] <- rec
    }
    content
}

# pending edit files not yet in the registry
.gdsrc_pending_files <- function(p, applied = .gdsrc_registry_applied(p)) {
    d <- .gdsrc_json_pending_dir(p)
    if (!dir.exists(d)) return(character())
    f <- list.files(d, pattern = "\\.json$")
    legacy <- f[!grepl(.GDSRC_EDIT_RE, f)]
    if (length(legacy)) {
        stop("[gDirSource] pending edits written by an older GiottoDisk in:\n  ",
            d, "\nRun sourceUpdate(\"", p, "\") with no other process using ",
            "the project.\n", call. = FALSE)
    }
    setdiff(f, applied)
}

#' @describeIn giotto_json Apply pending edits and commit the result as the
#' next registry version, retrying when another writer commits first.
#' Returns the registry as `.gdsrc_json_read()` does. A project that cannot
#' be written is returned with the edits applied, uncommitted.
.gdsrc_json_consolidate <- function(p, check = TRUE) {
    d <- .gdsrc_json_pending_dir(p)
    for (attempt in seq_len(.GDSRC_COMMIT_ATTEMPTS)) {
        reg <- .gdsrc_registry_load(p)
        files <- .gdsrc_pending_files(p)
        edits <- lapply(file.path(d, files), function(f) {
            tryCatch(.json_read(f), error = function(e) NULL)
        })
        if (!length(files)) {
            # every pending edit is in a version above the one read: read again
            if (reg$top < max(.gdsrc_registry_versions(p))) next
            return(reg)
        }
        # a file gone between listing and reading was applied and removed by
        # a commit after `reg`: read the registry again
        ok <- !vapply(edits, is.null, NA)
        if (!any(ok)) next
        content <- Reduce(.gdsrc_apply_edit, edits[ok], reg$content)
        if (check) .gdsrc_check_schema(p)
        committed <- tryCatch(
            .gdsrc_registry_commit(p, reg$top, content, applied = files[ok]),
            gdsrc_readonly = function(e) NA
        )
        if (is.na(committed)) {
            reg$content <- content
            return(reg)
        }
        if (committed) return(.gdsrc_registry_load(p))
        .gdsrc_backoff(attempt)
    }
    stop("[gDirSource] could not commit the registry in ",
        .GDSRC_COMMIT_ATTEMPTS, " attempts: other processes kept committing ",
        "first\n", call. = FALSE)
}

.gdsrc_backoff <- function(attempt) {
    Sys.sleep(stats::runif(1L, 0, min(0.5, 0.005 * 2^min(attempt, 7L))))
}

## older layouts ####

# Before schema v3 the registry was giottodir.json itself, with pending edits
# as `<uid>.json` merged field by field. Kept for sourceUpdate() and for
# reading with the schema check off.

.gdsrc_legacy_read <- function(p, consolidate = FALSE) {
    f <- .gdsrc_json_path(p)
    if (!file.exists(f)) {
        stop("[gDirSource] giottodir.json not found in:\n", p, "\n",
            call. = FALSE)
    }
    data <- .json_read(f)
    if (consolidate) data$content <- .gdsrc_legacy_apply(p, data$content)
    data
}

.gdsrc_legacy_apply <- function(p, content) {
    d <- .gdsrc_json_pending_dir(p)
    edits <- list.files(d, pattern = "\\.json$", full.names = TRUE)
    edits <- edits[!grepl(.GDSRC_EDIT_RE, basename(edits))]
    for (f in edits) {
        uid <- sub("\\.json$", "", basename(f))
        content[[uid]] <- utils::modifyList(content[[uid]] %||% list(),
            .json_read(f)$content)
    }
    content
}

# fold `<uid>.json` edits into giottodir.json; run by sourceUpdate() only
.gdsrc_legacy_consolidate <- function(p) {
    d <- .gdsrc_json_pending_dir(p)
    edits <- list.files(d, pattern = "\\.json$", full.names = TRUE)
    edits <- edits[!grepl(.GDSRC_EDIT_RE, basename(edits))]
    if (!length(edits)) return(invisible(TRUE))
    data <- .json_read(.gdsrc_json_path(p))
    .gdsrc_legacy_write(p, .gdsrc_legacy_apply(p, data$content),
        version = data$header$versions$manifest %||% .GDSRC_SCHEMA_MIN)
    unlink(edits, force = TRUE)
    invisible(TRUE)
}

.gdsrc_legacy_write <- function(p, x, version) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
    .json_atomic_write(
        list(header = .gdsrc_json_header(version), content = x),
        .gdsrc_json_path(p)
    )
}

## snapshot layout ####

# giottosave/<name>.<id>/gobject.<ext>; the folder name is the tag. See
# design.Rmd "Snapshots and tags".

.GDSRC_SNAPSHOT_DIR <- "^(.+)\\.([a-z0-9]{8})$"
.GDSRC_SNAPSHOT_OBJ <- "^gobject\\.(rds|qs|qs2)$"

.gdsrc_snapshot_tag <- function(name) {
    id <- tolower(.make_uid(n = 8L, include_node = FALSE, include_pid = FALSE))
    paste0(name, ".", id)
}

.gdsrc_snapshot_objfile <- function(method) {
    switch(tolower(method), "rds" = "gobject.RDS", paste0("gobject.", method))
}

# every snapshot in a giottosave dir: tag, name, object file, and mtime. A
# folder counts only once its object file is in place, so a half-written save
# protects nothing; dot-prefixed folders are staging.
.gdsrc_snapshots <- function(gsdir) {
    none <- data.frame(tag = character(), name = character(),
        file = character(), mtime = as.POSIXct(character())
    )
    if (!dir.exists(gsdir)) return(none)
    tags <- list.dirs(gsdir, full.names = FALSE, recursive = FALSE)
    tags <- tags[!startsWith(tags, ".") & grepl(.GDSRC_SNAPSHOT_DIR, tags)]
    files <- vapply(tags, function(tg) {
        f <- list.files(file.path(gsdir, tg), pattern = .GDSRC_SNAPSHOT_OBJ,
            ignore.case = TRUE, full.names = TRUE
        )
        if (length(f) == 0L) NA_character_ else f[[1L]]
    }, character(1L), USE.NAMES = FALSE)
    keep <- !is.na(files)
    if (!any(keep)) return(none)
    data.frame(
        tag = tags[keep],
        name = sub(.GDSRC_SNAPSHOT_DIR, "\\1", tags[keep]),
        file = files[keep],
        mtime = file.info(files[keep])$mtime
    )
}

# tags that protect artifacts from sourcePrune
.gdsrc_snapshot_tags <- function(gsdir) .gdsrc_snapshots(gsdir)$tag

# user-facing snapshot names
.gdsrc_snapshot_names <- function(gsdir) unique(.gdsrc_snapshots(gsdir)$name)

# The snapshot to load for `name`, or the most recently written one if
# `name` is NULL; NULL if there is none. Recency is the object file's mtime,
# since a folder's own moves whenever anything inside it changes. More than
# one folder per name exists only between an overwrite's write and its
# cleanup, and the newer is the overwrite.
.gdsrc_find_snapshot <- function(gsdir, name = NULL) {
    snaps <- .gdsrc_snapshots(gsdir)
    if (!is.null(name)) snaps <- snaps[snaps$name == name, , drop = FALSE]
    if (nrow(snaps) == 0L) return(NULL)
    as.list(snaps[which.max(snaps$mtime), , drop = FALSE])
}

# A name becomes part of a folder name, so it must stay one path component
.gdsrc_assert_snapshot_name <- function(name, site = "[snapshotSave]") {
    checkmate::assert_string(name, min.chars = 1L)
    if (grepl("[/\\\\]", name) || startsWith(name, ".")) {
        stop(sprintf(
            "%s invalid snapshot name '%s': %s", site, name,
            "it may not contain a path separator or start with '.'"
        ), call. = FALSE)
    }
    invisible(name)
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
#' @exportS3Method base::as.data.frame
as.data.frame.gDirSource <- function(x, ...) {
    data <- as.list.gDirSource(x)
    # one row per giottosave tag; depends stays a single list cell
    data <- lapply(data, function(art) {
        if (!is.null(art$depends)) art$depends <- list(art$depends)
        art
    })
    data.table::rbindlist(data,
        fill = TRUE, idcol = "uid"
    )
}

## pruning ####

# Artifacts reachable from a live snapshot: tagged by one, or depended on,
# transitively, by one that is
.gdsrc_protected <- function(manifest, live_tags) {
    art_ids <- names(manifest)
    tagged <- vapply(manifest, function(art) {
        any(unlist(art$giottosave) %in% live_tags)
    }, FUN.VALUE = logical(1L))
    queue <- art_ids[tagged]
    visited <- character(0L)
    while (length(queue) > 0L) {
        current <- queue[[1L]]
        queue <- queue[-1L]
        if (current %in% visited) next
        visited <- c(visited, current)
        deps <- intersect(unlist(manifest[[current]]$depends), art_ids)
        queue <- c(queue, setdiff(deps, visited))
    }
    visited
}

# uids named by pending edits not yet in the registry
.gdsrc_pending_uids <- function(p) {
    d <- .gdsrc_json_pending_dir(p)
    unique(unlist(lapply(file.path(d, .gdsrc_pending_files(p)), function(f) {
        tryCatch(unlist(.json_read(f)$uid), error = function(e) NULL)
    })))
}

# Drop artifacts no live snapshot protects, and tags whose snapshot is gone.
# The registry change is committed first and the vault folders removed after,
# so an interruption leaves stray folders, never records of missing data. An
# artifact with an edit queued after the registry was read is left alone.
# Not safe against a snapshot being saved while it runs (adr/0023).
.gdsrc_artifact_prune <- function(p) {
    .gdsrc_check_schema(p, site = "[sourcePrune]")
    gsave_dir <- .gdsrc_giottosave_dir(p)
    for (attempt in seq_len(.GDSRC_COMMIT_ATTEMPTS)) {
        reg <- .gdsrc_json_consolidate(p)
        manifest <- reg$content
        live <- .gdsrc_snapshot_tags(gsave_dir)
        protected <- c(.gdsrc_protected(manifest, live), .gdsrc_pending_uids(p))

        vault_ids <- list.files(.gdsrc_vault_dir(p), recursive = FALSE,
            full.names = FALSE)
        to_prune <- setdiff(vault_ids, protected)

        kept <- manifest[setdiff(names(manifest), to_prune)]
        kept <- lapply(kept, function(art) {
            tags <- unlist(art$giottosave)
            tags <- tags[!is.na(tags) & tags %in% live]
            art$giottosave <- if (length(tags)) tags else NA
            art
        })
        if (identical(kept, manifest)) {
            committed <- TRUE
        } else {
            .gdsrc_hook("prune")
            committed <- tryCatch(.gdsrc_registry_commit(p, reg$top, kept),
                gdsrc_readonly = function(e) stop(conditionMessage(e),
                    call. = FALSE))
        }
        if (committed) {
            for (uid in to_prune) {
                unlink(.gdsrc_artifact_dir(p, uid), recursive = TRUE,
                    force = TRUE)
            }
            message("[GiottoDisk] Pruned ", length(to_prune),
                " unprotected artifacts")
            return(invisible(to_prune))
        }
        .gdsrc_backoff(attempt)
    }
    stop("[sourcePrune] could not commit the registry in ",
        .GDSRC_COMMIT_ATTEMPTS, " attempts: other processes kept committing ",
        "first\n", call. = FALSE)
}
