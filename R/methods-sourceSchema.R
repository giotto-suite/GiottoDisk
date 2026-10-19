#' @include AllGenerics.R
NULL

# Schema versions of gDirSource project directories: see design.Rmd
# "Schema versions". A layout change bumps .GDSRC_SCHEMA_VERSION and adds its
# upgrade to .GDSRC_SCHEMA_UPGRADES.

.GDSRC_SCHEMA_VERSION <- 2L
.GDSRC_SCHEMA_MIN <- 2L

# per-session cache of on-disk versions, keyed by project path. The version
# only changes through sourceUpdate(), which refreshes it.
.gdsrc_schema_env <- new.env(parent = emptyenv())

.gdsrc_schema_key <- function(p) normalizePath(p, mustWork = FALSE)

# on-disk schema version, NA for a directory without a giottodir.json
.gdsrc_schema_ondisk <- function(p, cache = TRUE) {
    key <- .gdsrc_schema_key(p)
    if (cache && !is.null(v <- .gdsrc_schema_env[[key]])) return(v)
    f <- .gdsrc_json_path(p)
    if (!file.exists(f)) return(NA_integer_)
    v <- as.integer(.json_read(f)$header$versions$manifest %||% NA_integer_)
    .gdsrc_schema_env[[key]] <- v
    v
}

.gdsrc_schema_forget <- function(p) {
    key <- .gdsrc_schema_key(p)
    if (exists(key, envir = .gdsrc_schema_env, inherits = FALSE)) {
        rm(list = key, envir = .gdsrc_schema_env)
    }
    invisible(NULL)
}

# Refuse a project whose layout this GiottoDisk does not read or write.
.gdsrc_check_schema <- function(p, site = "[GiottoDisk]") {
    v <- .gdsrc_schema_ondisk(p)
    if (is.na(v)) return(invisible(TRUE)) # new project, stamped on first write
    if (v >= .GDSRC_SCHEMA_MIN && v <= .GDSRC_SCHEMA_VERSION) {
        return(invisible(TRUE))
    }

    msg <- if (v < .GDSRC_SCHEMA_MIN) {
        sprintf(paste0(
            "%s project directory schema v%d is older than this GiottoDisk ",
            "supports (v%d).\n  Run sourceUpdate(\"%s\") to upgrade it."
        ), site, v, .GDSRC_SCHEMA_MIN, p)
    } else {
        sprintf(paste0(
            "%s project directory schema v%d was written by a newer ",
            "GiottoDisk (this one writes v%d).\n  Update GiottoDisk."
        ), site, v, .GDSRC_SCHEMA_VERSION)
    }
    if (isTRUE(getOption("giottodisk.schema_check", TRUE))) {
        stop(msg, "\n  options(giottodisk.schema_check = FALSE) proceeds ",
            "anyway, at your own risk.", call. = FALSE)
    }
    warning(msg, "\n  Proceeding: giottodisk.schema_check is FALSE.",
        call. = FALSE)
    invisible(FALSE)
}

# sourceSchema ####

#' @name sourceSchema
#' @title Schema Version of a Project Source
#' @description
#' Report the layout version of a [gsource] backend and whether this
#' version of GiottoDisk can work with it. A project older than the minimum
#' supported version is upgraded with [sourceUpdate()].
#' @param src `gsource` object or `character` filepath to a project
#'   directory.
#' @param ... additional params to pass (none implemented)
#' @returns named list:
#'   * `version` - schema version on disk (`NA` if the directory holds no
#'     project yet)
#'   * `current` - schema version this GiottoDisk writes
#'   * `min` - oldest schema version this GiottoDisk works with
#'   * `supported` - whether `version` is within `min` and `current`
#' @examples
#' src <- gDirSource(file.path(tempdir(), "schema_example"))
#' sourceSchema(src)$version
#' @export
setMethod("sourceSchema", signature("character"), function(src, ...) {
    if (!dir.exists(src)) {
        stop("[sourceSchema] not an existing project directory\n",
            call. = FALSE)
    }
    sourceSchema(sourceCreate(src, type = "gDirSource"), ...)
})

#' @rdname sourceSchema
#' @export
setMethod("sourceSchema", signature("gDirSource"), function(src, ...) {
    v <- .gdsrc_schema_ondisk(src@path, cache = FALSE)
    list(
        version = v,
        current = .GDSRC_SCHEMA_VERSION,
        min = .GDSRC_SCHEMA_MIN,
        supported = !is.na(v) &&
            v >= .GDSRC_SCHEMA_MIN && v <= .GDSRC_SCHEMA_VERSION
    )
})

# sourceUpdate ####

#' @name sourceUpdate
#' @title Upgrade a Project Source to the Current Schema
#' @description
#' Migrate a [gsource] backend written by an older GiottoDisk to the layout
#' this version reads and writes, then stamp the new schema version. A
#' project already at the current version is left untouched.
#'
#' Migrations run in order, one version at a time, and each rewrites the
#' manifest only after the files it moves are in place, so an interrupted
#' upgrade can be rerun.
#' @inheritParams sourceSchema
#' @param verbose verbosity
#' @returns the `gsource`, invisibly
#' @export
setMethod("sourceUpdate", signature("character"), function(src, ...) {
    if (!dir.exists(src)) {
        stop("[sourceUpdate] not an existing project directory\n",
            call. = FALSE)
    }
    sourceUpdate(sourceCreate(src, type = "gDirSource"), ...)
})

#' @rdname sourceUpdate
#' @export
setMethod("sourceUpdate", signature("gDirSource"), function(src,
    verbose = NULL, ...) {
    p <- src@path
    v <- .gdsrc_schema_ondisk(p, cache = FALSE)
    if (is.na(v) || v == .GDSRC_SCHEMA_VERSION) return(invisible(src))
    if (v > .GDSRC_SCHEMA_VERSION) {
        stop(sprintf(paste0(
            "[sourceUpdate] project directory schema v%d was written by a ",
            "newer GiottoDisk (this one writes v%d)\n"
        ), v, .GDSRC_SCHEMA_VERSION), call. = FALSE)
    }

    while (v < .GDSRC_SCHEMA_VERSION) {
        upgrade <- .GDSRC_SCHEMA_UPGRADES[[as.character(v)]]
        if (is.null(upgrade)) {
            stop(sprintf("[sourceUpdate] no upgrade from schema v%d\n", v),
                call. = FALSE)
        }
        vmsg(.v = verbose, sprintf(
            "[GiottoDisk] upgrading project schema v%d -> v%d", v, v + 1L
        ))
        upgrade(p)
        v <- v + 1L
    }
    invisible(src)
})

## upgrades ####

# v1 -> v2: flat snapshot files become folders with a fixed object file name,
# and each artifact's `<name>` tag becomes the folder's `<name>.<id>`.
# `.qs` snapshots (qs is no longer read) are moved too, so their artifacts
# stay protected from sourcePrune.
.gdsrc_upgrade_1_to_2 <- function(p) {
    # pending edits carry v1 tags; fold them in before the rename
    .gdsrc_json_consolidate(p, check = FALSE)

    gsdir <- .gdsrc_giottosave_dir(p)
    if (dir.exists(gsdir)) {
        flat <- list.files(gsdir,
            pattern = "\\.(rds|qs|qs2)$", ignore.case = TRUE,
            full.names = TRUE
        )
        flat <- flat[!dir.exists(flat)]
        for (f in flat) {
            ext <- tolower(file_extension(f)[[1L]])
            name <- sub("\\.(rds|qs|qs2)$", "", basename(f),
                ignore.case = TRUE
            )
            dest <- file.path(gsdir, .gdsrc_snapshot_tag(name))
            dir.create(dest)
            if (!file.rename(f, file.path(dest, .gdsrc_snapshot_objfile(ext)))) {
                unlink(dest, recursive = TRUE)
                stop("[sourceUpdate] could not move snapshot ", f, "\n",
                    call. = FALSE)
            }
        }
    }

    # Map tags from every snapshot folder, not only the ones moved above, so a
    # rerun after an interruption still finishes the retag.
    snaps <- .gdsrc_snapshots(gsdir)
    retag <- stats::setNames(snaps$tag, snaps$name)
    manifest <- .gdsrc_json_read(p)$content
    for (uid in names(manifest)) {
        tags <- unlist(manifest[[uid]]$giottosave)
        hit <- !is.na(tags) & tags %in% names(retag)
        if (!any(hit)) next
        tags[hit] <- retag[tags[hit]]
        manifest[[uid]]$giottosave <- tags
    }
    .gdsrc_json_write(p, manifest, version = 2L)
    invisible(TRUE)
}

.GDSRC_SCHEMA_UPGRADES <- list(
    "1" = .gdsrc_upgrade_1_to_2
)
