#' @name artifact_dump
#' @title GiottoDisk artifact dump management
#' @description
#' Utilities to control the location of automatically created
#' GiottoDisk artifacts. Some GiottoDisk processes perform operation
#' that materialize to disk. The dump is used as the default location
#' artifacts will be written to. When setting objects into the
#' Giotto object, these objects may be reclaimed from the dump and
#' placed under control of the `gsource` managing the object.
#'
#' Passing a `gDirSource` object points the dump directly at the
#' source's vault (`artifacts/`). Artifacts are then written to their
#' final vault location immediately, making adoption a pure registration
#' step with no file movement.
#' @param x `character` path to directory, or a `gDirSource` object.
#'   When a `gDirSource` is supplied, the dump is set to the vault
#'   directory (`<x@path>/artifacts/`). When missing, resets to the
#'   default session temp location.
#' @param verbose verbosity
#' @returns artifact dump directory path (invisibly)
NULL

#' @rdname artifact_dump
#' @export
setArtifactDumpDir <- function(x, verbose = NULL) {
    dir <- if (missing(x)) {
        file.path(tempdir(), "gdisk_dump")
    } else if (inherits(x, "gDirSource")) {
        .gdsrc_vault_dir(x@path)
    } else {
        x
    }
    vmsg(.v = verbose, "[GiottoDisk] Setting artifact dump:\n", dir)
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }
    options("giottodisk.artifact_dump" = dir)
    invisible(dir)
}

#' @rdname artifact_dump
#' @export
getArtifactDumpDir <- function() {
    getOption("giottodisk.artifact_dump", file.path(tempdir(), "gdisk_dump"))
}

.dump_tempfile <- function() {
    dir <- getArtifactDumpDir()
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    tempfile(tmpdir = dir)
}

# get number of rows in data
# as.numeric is enforced for uniform output
.dplyr_nrow <- function(data) {
    nr <- nrow(data)
    if (!is.na(nr) && !is.null(nr)) return(as.numeric(nr))
    # fallback
    data |>
        dplyr::summarize(n = n()) |>
        dplyr::pull(n, as_vector = TRUE) |>
        as.numeric()
}

# determine if a function has no content
.is_empty_fun <- function(f) {
    body_code <- deparse(body(f))
    return(identical(body_code, "NULL") ||
               identical(body_code, "{\n    NULL\n}"))
}

#' @name .dt_set_row_index
#' @title data.table Utility: Add row index
#' @description
#' Internal that adds a sequential row index column to a data.frame using
#' data.table's by-reference modification. This modifies the input object
#' in-place for memory efficiency with large in-memory datasets.
#'
#' This is the in-memory counterpart to [.arrow_add_row_index()] for use with
#' materialized data.frames rather than lazy arrow queries.
#' @param x `data.frame` or `data.table`. Will be converted to data.table if
#'   not already.
#' @param offset `integer`. Value to add to the row index. Default is `1`,
#'   meaning the first row will be assigned index `1`. Use higher offsets when
#'   processing data chunks that continue from previous batches.
#' @param col `character`. Name of the column to create for the row index.
#'   Default is `"row_index"`.
#' @return The modified data.table with the new index column added. Note that
#'   the input `x` is also modified by reference.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Add row_index starting from 1
#' df <- data.frame(x = 1:100, y = letters[1:100])
#' .dt_set_row_index(df)
#'
#' # Add row_index starting from 1000 (for tile/chunk continuation)
#' .dt_set_row_index(df, offset = 1000)
#'
#' # Custom column name
#' .dt_set_row_index(df, col = "id")
#' }
.dt_set_row_index <- function(x, offset = 1L, col = "row_index") {
    checkmate::assert_data_frame(x)
    x <- data.table::setDT(x)
    x[, (col) := as.integer(.I + offset)]
    x
}

#' @name artifact_uid
#' @aliases uid
#' @title Artifact Unique Identifier
#' @description
#' Artifacts are tracked in [gDirSource-class] using an unique identifier.
#' These are randomized alphanumeric strings of length 8 that can
#' additionally be tagged with the process ID (pid) and node ID
#' that generated it. Only the pid is tagged by default, but the
#' node can also be tagged to help prevent collisions on multi-
#' node setups by setting the appropriate option below:
#' 
#' * `giottodisk.uid_include_pid` (default = TRUE)
#' * `giottodisk.uid_include_node` (default = FALSE)
#' 
#' These uids are not affected by user seed setting and instead use
#' a temporary random seed based on the time and a process-specific
#' counter.
NULL

.make_uid <- function(n = 8L, include_node = NULL, include_pid = NULL) {
    include_node <- include_node %||% 
        getOption("giottodisk.uid_include_node", FALSE)
    include_pid <- include_pid %||% 
        getOption("giottodisk.uid_include_pid", TRUE)
    count <- getOption("giottodisk.uid_count", 1L)
    on.exit({
        options("giottodisk.uid_count" = count + 1L)
    })
    GiottoUtils::gwith_seed(seed = Sys.time() + count, {
        sampleset <- sample(c(letters, LETTERS, as.character(0:9)), 
            size = n,
            replace = TRUE
        )
        rand <- paste(sampleset, collapse = "")
    })

    parts <- rand
    if (include_pid) {
        pid <- Sys.getpid()
        parts <- c(pid, parts)
    }
    if (include_node) {
        host <- substr(Sys.info()[["nodename"]], 1, 8)
        parts <- c(host, parts)
    }
    paste0(parts, collapse = "_")
}

# generate geom_id from available special cols (source_id, tile_index, row_index).
# called at write time before special cols are dropped, so all components are present.
.auto_geom_id <- function(data, prefix) {
    paste(prefix, data$row_index, sep = "_")
}

# generate a character timestamp
.timestamp <- function() {
    as.character(format(Sys.time()))
}

# generate hash
.hash <- function(x) {
    digest::digest(x)
}

# Move a path (file or directory) to a new location.
# Prefers file.rename (atomic, instant, same-fs). Falls back to
# recursive copy + delete when rename fails (e.g. cross-device moves).
# Use this over bare file.rename whenever the source and destination
# may be on different filesystems.
.move_path <- function(from, to) {
    if (isTRUE(file.rename(from, to))) return(invisible(to))
    # cross-device fallback
    if (dir.exists(from)) {
        # file.copy(dir, parent, recursive=TRUE) copies dir INTO parent as basename(from)
        to_parent <- dirname(to)
        from_name <- basename(from)
        to_name <- basename(to)
        file.copy(from, to_parent, recursive = TRUE)
        if (from_name != to_name) {
            file.rename(file.path(to_parent, from_name), to)
        }
    } else {
        file.copy(from, to, overwrite = FALSE)
    }
    unlink(from, recursive = TRUE, force = TRUE)
    invisible(to)
}
