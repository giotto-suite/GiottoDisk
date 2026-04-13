# docs ####

#' @name sourceWrite
#' @title Write Data to a Source
#' @description
#' Write data to a [gsource] that manages data formats and file
#' organization. Docs for specific `gsource` implementations have further
#' information.
#'
#' * [sourceWrite-gDirSource]
#'
#' @param src `gsource` object. Used to provide default save locations
#' and save formats for a managed Giotto backend.
#' @param data data to write
#' @param meta `list`. Additional metadata to attach to this Giotto backend
#' @param store_type `character`. Store type to write as
#' @param ... additional params passed to `storeWrite` method
#' managed artifact.
#' @returns written `store` object
#' @family sourceWrite methods
NULL

#' @name sourceWrite-gDirSource
#' @title Write Data to a `gDirSource`
#' @inheritParams sourceWrite
#' @description
#' Write data to a [gDirSource] managed project directory. Defaults for
#' different data types are settable via global options:
#' 
#' * *sparse matrices:* `giotto.gdsrc_sparsematrix_format` (default = "bpcells")
#' * *dense matrices:* `"giotto.gdsrc_densematrix_format"` (default = "h5")
#' * *dataframes:* `"giotto.gdsrc_dataframe_format"` (default = "parquet")
#' * *spatvector:* `"giotto.gdsrc_spatvector_format"` (default = "parquetGeom")
#'
#' @examples
#' testdir <- file.path(tempdir(), "testdirsource")
#' gsrc <- sourceCreate(path = testdir)
#' @returns written `store` object
#' @family sourceWrite methods
NULL

# definitions ####

#' @rdname sourceWrite-gDirSource
#' @examples
#' # arrays ----------------------------------------- #
#' m <- matrix(1:9, nrow = 3)
#' m_written <- sourceWrite(gsrc, m)
#' storeRead(m_written)
#' @export
setMethod("sourceWrite", signature("gDirSource", "memoryMatrix"),
    function(src, data, meta = NULL, store_type = NULL, ...) {
        if (is.null(store_type)) {
            store_type <- if (.is_sparse_matrix(data)) {
                getOption("giotto.gdsrc_sparsematrix_format", "bpcells")
            } else {
                getOption("giotto.gdsrc_densematrix_format", "h5")
            }
        }

        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })

.is_sparse_matrix <- function(mat) {
    if (inherits(mat, "sparseMatrix")) return(TRUE)
    if (inherits(mat, "matrix")) return(sum(mat != 0) / length(mat) < 0.1)
    FALSE
}

#' @rdname sourceWrite-gDirSource
#' @examples
#' # spatvector ----------------------------------------- #
#' sv <- terra::vect(system.file("ex/lux.shp", package="terra"))
#' sv_written <- sourceWrite(gsrc, sv)
#' storeRead(sv_written) # default is output = "query"
#' storeRead(sv_written, output = "sf")
#' storeRead(sv_written, output = "terra")
#' storeRead(sv_written, output = "tibble")
#' @export
setMethod("sourceWrite", signature("gDirSource", "SpatVector"),
    function(src, data, meta = NULL, store_type = getOption("giotto.gdsrc_spatvector_format", "parquetGeom"), ...) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })

#' @rdname sourceWrite-gDirSource
#' @examples
#' # tables ----------------------------------------- #
#' store <- sourceWrite(gsrc, iris)
#' force(store)
#' storeRead(store) # default output is an arrow query
#' storeRead(store, output = "tibble") # pull into memory as tibble
#' @export
setMethod("sourceWrite", signature("gDirSource", "data.frame"),
    function(src, data, meta = NULL, store_type = getOption("giotto.gdsrc_dataframe_format", "parquet"), ...) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })

setMethod("sourceWrite", signature("gDirSource", "fileStore"),
    function(src, data, meta = NULL, store_type, ...) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })


# sourceAdopt ####

#' @name sourceAdopt
#' @title Adopt a Store into a Source
#' @description
#' Move an existing store into the managed artifact vault of a [gsource],
#' registering it in the manifest.
#'
#' * `fileStore`: the store's `@uid` is preserved so existing `source_id`
#'   partition references remain valid. Store must already be written to disk.
#' * `SpatRaster`: in-memory rasters are written as COG; on-disk rasters are
#'   moved into the vault. A fresh uid is assigned. Returns the updated
#'   `SpatRaster` pointing to the vault path.
#' * Union stores: each substore is adopted independently.
#' @param src `gsource` object
#' @param store object to adopt (`fileStore`, `SpatRaster`, or union store)
#' @param meta `list` (optional). Additional metadata to attach to the artifact entry.
#' @param giottosave `character` (optional). Giottosave name to tag the artifact with.
#' @param ... additional params (unused)
#' @returns updated store/raster object pointing to the managed vault location
#' @export
setMethod("sourceAdopt", signature("gDirSource", "fileStore"),
    function(src, store, meta = NULL, giottosave = NULL, ...) {
    if (!storeExists(store)) {
        stop("[sourceAdopt] store does not exist at path:\n  ", store@path,
            call. = FALSE)
    }
    old_path <- store@path
    uid <- store@uid
    new_path <- .gdsrc_allocate_artifact_dir(src@path, uid = uid, create = TRUE)
    .move_path(old_path, new_path)
    store@path <- new_path
    hash <- .hash(storeRead(store))
    .gdsrc_json_add_artifact(src@path,
        store_type = class(store),
        uid = uid,
        hash = hash,
        meta = meta,
        giottosave = giottosave %||% NA_character_
    )
    store
})

#' @rdname sourceAdopt
#' @export
setMethod("sourceAdopt", signature("gDirSource", "SpatRaster"),
    function(src, store, meta = NULL, giottosave = NULL, ...) {
    p <- src@path
    r <- store
    vdir <- normalizePath(.gdsrc_vault_dir(p))

    f <- normalizePath(terra::sources(r), mustWork = FALSE)
    is_in_memory <- all(!nzchar(f))

    # already in vault — nothing to do
    if (!is_in_memory && all(startsWith(f, paste0(vdir, "/")))) return(r)

    savepath <- .gdsrc_allocate_artifact_dir(p, create = TRUE)
    uid <- names(savepath)
    names(savepath) <- NULL

    if (is_in_memory) {
        r <- terra::writeRaster(r,
            filename = savepath,
            filetype = "COG",
            NAflag = NA,
            overwrite = FALSE
        )
    } else {
        if (length(f) != 1L) {
            stop("[sourceAdopt] multi-source SpatRasters are not supported",
                call. = FALSE)
        }
        .move_path(f, savepath)
        r <- terra::rast(savepath)
    }

    capture.output(show(r)) # ping before hashing — terra hash changes on first access
    hash <- .hash(r)
    .gdsrc_json_add_artifact(p,
        store_type = "IMAGE",
        uid = uid,
        hash = hash,
        meta = meta,
        giottosave = giottosave %||% NA_character_
    )
    r
})

#' @rdname sourceAdopt
#' @export
setMethod("sourceAdopt", signature("gDirSource", "unionParquetStore"),
    function(src, store, meta = NULL, giottosave = NULL, ...) {
    store@stores <- lapply(store@stores, function(s) {
        sourceAdopt(src, s, meta = meta, giottosave = giottosave, ...)
    })
    store
})

# internals ####

# return allocated artifact save path.
# value is named with the uid
.gdsrc_allocate_artifact_dir <- function(p, uid = NULL, basename = "data", create = TRUE) {
    checkmate::assert_character(p)
    checkmate::assert_character(uid, null.ok = TRUE)
    # setup save info
    if (is.null(uid)) {
        uid <- .make_uid()
    }
    savedir <- .gdsrc_artifact_dir(p, uid = uid)
    if (create && !dir.exists(savedir)) dir.create(savedir, recursive = TRUE)
    savepath <- normalizePath(file.path(savedir, basename), mustWork = FALSE)
    stats::setNames(savepath, uid)
}

# Move a path to a new location. Tries file.rename first (atomic, same-fs).
# Falls back to copy + delete for cross-device moves.
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

.gdsrc_write_artifact <- function(p, data, store_type, meta = NULL, ...) {
    dstore <- storeCreate(type = store_type)
    uid <- dstore@uid
    savepath <- .gdsrc_allocate_artifact_dir(p, uid = uid, create = TRUE)
    dstore@path <- savepath # update save path
    dstore <- storeWrite(dstore, data, ...)
    # record hash of delayed representation
    hash <- .hash(storeRead(dstore))
    # record artifact entry
    .gdsrc_json_add_artifact(p,
        store_type = store_type,
        uid = uid,
        hash = hash,
        meta = meta
    )
    dstore
}
