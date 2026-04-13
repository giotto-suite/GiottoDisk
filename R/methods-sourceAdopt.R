#' @include methods-sourceWrite.R
NULL

#' @name sourceContains
#' @title Test if a Store is Managed by a Source
#' @description
#' Returns `TRUE` if the store is already registered in the artifact vault
#' of a [gsource].
#'
#' * `fileStore`: checks whether the store's uid is present in the manifest.
#' * `SpatRaster`: checks whether all source file paths are inside the vault
#'   directory. In-memory rasters always return `FALSE`.
#' * Union stores: `TRUE` only if all substores are contained.
#' @param src `gsource` object
#' @param store object to test (`fileStore`, `SpatRaster`, or union store)
#' @param ... additional params (unused)
#' @returns `logical(1)`
#' @export
setMethod("sourceContains", signature("gDirSource", "fileStore"),
    function(src, store, ...) {
    !is.null(src[store@uid])
})

#' @rdname sourceContains
#' @export
setMethod("sourceContains", signature("gDirSource", "SpatRaster"),
    function(src, store, ...) {
    f <- normalizePath(terra::sources(store), mustWork = FALSE)
    if (all(!nzchar(f))) return(FALSE) # in-memory
    vdir <- normalizePath(.gdsrc_vault_dir(src@path))
    all(startsWith(f, paste0(vdir, "/")))
})

#' @rdname sourceContains
#' @export
setMethod("sourceContains", signature("gDirSource", "unionParquetStore"),
    function(src, store, ...) {
    all(vapply(store@stores, function(s) sourceContains(src, s),
        FUN.VALUE = logical(1L)))
})

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

