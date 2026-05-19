#' @include methods-sourceWrite.R
NULL

# session-scoped cache: maps old_path -> new_path for adopted leaf dirs
# reset at the start of each snapshotSave to avoid stale entries across saves
.adopt_session_map <- new.env(parent = emptyenv())

.adopt_session_lookup <- function(old_path) .adopt_session_map[[old_path]]
.adopt_session_record <- function(old_path, new_path) {
    .adopt_session_map[[old_path]] <- new_path
}
.adopt_session_reset <- function() {
    rm(list = ls(.adopt_session_map), envir = .adopt_session_map)
}

# IterableMatrix leaf traversal primitives ####

# returns all leaf @dir paths from any IterableMatrix structure
.im_leaf_dirs <- function(x) {
    if ("dir" %in% slotNames(x)) {
        slot(x, "dir")
    } else if ("matrix_list" %in% slotNames(x)) {
        unlist(lapply(slot(x, "matrix_list"), .im_leaf_dirs))
    } else if ("matrix" %in% slotNames(x)) {
        .im_leaf_dirs(slot(x, "matrix"))
    } else {
        character(0L)
    }
}

# applies f to each leaf IterableMatrix (has @dir), rebuilds the structure
.im_map_leaves <- function(x, f) {
    if ("dir" %in% slotNames(x)) {
        f(x)
    } else if ("matrix_list" %in% slotNames(x)) {
        slot(x, "matrix_list") <- lapply(slot(x, "matrix_list"),
            function(m) .im_map_leaves(m, f))
        x
    } else if ("matrix" %in% slotNames(x)) {
        slot(x, "matrix") <- .im_map_leaves(slot(x, "matrix"), f)
        x
    } else {
        warning("[sourceAdopt] unknown compound IterableMatrix type: ",
            class(x), call. = FALSE)
        x
    }
}

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
    function(src, store, meta = NULL, giottosave = NULL, depends = NULL, ...) {
    if (!storeExists(store)) {
        stop("[sourceAdopt] store does not exist at path:\n  ", store@path,
            call. = FALSE)
    }
    old_path <- store@path
    uid <- store@uid
    new_path <- .gdsrc_allocate_artifact_dir(src@path, uid = uid, create = TRUE)
    .move_path(old_path, new_path)
    store@path <- new_path
    hash <- .hash(storeRead(storeBase(store)))
    .gdsrc_json_add_artifact(src@path,
        store_type = class(store),
        uid = uid,
        hash = hash,
        meta = meta,
        giottosave = giottosave %||% NA_character_,
        depends = depends
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

    # already in vault -- nothing to do
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
        file.copy(f, savepath)
        r <- terra::rast(savepath)
    }

    capture.output(show(r)) # ping before hashing -- terra hash changes on first access
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

# * ANY fallbacks ####

#' @rdname sourceContains
#' @export
setMethod("sourceContains", signature("gDirSource", "ANY"),
    function(src, store, ...) {
    if (inherits(store, "IterableMatrix")) {
        dirs <- normalizePath(.im_leaf_dirs(store), mustWork = FALSE)
        if (length(dirs) == 0L) return(FALSE)
        vdir <- normalizePath(.gdsrc_vault_dir(src@path))
        all(startsWith(dirs, paste0(vdir, "/")))
    } else {
        FALSE
    }
})

#' @rdname sourceAdopt
#' @export
setMethod("sourceAdopt", signature("gDirSource", "ANY"),
    function(src, store, meta = NULL, giottosave = NULL, ...) {
    if (inherits(store, "IterableMatrix")) {
        .im_map_leaves(store, function(leaf) {
            old_path <- normalizePath(slot(leaf, "dir"), mustWork = FALSE)
            if (.any_adopt_in_vault(src, old_path)) {
                uid <- basename(dirname(old_path))
                if (is.null(src[uid])) {
                    .gdsrc_json_add_artifact(src@path,
                        store_type = class(leaf),
                        uid = uid,
                        hash = .hash(BPCells::open_matrix_dir(old_path)),
                        meta = meta,
                        giottosave = giottosave %||% NA_character_
                    )
                }
                return(leaf)
            }
            if (!is.null(cached <- .adopt_session_lookup(old_path))) {
                slot(leaf, "dir") <- cached
                return(leaf)
            }
            if (!.any_adopt_external_ok(old_path)) return(leaf)
            new_path <- .gdsrc_allocate_artifact_dir(src@path, create = TRUE)
            uid <- names(new_path)
            names(new_path) <- NULL
            .move_path(old_path, new_path)
            .adopt_session_record(old_path, new_path)
            slot(leaf, "dir") <- new_path
            .gdsrc_json_add_artifact(src@path,
                store_type = class(leaf),
                uid = uid,
                hash = .hash(BPCells::open_matrix_dir(new_path)),
                meta = meta,
                giottosave = giottosave %||% NA_character_
            )
            leaf
        })
    } else {
        stop("[sourceAdopt] no method for class: ", class(store), call. = FALSE)
    }
})

.any_adopt_in_vault <- function(src, path) {
    vdir <- normalizePath(.gdsrc_vault_dir(src@path))
    startsWith(path, paste0(vdir, "/"))
}

.any_adopt_external_ok <- function(path) {
    dump_dir <- normalizePath(getOption("giottodisk.artifact_dump"), mustWork = FALSE)
    if (startsWith(path, paste0(dump_dir, "/"))) return(TRUE)
    if (isTRUE(getOption("giottodisk.adopt_external", FALSE))) return(TRUE)
    warning(
        "[sourceAdopt] skipping external store at:\n  ", path, "\n",
        "Set options(giottodisk.adopt_external = TRUE) to adopt.",
        call. = FALSE
    )
    FALSE
}
