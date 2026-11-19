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
#' * `giottoLargeImage`: as `SpatRaster`, for `@file_path` (for a JPEG-2000
#'   image, the `.ome.tif` rather than the VRT GiottoClass reads it through).
#'   An image made from a `SpatRaster` has no `@file_path`; its raster is used.
#' * Union stores: `TRUE` only if all substores are contained.
#' @param src `gsource` object
#' @param store object to test (`fileStore`, `SpatRaster`, `giottoLargeImage`,
#'   or union store)
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
    .gdsrc_in_vault(src@path, terra::sources(store))
})

#' @rdname sourceContains
#' @export
setMethod("sourceContains", signature("gDirSource", "giottoLargeImage"),
    function(src, store, ...) {
    .gdsrc_in_vault(src@path, .img_source_file(store))
})

# FALSE for an in-memory raster (no source file)
.gdsrc_in_vault <- function(p, f) {
    f <- normalizePath(f, mustWork = FALSE)
    if (all(!nzchar(f))) return(FALSE)
    vdir <- normalizePath(.gdsrc_vault_dir(p), mustWork = FALSE)
    all(startsWith(f, paste0(vdir, "/")))
}

# The file an image's pixels come from: @file_path. Only an image made from a
# SpatRaster has none, and then it is the raster's own file ("" in memory).
# terra::sources() is not the image when @file_path is set: GiottoClass reads
# a JPEG-2000 TIFF through a VRT it writes under tempdir() (GiottoClass
# adr/0005) and rebuilds on reconnect().
.img_source_file <- function(img) {
    if (.img_has_file_path(img)) {
        return(normalizePath(img@file_path, mustWork = FALSE))
    }
    terra::sources(img@raster_object)
}

.img_has_file_path <- function(img) {
    fp <- img@file_path
    is.character(fp) && length(fp) == 1L && !is.na(fp) && nzchar(fp)
}

#' @rdname sourceContains
#' @export
setMethod("sourceContains", signature("gDirSource", "unionParquetStore"),
    function(src, store, ...) {
    all(vapply(store@stores, function(s) sourceContains(src, s),
        FUN.VALUE = logical(1L)))
})

#' @rdname sourceContains
#' @details For a `giotto` object: whether all the stores it holds are in the
#'   vault (expression matrices, points, polygons, overlaps and networks), so
#'   that [snapshotSave()] would move none of them into it. In-memory data is
#'   written into the snapshot's object file and images are copied, so neither
#'   counts. When `FALSE`, the attribute `external` names the stores outside.
#' @export
setMethod("sourceContains", signature("gDirSource", "giotto"),
    function(src, store, ...) {
    external <- .gdsrc_external_stores(src, store)
    if (!length(external)) return(TRUE)
    structure(FALSE, external = external)
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
#'   copied into the vault. A fresh uid is assigned. Returns the updated
#'   `SpatRaster` pointing to the vault path. A raster read through a VRT
#'   errors: the VRT points at image files rather than holding them.
#' * `giottoLargeImage`: copies `@file_path`, so a JPEG-2000 image is copied
#'   as the original `.ome.tif`, not the VRT it is read through. An image made
#'   from a `SpatRaster` (no `@file_path`) adopts the raster. Returns the image
#'   with `@file_path` naming the vault copy, reopened through
#'   [GiottoClass::reconnect()]. Errors when `@file_path` is a VRT outside the
#'   vault.
#' * Union stores: each substore is adopted independently.
#' @param src `gsource` object
#' @param store object to adopt (`fileStore`, `SpatRaster`, `giottoLargeImage`,
#'   or union store)
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
    hash <- .hash(storeRead(.store_nostate(store)))
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
# parquetEdgeStore: parent fileStore method moves <root>/ — that
# physically relocates BOTH the edges/ and nodes/ subdirs since
# they live under the moved parent. The in-memory @nodes handle
# (a parquetStore at <old_root>/nodes/) is stale after the move; we
# re-point it at the new path to keep storeRead working. Re-pointing rather
# than rebuilding keeps its @uid, so the network's fingerprint survives.
setMethod("sourceAdopt", signature("gDirSource", "parquetEdgeStore"),
    function(src, store, meta = NULL, giottosave = NULL, depends = NULL, ...) {
    store <- callNextMethod() # fileStore: moves files, updates @path
    nodes_path <- file.path(store@path, "nodes")
    if (is.null(store@nodes)) {
        store@nodes <- parquetStore(path = nodes_path)
    } else {
        store@nodes@path <- nodes_path
    }
    store
})

#' @rdname sourceAdopt
#' @export
setMethod("sourceAdopt", signature("gDirSource", "SpatRaster"),
    function(src, store, meta = NULL, giottosave = NULL, ...) {
    p <- src@path
    r <- store

    f <- normalizePath(terra::sources(r), mustWork = FALSE)
    is_in_memory <- all(!nzchar(f))

    # already in vault -- nothing to do
    if (.gdsrc_in_vault(p, f)) return(r)
    # copying a VRT copies a pointer to the image, not the image
    if (any(grepl("\\.vrt$", f, ignore.case = TRUE))) {
        stop("[sourceAdopt] this SpatRaster is read through a VRT (", f[[1L]],
            "), which only points at image files and cannot be copied into ",
            "the vault. For an image GiottoClass reads through a VRT, adopt ",
            "the giottoLargeImage: its @file_path names the image file.",
            call. = FALSE)
    }

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
        r <- .create_terra_spatraster(savepath)
    }

    .gdsrc_json_add_artifact(p,
        store_type = "IMAGE",
        uid = uid,
        hash = fingerprint_data(r),
        meta = meta,
        giottosave = giottosave %||% NA_character_
    )
    r
})

#' @rdname sourceAdopt
#' @export
setMethod("sourceAdopt", signature("gDirSource", "giottoLargeImage"),
    function(src, store, meta = NULL, giottosave = NULL, ...) {
    img <- store
    r <- img@raster_object
    f <- .img_source_file(img)
    # reconnect() reapplies the tracked extent to the reopened file
    img@extent <- terra::ext(r)[]

    if (!.img_has_file_path(img)) {
        # made from a SpatRaster: the raster is the image
        r <- sourceAdopt(src, r, meta = meta, giottosave = giottosave)
        img@file_path <- terra::sources(r)
        return(GiottoClass::reconnect(img))
    }

    p <- src@path
    # a VRT in the vault is already broken: the project needs the image redone
    if (.gdsrc_in_vault(p, f)) return(img)
    if (grepl("\\.vrt$", f, ignore.case = TRUE)) {
        stop("[sourceAdopt] image '", GiottoClass::objName(img),
            "' has a VRT as its @file_path (", f, "). A VRT only points at ",
            "image files, so it cannot be copied into the vault: create the ",
            "image from the image file.", call. = FALSE)
    }
    savepath <- .gdsrc_allocate_artifact_dir(p, create = TRUE)
    uid <- names(savepath)
    names(savepath) <- NULL
    file.copy(f, savepath)
    img@file_path <- savepath
    img <- GiottoClass::reconnect(img)

    .gdsrc_json_add_artifact(p,
        store_type = "IMAGE",
        uid = uid,
        hash = fingerprint_data(img@raster_object),
        meta = meta,
        giottosave = giottosave %||% NA_character_
    )
    img
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

# The stores of giotto object `x` that snapshotSave() would adopt into
# `src`, named by kind and object. Mirrors the .ss_gdsrc_register_external_*
# steps: keep the two in step.
.gdsrc_external_stores <- function(src, x) {
    out <- character()
    outside <- function(store) !sourceContains(src, store)
    name <- function(kind, obj) sprintf("%s '%s'", kind, GiottoClass::objName(obj))
    for (e in x[["expression"]]) {
        if (inherits(e[], "IterableMatrix") && outside(e[]))
            out <- c(out, name("expression", e))
    }
    for (p in x[["feat_info"]]) {
        stores <- .feat_info_stores(p)
        if (length(stores) && any(vapply(stores, outside, logical(1L))))
            out <- c(out, name("points", p))
    }
    for (p in x[["spatial_info"]]) {
        if (inherits(p[], "dataStore") && outside(p[]))
            out <- c(out, name("polygons", p))
        for (f in setdiff(names(p@overlaps), "intensity")) {
            o <- p@overlaps[[f]]
            if (inherits(o, "overlapPointDisk") && outside(o@data))
                out <- c(out, sprintf("overlaps '%s/%s'", GiottoClass::objName(p), f))
        }
    }
    for (slot in c("nn_network", "spatial_network")) {
        for (n in x[[slot]]) {
            if (inherits(n@network, "dataStore") && outside(n@network))
                out <- c(out, name("network", n))
        }
    }
    out
}
