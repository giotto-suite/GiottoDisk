# Image files that tile readers (GiottoLens, Viv, QuPath) cannot use as they
# are -- large, and stored in strips or without overviews, like MERSCOPE's
# mosaics -- converted once to a COG (tiled, with overviews) in the project
# vault. The COG is an IMAGE artifact like any adopted image: a snapshot finds
# it in the vault and tags it, and it is reused while its source file is
# unchanged.

# TRUE when a TIFF can be read by tile at any zoom: small enough to read
# whole, or tiled with overviews. A file GDAL cannot describe is left as is.
.tif_tile_readable <- function(path, max_whole = 4096) {
    d <- tryCatch(terra::describe(path), error = function(e) character())
    num <- function(pattern) {
        m <- regmatches(d, regexec(pattern, d))
        as.numeric(unlist(Filter(length, m))[2:3])
    }
    size <- num("^Size is ([0-9]+), ([0-9]+)")
    if (anyNA(size)) return(TRUE)
    if (max(size) <= max_whole) return(TRUE)
    block <- num("Block=([0-9]+)x([0-9]+)")
    # A strip spans the full width
    tiled <- !anyNA(block) && block[[1L]] < size[[1L]]
    tiled && any(grepl("^\\s*Overviews:", d))
}

# Path of a COG of `path` in the vault of `gsource`, written on first use:
# deflate with a predictor, 512-pixel tiles, the source's pixel type and
# values and row order. Keyed by the source's path, size and modification
# time, so a conversion is reused only while it matches its source.
.image_cog <- function(path, gsource, verbose = NULL) {
    p <- gsource@path
    info <- file.info(path)
    src <- normalizePath(path, winslash = "/", mustWork = TRUE)
    key <- .hash(list(
        path = src,
        size = as.numeric(info$size),
        mtime = as.numeric(info$mtime),
        what = "image-cog-v1"
    ))

    manifest <- .gdsrc_json_read(p, consolidate = TRUE)$content
    for (uid in names(manifest)) {
        a <- manifest[[uid]]
        if (!identical(a$store, "IMAGE") || !identical(a$source_key, key)) next
        f <- file.path(.gdsrc_artifact_dir(p, uid), "data")
        if (file.exists(f)) {
            GiottoUtils::vmsg(.v = verbose, "[image_cog] using converted image:", f)
            return(normalizePath(f))
        }
    }

    GiottoUtils::vmsg(.v = verbose,
        "[image_cog] converting to a tiled, pyramidal COG (once):", src)
    savepath <- .gdsrc_allocate_artifact_dir(p, create = TRUE)
    uid <- names(savepath)
    savepath <- unname(savepath)
    # written aside and renamed, so an interrupted conversion is never reused
    tmp <- paste0(savepath, ".writing.tif")
    on.exit(unlink(tmp), add = TRUE)
    # Opened as Giotto opens images (noflip): a TIFF without a geotransform,
    # like MERSCOPE's, is otherwise flipped by terra, and the north-up COG
    # would then hold the source's rows reversed.
    r <- suppressWarnings(.create_terra_spatraster(src)) # "unknown extent"
    terra::writeRaster(r, tmp, filetype = "COG", datatype = terra::datatype(r)[[1L]],
        NAflag = NA, overwrite = TRUE, progress = 0,
        gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "BLOCKSIZE=512",
                 "NUM_THREADS=ALL_CPUS", "BIGTIFF=IF_SAFER"))
    file.rename(tmp, savepath)
    .gdsrc_json_add_artifact(p,
        store_type = "IMAGE",
        uid = uid,
        hash = key,
        meta = list(source = src, source_key = key)
    )
    savepath
}
