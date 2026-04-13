#' @name snapshotSave
#' @title Write Giotto Snapshot
#' @description
#' Register a gobject snapshot into the source
#' @param src `gsource` object. Used to provide default save locations
#' and save formats for a managed Giotto backend.
#' @param x `giotto` object to write
#' @param method `character` (default = "rds").
#' @param method_params `list` (optional) list of additional named params for
#'   save method.
#' @param overwrite `logical` (default = FALSE). Whether to overwrite if a
#'   snapshot of the same name exists.
#' @param verbose verbosity
#' @param ... additional params to pass (none implemented)
#' @returns TRUE if save completed
NULL

#' @rdname snapshotSave
#' @param export_image `logical` (default = TRUE) Whether to make a copy of
#'   external images in the project directory
#' @export
setMethod("snapshotSave", signature("gDirSource", "giotto"), function(src, x,
    name = format(Sys.time(), "%Y%m%d_giottosave"),
    method = c("rds", "qs"),
    method_params = list(),
    overwrite = FALSE,
    export_image = TRUE,
    verbose = NULL,
    ...) {
    method <- match.arg(tolower(method), choices = c("rds", "qs"))
    checkmate::assert_list(method_params)
    checkmate::assert_flag(overwrite)
      
    p <- x@source@path
    vmsg(.v = verbose, "[GiottoDisk] Saving giotto...")
    print_list(list(
        backend = sprintf("`%s`", class(src)),
        path = p
    ), pre = "  ")
    gsdir <- .gdsrc_giottosave_dir(p)
    if (!dir.exists(gsdir)) dir.create(gsdir, recursive = TRUE)
    if (!methods::hasArg(name)) {
        name <- .ss_gdsrc_find_autoname(gsdir, name)
    }
    if (name %in% .gdsrc_detect_gsavename(gsdir) && !overwrite) {
       stop("[snapshotSave] giotto snapshot of this name already exists.\n",
          "Change 'name' or set 'overwrite = TRUE'", call. = FALSE)
    }
  
    vmsg(.v = verbose, "[GiottoDisk] checking images...")
    x <- .ss_gdsrc_register_external_images(x, 
        giottosave = name,
        export_image = export_image,
        verbose = verbose
    )
  
    vmsg(.v = verbose, "[GiottoDisk] writing snapshot")
    temp <- .dump_tempfile() # temp location for atomic writes
    switch(method,
        "rds" = {
            fullpath <- file.path(gsdir, paste0(name, ".rds"))
            a <- c(list(object = x, file = temp), method_params)
            do.call(saveRDS, a)
        },
        "qs" = {
            package_check(pkg_name = "qs", repository = "CRAN")
            qsave_fun <- get("qsave", asNamespace("qs"))
            fullpath <- file.path(gsdir, paste0(name, ".qs"))
            a <- c(list(x = x, file = temp), method_params)
            do.call(qsave_fun, a)
        }
    )
    
    res <- file.rename(from = temp, to = fullpath)
    if (!res) stop("[snapshotSave] save failed\n", call. = FALSE)
  
    # tagging --------------------------------------------------- #
    vmsg(.v = verbose, "[GiottoDisk] tagging snapshot artifacts...")
    uids <- .ss_gdsrc_detect_uid(x)
    if (length(uids) == 0L) return(invisible(TRUE))
    
    manifest <- as.data.frame(src)
    for (uid_to_tag in uids) {
        content <- manifest[uid == uid_to_tag, giottosave]
        content <- c(content, name)
        content <- unique(content[!is.na(content)])
        src[uid_to_tag, "giottosave"] <- content
    }
  
    vmsg(.v = verbose, "[GiottoDisk] done")
    invisible(TRUE)
})

# internals ####

# copied from GiottoClass
.create_terra_spatraster <- function(image_path) {
    raster_object <- try(terra::rast(x = image_path, noflip = TRUE))
    if (inherits(raster_object, "try-error")) {
        stop(raster_object, " can not be read by terra::rast() \n")
    }
    return(raster_object)
}

.ss_gdsrc_register_external_images <- function(gobject, giottosave,
    export_image = TRUE,
    verbose = NULL) {
    src <- gobject@source
    img_list <- gobject[["images"]]
    for (img in img_list) {
        r <- img@raster_object
        if (is.null(r)) next
        if (sourceContains(src, r)) next

        f <- normalizePath(terra::sources(r), mustWork = FALSE)
        if (!export_image && nzchar(f)) next # skip external on-disk if not exporting

        vmsg(.v = verbose, sprintf(
            "[GiottoDisk] processing external image '%s'",
            GiottoClass::objName(img)
        ))

        # save extent info (needed for non-COG formats)
        img@extent <- terra::ext(r)[]

        r <- sourceAdopt(src, r, giottosave = giottosave)

        img@file_path <- terra::sources(r)
        gobject <- GiottoClass::setGiotto(gobject, img, verbose = FALSE)
    }
    gobject
}

.ss_gdsrc_detect_uid <- function(gobject) {
    manifest <- as.data.frame(gobject@source)
    uids <- c()
    # matrices
    uids <- c(uids, .ss_gdsrc_detect_uid_matrices(gobject, manifest))
    # points
    uids <- c(uids, .ss_gdsrc_detect_uid_spatial_points(gobject))
    # polys
    uids <- c(uids, .ss_gdsrc_detect_uid_spatial_polygons(gobject))
    uids
}

# detect spatial store uids
.ss_gdsrc_detect_uid_spatial_points <- function(gobject) {
    pts_list <- gobject[["feat_info"]]
    is_tracked_class <- vapply(pts_list,
        function(x) inherits(x[], "dataStore"),
        FUN.VALUE = logical(1L)
    )
    pts_list <- pts_list[is_tracked_class]
    if (length(pts_list) == 0L) return(c())
  
    vapply(pts_list,
        function(x) x[]@uid,
        FUN.VALUE = character(1L)
    )
}
.ss_gdsrc_detect_uid_spatial_polygons <- function(gobject) {
    polys_list <- gobject[["spatial_info"]]
    is_tracked_class <- vapply(polys_list,
        function(x) inherits(x[], "dataStore"),
        FUN.VALUE = logical(1L)
    )
    polys_list <- polys_list[is_tracked_class]
    if (length(polys_list) == 0L) return(c())
  
    vapply(polys_list,
        function(x) x[]@uid,
        FUN.VALUE = character(1L)
    )
}

# detect matrix uid based on hash lookup in manifest
.ss_gdsrc_detect_uid_matrices <- function(gobject, manifest) {
    tracked_classes <- c("IterableMatrix", "DelayedArray")
    mat_list <- gobject[["expression"]]
  
    is_tracked_class <- vapply(mat_list,
        function(x) inherits(x[], tracked_classes),
        FUN.VALUE = logical(1L)
    )
    mat_list <- mat_list[is_tracked_class]
    if (length(mat_list) == 0L) return(c())
  
    protected_hash <- vapply(mat_list,
        function(x) .hash(x[]), 
        FUN.VALUE = character(1L)
    )
    manifest <- manifest[!duplicated(hash)]
    manifest[hash %in% protected_hash, uid] # get uids matching hashes
}

# add numerical suffix to prevent file naming collision
# p is the gsave subdirectory
# name is the proposed gsave name
.ss_gdsrc_find_autoname <- function(p, name) {
    ntest <- name
    counter <- 1L
    existing_names <- .gdsrc_detect_gsavename(p)
    while (ntest %in% existing_names) {
        ntest <- sprintf("%s_%02d", name, counter)
        counter <- counter + 1L
        if (counter > 100L) {
            stop("[snapshotSave] automatic 'name' setting failed\n", call. = FALSE)
        }
    }
    ntest
}
