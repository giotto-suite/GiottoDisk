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
#' @returns the modified gobject (invisible). Mutated by the internal
#'   adoption pass — captured by `snapshotSave(gDirSource, giottoMulti)`
#'   so the multi-level `.rds` sees post-adoption file handles in each
#'   child.
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
    .adopt_session_reset()
      
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

    vmsg(.v = verbose, "[GiottoDisk] checking expression values...")
    x <- .ss_gdsrc_register_external_expr(x,
        giottosave = name,
        verbose = verbose
    )

    vmsg(.v = verbose, "[GiottoDisk] checking spatial stores...")
    x <- .ss_gdsrc_register_external_geom(x,
        giottosave = name,
        verbose = verbose
    )

    vmsg(.v = verbose, "[GiottoDisk] checking overlap stores...")
    x <- .ss_gdsrc_register_external_overlap(x,
        giottosave = name,
        verbose = verbose
    )

    vmsg(.v = verbose, "[GiottoDisk] checking network stores...")
    x <- .ss_gdsrc_register_external_network(x,
        giottosave = name,
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
    if (length(uids) > 0L) {
        manifest <- as.data.frame(src)
        for (uid_to_tag in uids) {
            content <- manifest[uid == uid_to_tag, giottosave]
            content <- c(content, name)
            content <- unique(content[!is.na(content)])
            src[uid_to_tag, "giottosave"] <- content
        }
    }

    vmsg(.v = verbose, "[GiottoDisk] done")
    invisible(x)
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

.ss_gdsrc_register_external_expr <- function(gobject, giottosave, verbose = NULL) {
    src <- gobject@source
    data_list <- gobject[["expression"]]
    for (x_i in data_list) {
        mat <- x_i[]
        if (!inherits(mat, "IterableMatrix")) next
        if (sourceContains(src, mat)) next

        vmsg(.v = verbose, sprintf(
            "[GiottoDisk] adopting external matrix '%s'",
            GiottoClass::objName(x_i)
        ))

        mat <- sourceAdopt(src, mat, giottosave = giottosave)
        x_i[] <- mat
        gobject <- GiottoClass::setGiotto(gobject, x_i, verbose = FALSE)
    }
    gobject
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
    # overlaps
    uids <- c(uids, .ss_gdsrc_detect_uid_overlaps(gobject))
    # networks (nn + spatial)
    uids <- c(uids, .ss_gdsrc_detect_uid_networks(gobject))
    uids
}

# detect spatial store uids — uses .ss_store_uids so union stores get
# expanded to their substore uids (a union has no @uid of its own).
.ss_gdsrc_detect_uid_spatial_points <- function(gobject) {
    pts_list <- gobject[["feat_info"]]
    is_tracked_class <- vapply(pts_list,
        function(x) inherits(x[], "dataStore"),
        FUN.VALUE = logical(1L)
    )
    pts_list <- pts_list[is_tracked_class]
    if (length(pts_list) == 0L) return(c())
    unlist(lapply(pts_list, function(x) .ss_store_uids(x[])))
}
.ss_gdsrc_detect_uid_spatial_polygons <- function(gobject) {
    polys_list <- gobject[["spatial_info"]]
    is_tracked_class <- vapply(polys_list,
        function(x) inherits(x[], "dataStore"),
        FUN.VALUE = logical(1L)
    )
    polys_list <- polys_list[is_tracked_class]
    if (length(polys_list) == 0L) return(c())
    unlist(lapply(polys_list, function(x) .ss_store_uids(x[])))
}

.ss_gdsrc_detect_uid_overlaps <- function(gobject) {
    poly_list <- gobject[["spatial_info"]]
    uids <- c()
    for (poly_obj in poly_list) {
        ovlps <- poly_obj@overlaps
        if (is.null(ovlps)) next
        for (feat_name in names(ovlps)) {
            if (feat_name == "intensity") next
            ovlp <- ovlps[[feat_name]]
            if (!inherits(ovlp, "overlapPointDisk")) next
            uids <- c(uids, ovlp@data@uid)
        }
    }
    uids
}

# detect network store uids — networks land in @network on nnNetObj /
# spatialNetworkObj. In-memory igraphs are skipped; only dataStore-backed
# networks (typically parquetEdgeStore after a setter auto-write) carry
# vault uids and need snapshot tagging to survive sourcePrune.
.ss_gdsrc_detect_uid_networks <- function(gobject) {
    uids <- c()
    for (slot in c("nn_network", "spatial_network")) {
        net_list <- gobject[[slot]]
        for (net_obj in net_list) {
            net <- net_obj@network
            if (!inherits(net, "dataStore")) next
            uids <- c(uids, .ss_store_uids(net))
        }
    }
    uids
}

.ss_gdsrc_register_external_geom <- function(gobject, giottosave, verbose = NULL) {
    src <- gobject@source

    pts_list <- gobject[["feat_info"]]
    for (pts_obj in pts_list) {
        store <- pts_obj[]
        if (!inherits(store, "dataStore")) next
        if (sourceContains(src, store)) next

        vmsg(.v = verbose, sprintf(
            "[GiottoDisk] adopting external points '%s'",
            GiottoClass::objName(pts_obj)
        ))

        store <- sourceAdopt(src, store, giottosave = giottosave)
        pts_obj[] <- store
        gobject <- GiottoClass::setGiotto(gobject, pts_obj, verbose = FALSE)
    }

    poly_list <- gobject[["spatial_info"]]
    for (poly_obj in poly_list) {
        store <- poly_obj[]
        if (!inherits(store, "dataStore")) next
        if (sourceContains(src, store)) next

        vmsg(.v = verbose, sprintf(
            "[GiottoDisk] adopting external polygons '%s'",
            GiottoClass::objName(poly_obj)
        ))

        store <- sourceAdopt(src, store, giottosave = giottosave)
        poly_obj[] <- store
        gobject <- GiottoClass::setGiotto(gobject, poly_obj, verbose = FALSE)
    }

    gobject
}

.ss_gdsrc_register_external_overlap <- function(gobject, giottosave, verbose = NULL) {
    src <- gobject@source

    poly_list <- gobject[["spatial_info"]]
    for (poly_obj in poly_list) {
        ovlps <- poly_obj@overlaps
        if (is.null(ovlps)) next
        modified <- FALSE

        for (feat_name in names(ovlps)) {
            if (feat_name == "intensity") next
            ovlp <- ovlps[[feat_name]]
            if (!inherits(ovlp, "overlapPointDisk")) next
            if (sourceContains(src, ovlp@data)) next

            vmsg(.v = verbose, sprintf(
                "[GiottoDisk] adopting external overlap '%s/%s'",
                GiottoClass::objName(poly_obj), feat_name
            ))

            depends <- c(ovlp@poly_uids, ovlp@feat_uids)
            ovlp@data <- sourceAdopt(src, ovlp@data,
                giottosave = giottosave, depends = depends
            )
            poly_obj@overlaps[[feat_name]] <- ovlp
            modified <- TRUE
        }

        if (modified) {
            gobject <- GiottoClass::setGiotto(gobject, poly_obj, verbose = FALSE)
        }
    }

    gobject
}

# Mirror of the geom/overlap registrants for networks. The common path
# is: user calls setNearestNetwork / setSpatialNetwork with an in-memory
# igraph on a backed gobject → setter auto-writes via
# sourceWrite(gDirSource, igraph) → network is already vault-resident
# before snapshotSave runs. This helper covers the rarer case where a
# parquetEdgeStore was attached externally (e.g. constructed manually,
# or moved between projects). Walks both @nn_network and
# @spatial_network slots uniformly.
.ss_gdsrc_register_external_network <- function(gobject, giottosave, verbose = NULL) {
    src <- gobject@source

    for (slot in c("nn_network", "spatial_network")) {
        net_list <- gobject[[slot]]
        for (net_obj in net_list) {
            net <- net_obj@network
            if (!inherits(net, "dataStore")) next
            if (sourceContains(src, net)) next

            vmsg(.v = verbose, sprintf(
                "[GiottoDisk] adopting external network '%s'",
                GiottoClass::objName(net_obj)
            ))

            net_obj@network <- sourceAdopt(src, net, giottosave = giottosave)
            gobject <- GiottoClass::setGiotto(gobject, net_obj, verbose = FALSE)
        }
    }

    gobject
}

# canonical hash(es) of the underlying storage, stripping lazy ops /
# view state. Returns a character vector — one entry per leaf / substore
# so each manifest entry can be matched independently. This handles
# compound stores (multi-leaf IterableMatrix, union*Store) without the
# concat-collapse trick that never matched any single manifest hash.
.ss_hash_expr_base <- function(mat) {
    if (inherits(mat, "IterableMatrix")) {
        vapply(.im_leaf_dirs(mat),
            function(d) .hash(BPCells::open_matrix_dir(d)),
            FUN.VALUE = character(1L)
        )
    } else if (inherits(mat, "HDF5Array")) {
        .hash(HDF5Array::HDF5Array(HDF5Array::path(mat), HDF5Array::name(mat)))
    } else if (inherits(mat, "unionParquetStore") ||
               inherits(mat, "unionParquetExprStore")) {
        unlist(lapply(mat@stores, .ss_hash_expr_base))
    } else if (inherits(mat, "dataStore")) {
        .hash(storeRead(.store_nostate(mat)))
    } else {
        .hash(mat)
    }
}

# Returns the uids of all leaf / substore artifacts inside a (possibly
# compound) store. Used by the spatial-store detectors where direct
# @uid access fails for union stores (which have @stores, not @uid).
.ss_store_uids <- function(x) {
    if (inherits(x, "unionParquetStore") ||
        inherits(x, "unionParquetExprStore")) {
        return(unlist(lapply(x@stores, .ss_store_uids)))
    }
    if (inherits(x, "fileStore")) return(x@uid)
    character(0L)
}

# detect matrix uid based on hash lookup in manifest
.ss_gdsrc_detect_uid_matrices <- function(gobject, manifest) {
    tracked_classes <- c("IterableMatrix", "HDF5Array",
                         "parquetExprStore", "unionParquetExprStore")
    mat_list <- gobject[["expression"]]

    is_tracked_class <- vapply(mat_list,
        function(x) inherits(x[], tracked_classes),
        FUN.VALUE = logical(1L)
    )
    mat_list <- mat_list[is_tracked_class]
    if (length(mat_list) == 0L) return(c())

    # unlist+lapply (not vapply) because compound stores return multiple
    # hashes per matrix.
    protected_hash <- unlist(lapply(mat_list,
        function(x) .ss_hash_expr_base(x[])
    ))
    manifest <- manifest[!duplicated(hash)]
    manifest[hash %in% protected_hash, uid]
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
