# Xenium ingest pipeline for `gDirSource`-managed projects.
#
# Parallels Giotto's in-memory `XeniumReader`, but each `load_*` call
# routes through `sourceWrite()` so the produced parquet stores live in
# the user's project vault. Path detection is reused from Giotto via
# the parent's `@paths` slot, so this file owns only the disk-side
# helpers and the subclass init that swaps in disk-flavored closures.



# CLASS ####



setClass(
    "XeniumDiskReader",
    contains = "XeniumReader",
    slots = list(
        backend = "ANY"
    ),
    prototype = list(
        backend = NULL
    )
)

# * init ####
setMethod(
    "initialize", signature("XeniumDiskReader"),
    function(.Object, ..., backend) {
        obj <- callNextMethod(.Object, ...)

        if (!missing(backend)) {
            # mirror createGiottoObject(backend = ...) handling: a
            # character input is treated as a directory path and coerced
            # to a gDirSource.
            if (is.character(backend)) {
                backend <- gDirSource(path = backend)
            }
            checkmate::assert_class(backend, "gsource")
            obj@backend <- backend
        }
        if (is.null(obj@backend)) {
            stop("[XeniumDiskReader] `backend` is required", call. = FALSE)
        }

        # Mirror @paths into the init frame so subclass closures can
        # reference path names directly via default-arg expressions
        # (same convention as the parent init). Each default evaluates
        # to a bare character string at call time, so nothing scope-
        # laddering back to this frame travels into parallel workers.
        list2env(obj@paths, envir = environment())
        gsrc <- obj@backend

        # transcripts (disk override)
        tx_fun <- function(
            path = tx_path,
            feat_type = c(
                "rna",
                "NegControlProbe",
                "UnassignedCodeword",
                "NegControlCodeword",
                "GenomicControl"
            ),
            split_keyword = list(
                "NegControlProbe",
                "UnassignedCodeword",
                "NegControlCodeword",
                "GenomicControl"
            ),
            flip_vertical = TRUE,
            dropcols = c(),
            qv_threshold = obj@qv,
            cores = GiottoUtils::determine_cores(),
            output = c("giottoPoints", "store"),
            verbose = NULL,
            ...
        ) {
            .xenium_transcript_disk(
                path = path,
                gsource = gsrc,
                feat_type = feat_type,
                split_keyword = split_keyword,
                flip_vertical = flip_vertical,
                dropcols = dropcols,
                qv_threshold = qv_threshold,
                cores = cores,
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_transcripts <- tx_fun

        # polygons (disk override). Signature mirrors parent's `load_polys`
        # so the inherited `create_gobject` can pass the same args.
        poly_fun <- function(
            path = cell_bound_path,
            name = "cell",
            part_col = NULL,
            split_geom = FALSE,
            split_geom_fmt = "%d",
            split_geom_sourcename = "cell_poly_id",
            flip_vertical = TRUE,
            calc_centroids = TRUE,
            cores = GiottoUtils::determine_cores(),
            output = c("giottoPolygon", "store"),
            verbose = NULL,
            ...
        ) {
            .xenium_poly_disk(
                path = path,
                gsource = gsrc,
                name = name,
                part_col = part_col,
                split_geom = split_geom,
                split_geom_fmt = split_geom_fmt,
                split_geom_sourcename = split_geom_sourcename,
                flip_vertical = flip_vertical,
                calc_centroids = calc_centroids,
                cores = cores,
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_polys <- poly_fun

        # expression (disk override). Signature mirrors parent's
        # `load_expression` so the inherited `gobject_fun` (and the
        # disk reader's own `create_gobject`) can pass the same args.
        ex_fun <- function(
            path = expr_path,
            gene_ids = "symbols",
            remove_zero_rows = TRUE,
            split_by_type = TRUE,
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            .xenium_expression_disk(
                path = path,
                gsource = gsrc,
                gene_ids = gene_ids,
                remove_zero_rows = remove_zero_rows,
                split_by_type = split_by_type,
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_expression <- ex_fun

        # create_gobject (disk variant). Mirrors parent's gobject_fun
        # signature so `createGiottoXeniumObject(backend = gsrc, ...)` can
        # plumb identical args through. Image / aligned-image attach
        # logic is duplicated from parent's gobject_fun so the on-disk vs
        # in-mem boundary can shift independently without coordinating
        # changes across both packages.
        gobject_fun <- function(
            transcript_path = tx_path,
            load_bounds = list(
                cell = "cell",
                nucleus = "nucleus"
            ),
            gene_panel_json_path = panel_meta_path,
            expression_path = expr_path,
            metadata_path = cell_meta_path,
            feat_type = c(
                "rna",
                "NegControlProbe",
                "UnassignedCodeword",
                "NegControlCodeword",
                "GenomicControl"
            ),
            split_keyword = list(
                "NegControlProbe",
                "UnassignedCodeword",
                "NegControlCodeword",
                "GenomicControl"
            ),
            load_images = "focus",
            load_aligned_images = NULL,
            load_transcripts = TRUE,
            load_expression = TRUE,
            load_cellmeta = FALSE,
            instructions = NULL,
            verbose = NULL) {
            load_transcripts <- as.logical(load_transcripts)
            load_expression <- as.logical(load_expression)
            load_cellmeta <- as.logical(load_cellmeta)

            if (!load_transcripts && !load_expression) {
                warning(GiottoUtils::wrap_txt(
                    "One of either transcripts or expression info should be loaded for a fully functioning object"
                ))
            }

            if (!is.null(load_aligned_images)) {
                checkmate::assert_list(load_aligned_images)
                if (is.null(names(load_aligned_images))) {
                    stop(GiottoUtils::wrap_txt(
                        "'load_aligned_images' must be a named list"
                    ))
                }
                if (any(lengths(load_aligned_images) != 2L) ||
                    any(!vapply(load_aligned_images, is.character,
                        FUN.VALUE = logical(1L)
                    ))) {
                    stop(GiottoUtils::wrap_txt(
                        "'load_aligned_images' must be character with length 2:
                        1. image path
                        2. alignment matrix path"
                    ))
                }
            }
            if (!is.null(load_bounds)) {
                checkmate::assert_list(load_bounds)
                if (is.null(names(load_bounds))) {
                    stop("'load_bounds' must be named list of filepaths\n")
                }
            }

            funs <- obj@calls

            # init gobject with disk backend
            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            # transcripts (disk)
            if (load_transcripts) {
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, tx_list, verbose = FALSE)
            }

            # polys (disk)
            if (!is.null(load_bounds)) {
                load_bounds[load_bounds == "cell"] <- cell_bound_path
                load_bounds[load_bounds == "nucleus"] <- nuc_bound_path

                blist <- list()
                bnames <- names(load_bounds)
                for (b_i in seq_along(load_bounds)) {
                    b_name <- bnames[[b_i]]
                    if (b_name == "nucleus") {
                        b <- funs$load_polys(
                            path = load_bounds[[b_i]],
                            name = b_name,
                            part_col = NULL,
                            split_geom = TRUE,
                            split_geom_fmt = "nucleus_%d",
                            split_geom_sourcename = "cell_poly_id",
                            verbose = verbose
                        )
                    } else {
                        b <- funs$load_polys(
                            path = load_bounds[[b_i]],
                            name = b_name,
                            part_col = NULL,
                            split_geom = FALSE,
                            verbose = verbose
                        )
                    }
                    blist <- c(blist, b)
                }
                # centroids_to_spatlocs = TRUE attaches centroids as
                # spatlocs at set time, computing them lazily from the
                # parquetGeomBase. Replaces the explicit
                # addSpatialCentroidLocations step at the end of the
                # in-mem reader's gobject_fun.
                g <- GiottoClass::setGiotto(g, blist,
                    centroids_to_spatlocs = TRUE,
                    verbose = FALSE)
            }

            # expression (disk; overridden closure)
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    gene_ids = "symbols",
                    remove_zero_rows = TRUE,
                    split_by_type = TRUE,
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, ex)
            }

            # feat metadata (in-mem; inherited closure)
            # optional: some Xenium-format exports ship no panel json, and
            # feature metadata is generated from the expression matrix when
            # it is absent. Mirrors the guard in Giotto's create path.
            if (length(gene_panel_json_path) > 0L &&
                nzchar(gene_panel_json_path[[1L]])) {
                fx <- funs$load_featmeta(
                    path = gene_panel_json_path,
                    gene_ids = "symbols",
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, fx, verbose = FALSE)
            }

            # cell metadata (in-mem; inherited closure)
            if (load_cellmeta) {
                cx <- funs$load_cellmeta(
                    path = metadata_path,
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, cx, verbose = FALSE)
            }

            # images (inherited closures; orchestration duplicated from
            # parent's gobject_fun)
            if (!is.null(load_images)) {
                load_images <- lapply(load_images, normalizePath, mustWork = FALSE)
                img_focus_path <- normalizePath(img_focus_path, mustWork = FALSE)

                # replace shortname
                load_images[load_images == "focus"] <- img_focus_path

                is_dir <- dir.exists(img_focus_path)
                is_focus <- load_images == img_focus_path
                is_focus_image <- is_focus & !is_dir
                is_focus_dir <- is_focus & is_dir

                # handle matches to single focus images instead of a directory
                names(load_images)[is_focus_image] <- "dapi"

                # [exception] handle focus image dir
                if (any(is_focus_dir)) {
                    # split the focus image dir away from other entries
                    load_images <- load_images[!is_focus_dir]
                    focus_dir <- img_focus_path
                    focus_files <- list.files(focus_dir, pattern = "\\.ome\\.tif$", full.names = TRUE)
                    # fallback: if no OME-TIFFs, try plain TIFs in the same folder
                    if (length(focus_files) == 0L) {
                        focus_files <- list.files(focus_dir, pattern = "\\.tif$", full.names = TRUE)
                        focus_files <- focus_files[!grepl("\\.ome\\.tif$", focus_files, ignore.case = TRUE)]
                        focus_files <- focus_files[!dir.exists(focus_files)]
                    }
                    # ignore matches to export dir (if it is a subdirectory)
                    focus_files <- focus_files[!dir.exists(focus_files)]
                    if (length(focus_files) > 0L) {
                        bn <- basename(focus_files)
                        marker <- ifelse(grepl("^ch\\d+_.*\\.ome\\.tif$", bn, ignore.case = TRUE),
                                         sub("^ch\\d+_(.*?)\\.ome\\.tif$", "\\1", bn, ignore.case = TRUE),
                                         NA_character_)
                        is_dapi <- grepl("dapi", bn, ignore.case = TRUE) |
                            (!is.na(marker) & grepl("^dapi$", marker, ignore.case = TRUE))

                        # put DAPI first if present
                        ord <- c(which(is_dapi), which(!is_dapi))
                        focus_files <- focus_files[ord]
                        marker <- marker[ord]
                        is_dapi <- is_dapi[ord]
                        fallback_idx <- which(!is_dapi & is.na(marker))
                        names_vec <- character(length(focus_files))
                        names_vec[is_dapi] <- "dapi"
                        names_vec[!is_dapi & !is.na(marker)] <- marker[!is_dapi & !is.na(marker)]
                        if (length(fallback_idx)) {
                            names_vec[fallback_idx] <- sprintf("bound%d", seq_along(fallback_idx))
                        }
                        names(focus_files) <- names_vec

                        # append to rest of entries
                        load_images <- c(load_images, focus_files)
                    }
                }

                # ensure that input is list
                checkmate::assert_list(load_images)
                if (is.null(names(load_images))) {
                    stop("'load_images' must be a named list of filepaths\n")
                }

                imglist <- list()
                imnames <- names(load_images)
                for (impath_i in seq_along(load_images)) {
                    im <- funs$load_image(
                        path = load_images[[impath_i]],
                        name = imnames[[impath_i]]
                    )
                    imglist <- c(imglist, im)
                }
                g <- GiottoClass::setGiotto(g, imglist)
            }

            # aligned images can be placed in random places and do not have
            # a standardized naming scheme. Cannot load with expected default.
            if (!is.null(load_aligned_images)) {
                aimglist <- list()
                aimnames <- names(load_aligned_images)
                for (aim_i in seq_along(load_aligned_images)) {
                    GiottoUtils::vmsg(.v = verbose, sprintf(
                        "loading aligned image as '%s'",
                        aimnames[[aim_i]]
                    ))
                    aim <- funs$load_aligned_image(
                        path = load_aligned_images[[aim_i]][1],
                        imagealignment_path = load_aligned_images[[aim_i]][2],
                        name = aimnames[[aim_i]]
                    )
                    aimglist <- c(aimglist, aim)
                }
                g <- GiottoClass::setGiotto(g, aimglist)
            }

            # centroids are attached as spatlocs at setGiotto time
            # via `centroids_to_spatlocs = TRUE` (see polys block above).
            # With no polygons loaded there are none, so fall back to the
            # centroids recorded in the cell metadata file. Coordinates stay
            # in-memory here, as cellmeta and featmeta already do.
            if (length(GiottoClass::list_spatial_info_names(g)) == 0L) {
                sl <- funs$load_spatlocs(verbose = verbose)
                if (!is.null(sl)) {
                    g <- GiottoClass::setGiotto(g, sl, verbose = FALSE)
                }
            }

            GiottoUtils::vmsg(.v = verbose, "done")

            return(g)
        }
        obj@calls$create_gobject <- gobject_fun

        obj
    }
)



# CREATE READER ####

#' @title Import a 10X Xenium assay (disk-backed)
#' @name importXeniumDisk
#' @description
#' Disk-backed counterpart to [Giotto::importXenium()]. Produces a
#' `XeniumDiskReader` whose `load_transcripts()` and `load_polys()` calls
#' write to a `gDirSource`-managed project vault as `parquetGeomTile`
#' stores. Other modalities (expression, featmeta, cellmeta, images)
#' remain in-memory via the inherited `XeniumReader` closures.
#' @param xenium_dir Xenium output directory
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#'   Naming matches [GiottoClass::createGiottoObject()]'s `backend` param.
#' @param qv_threshold minimum Phred-scaled quality score retained
#' @returns `XeniumDiskReader` object
#' @seealso [Giotto::importXenium()] for the in-memory variant
#' @export
importXeniumDisk <- function(xenium_dir = NULL, backend, qv_threshold = 20) {
    if (missing(backend)) {
        stop("[importXeniumDisk] `backend` is required", call. = FALSE)
    }
    a <- list(Class = "XeniumDiskReader", backend = backend, qv = qv_threshold)
    if (!is.null(xenium_dir)) a$xenium_dir <- xenium_dir
    do.call(new, args = a)
}



# MODULAR ####



## transcript ####

# Disk-backed Xenium transcript ingestion.
#
# Streams the source parquet through arrow with a patched schema (string
# casts for `transcript_id`, `cell_id`, `feature_name`) plus lazy qv
# filter and optional y flip, then writes directly to a
# `parquetGeomTileStore`. No `parquetStore` intermediate is needed since
# points don't require an explicit `row_index`.
.xenium_transcript_disk <- function(
    path,
    gsource,
    feat_type = c(
        "rna",
        "NegControlProbe",
        "UnassignedCodeword",
        "NegControlCodeword",
        "GenomicControl"
    ),
    split_keyword = list(
        "NegControlProbe",
        "UnassignedCodeword",
        "NegControlCodeword",
        "GenomicControl"
    ),
    flip_vertical = TRUE,
    dropcols = c(),
    qv_threshold = 20,
    cores = GiottoUtils::determine_cores(),
    output = c("giottoPoints", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path)) {
        stop("[xenium_transcript_disk] no path provided", call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("giottoPoints", "store"))
    GiottoUtils::package_check("arrow")
    GiottoUtils::package_check("dplyr")

    GiottoUtils::vmsg("[xenium_transcript_disk] loading transcripts", .v = verbose)

    # capture bare scalars only — read_fun must be safe to ship to workers
    qv_thr <- qv_threshold
    flip <- isTRUE(flip_vertical)
    drop <- as.character(dropcols)

    read_fun <- function(x, schema = NULL) {
        if (is.null(schema)) {
            sc <- arrow::schema(arrow::open_dataset(x))
            cast_cols <- intersect(
                c("transcript_id", "cell_id", "feature_name"),
                names(sc)
            )
            new_fields <- lapply(names(sc), function(nm) {
                if (nm %in% cast_cols) {
                    arrow::field(nm, arrow::string())
                } else {
                    sc$GetFieldByName(nm)
                }
            })
            schema <- do.call(arrow::schema, new_fields)
        }
        a <- arrow::open_dataset(sources = x, schema = schema)
        if (length(drop) > 0L) {
            a <- dplyr::select(a, -dplyr::any_of(drop))
        }
        if (!is.null(qv_thr)) {
            a <- dplyr::filter(a, qv >= qv_thr)
        }
        if (flip) {
            a <- dplyr::mutate(a, y_location = -y_location)
        }
        a
    }

    fs <- fileStore(path = path, read_fun = read_fun)
    qs <- as(fs, "queryableStore")

    tx_store <- sourceWrite(
        gsource, qs,
        store_type = "parquetGeomTile",
        type = "points",
        id_col = "feature_name",
        sdimx = "x_location",
        sdimy = "y_location",
        verbose = verbose,
        ...
    )

    if (output == "store") return(tx_store)

    gpoints <- GiottoClass::createGiottoPoints(
        x = tx_store,
        feat_type = feat_type,
        split_keyword = split_keyword
    )
    if (!inherits(gpoints, "list")) gpoints <- list(gpoints)
    gpoints
}



## polygon ####

# Disk-backed Xenium polygon ingestion.
#
# Writes a `parquetStore` intermediate first (polygons need `row_index`
# for stable vertex ordering), then a `parquetGeomTileStore`. The
# optional vertical flip is applied lazily via the source read_fun on
# the way into the intermediate.
.xenium_poly_disk <- function(
    path,
    gsource,
    name = "cell",
    part_col = NULL,
    split_geom = FALSE,
    split_geom_fmt = "%d",
    split_geom_sourcename = "cell_poly_id",
    flip_vertical = TRUE,
    calc_centroids = TRUE,
    cores = GiottoUtils::determine_cores(),
    output = c("giottoPolygon", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path)) {
        stop("[xenium_poly_disk] no path provided", call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    checkmate::assert_class(gsource, "gsource")
    checkmate::assert_character(name, len = 1L)
    output <- match.arg(output, choices = c("giottoPolygon", "store"))
    GiottoUtils::package_check("arrow")
    GiottoUtils::package_check("dplyr")

    GiottoUtils::vmsg(
        sprintf("[xenium_poly_disk] loading boundary '%s'", name),
        .v = verbose
    )

    # No dplyr mutations here — those push arrow into a threaded scan,
    # which delivers batches out of source order. That's harmless for
    # points (each row is independent), but for polygons it scatters a
    # cell's vertices across row_index values, splitting cells into
    # multiple non-contiguous runs in the intermediate parquet and
    # breaking the polygon-construction step. Flip is applied lazily on
    # the final geom store below, after row_index is locked in.
    flip <- isTRUE(flip_vertical)

    read_fun <- function(x, schema = NULL) {
        if (is.null(schema)) {
            sc <- arrow::schema(arrow::open_dataset(x))
            casts <- list(
                cell_id = arrow::string(),
                label_id = arrow::int32()
            )
            casts <- casts[names(casts) %in% names(sc)]
            new_fields <- lapply(names(sc), function(nm) {
                if (nm %in% names(casts)) {
                    arrow::field(nm, casts[[nm]])
                } else {
                    sc$GetFieldByName(nm)
                }
            })
            schema <- do.call(arrow::schema, new_fields)
        }
        arrow::open_dataset(sources = x, schema = schema)
    }

    # auto-detect part_col from source schema if not provided.
    # Same preference order + regex fallback as the in-mem .xenium_poly().
    src_cols <- names(arrow::schema(arrow::open_dataset(path)))
    if (is.null(part_col) || !part_col %in% src_cols) {
        prefer <- c(
            "label_id", "boundary_id",
            "polygon_id", "poly_id", "object_id",
            "part_id", "part",
            "ring_id", "segment_id", "region_id"
        )
        hit <- intersect(prefer, src_cols)
        if (!length(hit)) {
            rx <- "(label|boundary|poly(?:gon)?|object|part|ring|segment|region).*(_?id)?$"
            cand <- setdiff(src_cols, c("cell_id", "vertex_x", "vertex_y", "x", "y", "z"))
            hit <- grep(rx, cand, ignore.case = TRUE, value = TRUE)
        }
        if (length(hit)) {
            GiottoUtils::vmsg(
                sprintf("Using '%s' as polygon part column", hit[[1]]),
                .v = verbose
            )
            part_col <- hit[[1]]
        } else {
            GiottoUtils::vmsg(
                "No polygon part column detected; assuming single part per cell (part_col = NULL)",
                .v = verbose
            )
            part_col <- NULL
        }
    }

    fs <- fileStore(path = path, read_fun = read_fun)
    qs <- as(fs, "queryableStore")

    # parquet intermediate to materialize row_index
    poly_intermediate <- sourceWrite(
        gsource, qs,
        store_type = "parquet",
        verbose = verbose
    )

    # `flip_vertical` is applied via a read_fun wrap on the
    # intermediate inside the tile-write method. This bakes the flip
    # into the on-disk tile values (matching the tx-on-disk
    # convention) so anything reading the raw parquet directly
    # Safe because the intermediate has row_index
    # baked in by the prior `sourceWrite(store_type = "parquet")`.
    poly_store <- sourceWrite(
        gsource, poly_intermediate,
        store_type = "parquetGeomTile",
        type = "polygons",
        id_col = "cell_id",
        sdimx = "vertex_x",
        sdimy = "vertex_y",
        group_col = "cell_id",
        part_col = part_col,
        split_geom = split_geom,
        split_geom_fmt = split_geom_fmt,
        split_geom_sourcename = split_geom_sourcename,
        flip_vertical = flip,
        verbose = verbose,
        ...
    )

    if (output == "store") return(poly_store)

    # Centroids are intentionally not pre-computed on the gpoly here.
    # `setGiotto(g, gpoly, centroids_to_spatlocs = TRUE)` (used by
    # the disk reader's `create_gobject`) computes them lazily from the
    # `parquetGeomBase` at attach time. The `calc_centroids` arg is
    # retained for API parity with the in-mem reader.
    GiottoClass::createGiottoPolygon(
        x = poly_store,
        name = name
    )
}



## expression ####

# Disk-backed Xenium expression ingestion.
#
# Builds an `exprInput` marker (`mtxInput` for the 10x mtx triple,
# `tenxH5Input` for cell_feature_matrix.h5) and routes it through
# `sourceWrite(gsource, inp, store_type = "parquetExpr")` into the
# project vault. tar.gz inputs are unpacked under `tempdir()` and the
# resulting cell_feature_matrix/ directory feeds the mtx path.
#
# Xenium-specific post-processing matches Giotto's in-mem
# `.xenium_expression`: drop zero-detection features (Arrow distinct on
# col_id), split by feature class (features.tsv column 3, or
# /features/feature_type on h5), and rename to Giotto feat_type
# conventions ("Gene Expression" -> "rna", etc).
.xenium_expression_disk <- function(
    path,
    gsource,
    gene_ids = "symbols",
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    output = c("exprObj", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path)) {
        stop("[xenium_expression_disk] no path provided", call. = FALSE)
    }
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("exprObj", "store"))

    feature_id_col <- switch(gene_ids,
        "ensembl" = 1L,
        "symbols" = 2L,
        stop("[xenium_expression_disk] unknown gene_ids: ", gene_ids,
             call. = FALSE)
    )

    # Detect format
    if (dir.exists(path)) {
        fmt <- "mtx"
    } else if (grepl("\\.tar\\.gz$", path, ignore.case = TRUE)) {
        fmt <- "tar.gz"
    } else {
        ext <- tolower(tools::file_ext(path))
        fmt <- switch(ext,
            "h5"  = "h5",
            "mtx" = "mtx",
            "gz"  = "mtx",   # matrix.mtx.gz inside a 10x dir is handled below
            stop("[xenium_expression_disk] unsupported expression format: ",
                 path, call. = FALSE)
        )
    }

    GiottoUtils::vmsg("[xenium_expression_disk] format:", fmt, .v = verbose)

    # tar.gz -> unpack to tempdir, fall through to mtx path
    if (fmt == "tar.gz") {
        unpack_root <- file.path(
            tempdir(),
            paste0("xenium_expr_unpacked_",
                   tools::file_path_sans_ext(
                       tools::file_path_sans_ext(basename(path))))
        )
        dir.create(unpack_root, recursive = TRUE, showWarnings = FALSE)
        utils::untar(path, exdir = unpack_root)
        candidates <- list.dirs(unpack_root, recursive = FALSE)
        unpacked <- candidates[grepl("cell_feature_matrix$", candidates)][1L]
        if (is.na(unpacked) || !dir.exists(unpacked)) {
            stop("[xenium_expression_disk] tarball did not contain a ",
                 "cell_feature_matrix/ directory.", call. = FALSE)
        }
        path <- unpacked
        fmt  <- "mtx"
    }

    # Construct exprInput marker + capture feat_classes (col 3 of features.tsv
    # for mtx; /features/feature_type for h5). Must align 1:1 with inp@feat_ids
    # in the original (pre-disambiguation) order; disambiguation only appends
    # suffixes so positional alignment is preserved.
    if (fmt == "mtx") {
        inp <- mtxInput(path, feature_id_col = feature_id_col)
        feat_classes_vec <- .tenx_feat_classes_mtx(path)
    } else {
        inp <- tenxH5Input(path, feature_id_col = feature_id_col)
        feat_classes_vec <- .tenx_feat_classes_h5(path)
    }

    # Write to vault. sourceWrite(gDirSource, fileStore) is the dispatch
    # entry; exprInput extends fileStore so it routes here. storeWrite
    # internally consumes the iterator returned by storeRead(inp).
    pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                      verbose = verbose, ...)

    # remove_zero_rows: streaming Arrow distinct on col_id to find
    # detected features; subset the store positionally.
    if (isTRUE(remove_zero_rows)) {
        col_id <- NULL  # NSE binding
        ds <- storeRead(pe)
        present <- ds |>
            dplyr::distinct(col_id) |>
            dplyr::collect()
        keep_idx <- sort(as.integer(present$col_id))
        if (length(keep_idx) > 0L &&
            length(keep_idx) < length(feat_classes_vec)) {
            pe <- pe[keep_idx, , drop = FALSE]
            feat_classes_vec <- feat_classes_vec[keep_idx]
        }
    }

    # split_by_type: group feat indices by feature class. If only one
    # class is present, return a single-element list.
    feat_class_to_name <- function(fc) {
        switch(fc,
            "Gene Expression"           = "rna",
            "Protein Expression"        = "protein",
            "Negative Control Codeword" = "NegControlCodeword",
            "Negative Control Probe"    = "NegControlProbe",
            "Blank Codeword"            = "UnassignedCodeword",
            "Genomic Control"           = "GenomicControl",
            "Unassigned Codeword"       = "UnassignedCodeword",
            "Deprecated Codeword"       = "DeprecatedCodeword",
            gsub(" ", "_", fc)
        )
    }
    uniq_classes <- unique(feat_classes_vec)
    if (length(uniq_classes) > 1L && isTRUE(split_by_type)) {
        store_list <- lapply(uniq_classes, function(fc) {
            idx <- which(feat_classes_vec == fc)
            pe[idx, , drop = FALSE]
        })
        names(store_list) <- vapply(uniq_classes, feat_class_to_name,
                                    character(1L))
    } else {
        # No split — single store, name from the (sole) class or default rna
        nm <- if (length(uniq_classes) == 1L) feat_class_to_name(uniq_classes)
              else "rna"
        store_list <- stats::setNames(list(pe), nm)
    }

    if (output == "store") return(store_list)

    # Build exprObjs directly (bypass createExprObj / .evaluate_expr_matrix
    # since they don't recognize parquetExprStore).
    lapply(seq_along(store_list), function(i) {
        methods::new("exprObj",
            name       = "raw",
            exprMat    = store_list[[i]],
            spat_unit  = "cell",
            feat_type  = names(store_list)[[i]],
            provenance = "cell"
        )
    })
}


# Feature class extractor: 10x mtx triple. Reads column 3 of features.tsv
# (the canonical feature-type column written by 10x); falls back to
# "Gene Expression" when only ID + name are present. Shared across 10x-
# format readers (Xenium / VisiumHD).
.tenx_feat_classes_mtx <- function(path) {
    files_10X <- list.files(path)
    features_file <- grep(files_10X, pattern = "features|genes",
                          value = TRUE)[1L]
    if (is.na(features_file)) {
        stop("[tenx_feat_classes_mtx] features.tsv not found under ", path,
             call. = FALSE)
    }
    featuresDT <- data.table::fread(
        input  = file.path(path, features_file),
        header = FALSE
    )
    if (ncol(featuresDT) >= 3L) as.character(featuresDT$V3)
    else rep("Gene Expression", nrow(featuresDT))
}


# Feature class extractor: 10x h5. Reads /<root>/features/feature_type.
# Opens its own h5 handle (separate from tenxH5Input's read); cheap since
# feature_type is a small vector. Shared across 10x-format readers
# (Xenium / VisiumHD).
.tenx_feat_classes_h5 <- function(path) {
    GiottoUtils::package_check("hdf5r")
    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)
    root <- names(h5)[1L]
    as.character(h5[[paste0(root, "/features/feature_type")]][])
}