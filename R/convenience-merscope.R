# MERSCOPE (Vizgen MERFISH) ingest pipeline for `gDirSource`-managed projects.


# CLASS ####


setClass(
    "MerscopeDiskReader",
    slots = list(
        merscope_dir    = "character",
        backend         = "ANY",
        score_threshold = "ANY",
        polygon_format  = "character",
        calls           = "list",
        paths           = "list"
    ),
    prototype = list(
        merscope_dir    = NA_character_,
        backend         = NULL,
        score_threshold = NULL,
        polygon_format  = "parquet",
        calls           = list(),
        paths           = list()
    )
)

# * show ####
setMethod("show", signature("MerscopeDiskReader"), function(object) {
    cat(sprintf("Giotto <%s>\n", class(object)[1L]))
    cat(sprintf("dir        : %s\n", object@merscope_dir))
    cat(sprintf("poly format: %s\n", object@polygon_format))
    cat(sprintf("score cutof: %s\n",
        if (is.null(object@score_threshold)) "none" else object@score_threshold))
    cat(sprintf("backend    : %s\n",
        if (is.null(object@backend)) "none" else class(object@backend)[1L]))
    cat(sprintf("funs       : %s\n", paste(names(object@calls), collapse = ", ")))
    invisible(object)
})

# * accessors ####
#' @export
setMethod("$", signature("MerscopeDiskReader"), function(x, name) {
    basic_info <- c("merscope_dir", "polygon_format", "score_threshold",
                    "backend", "paths")
    if (name %in% basic_info) {
        return(methods::slot(x, name))
    }
    x@calls[[name]]
})

#' @export
setMethod("$<-", signature("MerscopeDiskReader"), function(x, name, value) {
    basic_info <- c("merscope_dir", "polygon_format", "score_threshold",
                    "backend")
    if (name %in% basic_info) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }
    x@calls[[name]] <- value
    x
})

# * init ####
setMethod(
    "initialize", signature("MerscopeDiskReader"),
    function(.Object, ..., merscope_dir, backend, score_threshold,
             polygon_format) {
        .Object <- callNextMethod(.Object, ...)

        if (!missing(merscope_dir)) {
            .Object@merscope_dir <- normalizePath(path.expand(merscope_dir),
                                                  mustWork = TRUE)
        }
        if (!missing(score_threshold)) .Object@score_threshold <- score_threshold
        if (!missing(polygon_format)) {
            .Object@polygon_format <- match.arg(polygon_format,
                                                c("parquet", "hdf5"))
        }
        if (!missing(backend)) {
            if (is.character(backend)) {
                backend <- gDirSource(path = backend)
            }
            checkmate::assert_class(backend, "gsource")
            .Object@backend <- backend
        }
        if (is.null(.Object@backend)) {
            stop("[MerscopeDiskReader] `backend` is required", call. = FALSE)
        }
        if (is.na(.Object@merscope_dir)) {
            stop("[MerscopeDiskReader] `merscope_dir` is required", call. = FALSE)
        }

        .Object@paths <- .merscope_detect_paths(.Object@merscope_dir)

        list2env(.Object@paths, envir = environment())
        gsrc     <- .Object@backend
        score_th <- .Object@score_threshold
        poly_fmt <- .Object@polygon_format

        # transcripts
        tx_fun <- function(
            path = tx_path,
            feat_type = c("rna", "neg_probe"),
            split_keyword = list("Blank"),
            FOVs = NULL,
            dropcols = c("x", "y", "barcode_id", "transcript_id"),
            score_threshold = score_th,
            flip_vertical = TRUE,
            output = c("giottoPoints", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_transcript_disk(
                path = path, gsource = gsrc, feat_type = feat_type,
                split_keyword = split_keyword, FOVs = FOVs,
                dropcols = dropcols, score_threshold = score_threshold,
                flip_vertical = flip_vertical,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_transcripts <- tx_fun

        # polygons
        poly_fun <- function(
            path = cell_bound_path,
            name = "cell",
            poly_z_indices = NULL,
            polygon_format = poly_fmt,
            FOVs = NULL,
            metadata_path = cell_meta_path,
            flip_vertical = TRUE,
            output = c("giottoPolygon", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_poly_disk(
                path = path, gsource = gsrc, name = name,
                poly_z_indices = poly_z_indices,
                polygon_format = polygon_format,
                FOVs = FOVs, metadata_path = metadata_path,
                flip_vertical = flip_vertical,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_polys <- poly_fun

        # expression
        ex_fun <- function(
            path = expr_path,
            spat_unit = "cell",
            split_by_type = TRUE,
            feat_type = c("rna", "neg_probe"),
            split_keyword = "Blank",
            remove_zero_rows = TRUE,
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_expression_disk(
                path = path, gsource = gsrc, spat_unit = spat_unit,
                split_by_type = split_by_type,
                feat_type = feat_type, split_keyword = split_keyword,
                remove_zero_rows = remove_zero_rows,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_expression <- ex_fun

        # cell metadata (in memory)
        cm_fun <- function(path = cell_meta_path, spat_unit = "cell",
                           FOVs = NULL, verbose = NULL, ...) {
            .merscope_cellmeta(path = path, spat_unit = spat_unit, FOVs = FOVs,
                               verbose = verbose, ...)
        }
        .Object@calls$load_cellmeta <- cm_fun

        # spatial locations (cell centers from cell metadata, in memory)
        sl_fun <- function(path = cell_meta_path, spat_unit = "cell",
                           FOVs = NULL, flip_vertical = TRUE,
                           verbose = NULL, ...) {
            .merscope_spatlocs(path = path, spat_unit = spat_unit, FOVs = FOVs,
                               flip_vertical = flip_vertical,
                               verbose = verbose, ...)
        }
        .Object@calls$load_spatlocs <- sl_fun

        # images (in memory, file-backed)
        img_fun <- function(path, name = "image",
                            transform_path = transform_path,
                            negative_y = TRUE, verbose = NULL, ...) {
            .merscope_image(path = path, name = name,
                            transform_path = transform_path,
                            negative_y = negative_y, verbose = verbose)
        }
        .Object@calls$load_image <- img_fun

        # orchestrator
        gobject_fun <- function(
            transcript_path = tx_path,
            boundary_path = cell_bound_path,
            expression_path = expr_path,
            metadata_path = cell_meta_path,
            data_to_use = c("subcellular", "aggregate"),
            feat_type = c("rna", "neg_probe"),
            split_keyword = list("Blank"),
            poly_z_indices = NULL,
            FOVs = NULL,
            dropcols = c("x", "y", "barcode_id", "transcript_id"),
            spat_unit = "cell",
            load_transcripts = TRUE,
            load_expression = TRUE,
            load_cellmeta = FALSE,
            load_images = NULL,
            image_transform_path = transform_path,
            instructions = NULL,
            verbose = NULL) {

            data_to_use <- match.arg(data_to_use, c("subcellular", "aggregate"))
            load_transcripts <- as.logical(load_transcripts)
            load_expression  <- as.logical(load_expression)
            load_cellmeta    <- as.logical(load_cellmeta)

            if (data_to_use == "aggregate") {
                load_transcripts <- FALSE
                load_expression  <- TRUE
            }

            if (!load_transcripts && !load_expression) {
                warning(GiottoUtils::wrap_txt(
                    "One of either transcripts or expression info should be loaded for a fully functioning object"
                ))
            }

            funs <- .Object@calls

            tx_list <- NULL
            if (load_transcripts) {
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    FOVs = FOVs,
                    dropcols = dropcols,
                    verbose = verbose
                )
            }

            b <- NULL
            if (data_to_use == "subcellular" && !is.na(boundary_path)) {
                b <- funs$load_polys(
                    path = boundary_path,
                    name = spat_unit,
                    poly_z_indices = poly_z_indices,
                    FOVs = FOVs,
                    metadata_path = metadata_path,
                    verbose = verbose
                )
                if (!is.null(b) && !inherits(b, "list")) b <- list(b)
            }

            ex <- NULL
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    spat_unit = spat_unit,
                    split_by_type = TRUE,
                    feat_type = feat_type,
                    split_keyword = if (is.list(split_keyword)) {
                        unlist(split_keyword)
                    } else split_keyword,
                    verbose = verbose
                )
            }

            cx <- NULL
            if (load_cellmeta) {
                cx <- funs$load_cellmeta(path = metadata_path,
                                         spat_unit = spat_unit, FOVs = FOVs,
                                         verbose = verbose)
            }

            sl <- NULL
            if (is.null(b) && !is.na(metadata_path)) {
                sl <- funs$load_spatlocs(path = metadata_path,
                                         spat_unit = spat_unit, FOVs = FOVs,
                                         verbose = verbose)
            }

            imglist <- NULL
            if (!is.null(load_images)) {
                img_paths <- .merscope_image_paths(load_images, image_dir)
                imglist <- list()
                for (img_i in seq_along(img_paths)) {
                    im <- funs$load_image(
                        path = img_paths[[img_i]],
                        name = names(img_paths)[[img_i]],
                        transform_path = image_transform_path,
                        verbose = verbose
                    )
                    imglist <- c(imglist, im)
                }
            }

            # init gobject with disk backend, then attach
            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            if (!is.null(tx_list)) {
                g <- GiottoClass::setGiotto(g, tx_list, verbose = FALSE)
            }
            if (!is.null(b)) {
                g <- GiottoClass::setGiotto(g, b,
                    centroids_to_spatlocs = TRUE, verbose = FALSE)
            }
            if (!is.null(ex)) g <- GiottoClass::setGiotto(g, ex)
            if (!is.null(cx)) g <- GiottoClass::setGiotto(g, cx, verbose = FALSE)
            if (!is.null(sl)) g <- GiottoClass::setGiotto(g, sl, verbose = FALSE)
            if (!is.null(imglist)) {
                g <- GiottoClass::setGiotto(g, imglist, verbose = FALSE)
            }

            GiottoUtils::vmsg(.v = verbose, "done")
            g
        }
        .Object@calls$create_gobject <- gobject_fun

        .Object
    }
)


# CREATE READER ####

#' @title Import a Vizgen MERSCOPE assay (disk-backed)
#' @name importMerscopeDisk
#' @description
#' Disk-backed MERSCOPE reader. Produces a `MerscopeDiskReader` whose
#' `load_transcripts()`, `load_polys()` and `load_expression()` calls write to a
#' `gDirSource`-managed project vault. Cell metadata remains in memory.
#' @param merscope_dir MERSCOPE region directory (the one containing
#'   `detected_transcripts.parquet`). `experiment.json` is read from the parent
#'   directory.
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#' @param score_threshold `numeric` or `NULL` (default). Minimum
#'   `transcript_score` retained; `NULL` keeps all transcripts.
#' @param polygon_format `character`. `"parquet"` (default) or `"hdf5"`.
#' @returns `MerscopeDiskReader` object
#' @seealso [importXeniumDisk()]
#' @export
importMerscopeDisk <- function(merscope_dir = NULL,
                               backend,
                               score_threshold = NULL,
                               polygon_format = c("parquet", "hdf5")) {
    if (missing(backend)) {
        stop("[importMerscopeDisk] `backend` is required", call. = FALSE)
    }
    polygon_format <- match.arg(polygon_format, c("parquet", "hdf5"))
    a <- list(Class = "MerscopeDiskReader", backend = backend,
              score_threshold = score_threshold,
              polygon_format = polygon_format)
    if (!is.null(merscope_dir)) a$merscope_dir <- merscope_dir
    do.call(new, args = a)
}


#' @title Create a disk-backed MERSCOPE Giotto object
#' @name createGiottoMerscopeObjectDisk
#' @description
#' Create a disk-backed `giotto` object from a Vizgen MERSCOPE region
#' directory. The y coordinates of transcripts, polygons and cell centers are
#' negated (`y` becomes `-y`).
#' @param merscope_dir MERSCOPE region directory
#' @param backend a `gsource` project backend
#' @param score_threshold minimum `transcript_score`, or `NULL` (default)
#' @param polygon_format `"parquet"` (default) or `"hdf5"`
#' @param data_to_use `"subcellular"` (default) builds from transcripts and
#'   polygons; `"aggregate"` uses only the vendor cell x gene matrix.
#' @param poly_z_indices which z-indices to use for polygons. `NULL` (default)
#'   uses the lowest plane when all planes hold identical boundaries and stops
#'   otherwise.
#' @param FOVs which FOVs to load. `NULL` (default) loads all.
#' @param dropcols transcript columns not written to the store. Default
#'   `c("x", "y", "barcode_id", "transcript_id")`; `character(0)` keeps all.
#' @param feat_type feature type names: the first for detections matching no
#'   keyword, then one per `split_keyword`. Default `c("rna", "neg_probe")`.
#' @param split_keyword keywords separating feature detections. Default
#'   `list("Blank")`, the MERSCOPE negative-control prefix.
#' @param spat_unit name for the cell spatial unit (default `"cell"`)
#' @param load_transcripts,load_expression,load_cellmeta logical
#' @param load_images images to load: mosaic names (e.g. `"DAPI_z3"` for
#'   `images/mosaic_DAPI_z3.tif`) or image file paths. List names are used as
#'   image names. Images are placed with
#'   `images/micron_to_mosaic_pixel_transform.csv`. `NULL` (default) loads no
#'   images.
#' @param instructions Giotto instructions
#' @param verbose verbosity
#' @returns a `giotto` object backed by `backend`
#' @export
createGiottoMerscopeObjectDisk <- function(
    merscope_dir,
    backend,
    score_threshold = NULL,
    polygon_format = c("parquet", "hdf5"),
    data_to_use = c("subcellular", "aggregate"),
    poly_z_indices = NULL,
    FOVs = NULL,
    dropcols = c("x", "y", "barcode_id", "transcript_id"),
    feat_type = c("rna", "neg_probe"),
    split_keyword = list("Blank"),
    spat_unit = "cell",
    load_transcripts = TRUE,
    load_expression = TRUE,
    load_cellmeta = FALSE,
    load_images = NULL,
    instructions = NULL,
    verbose = NULL
) {
    x <- importMerscopeDisk(
        merscope_dir = merscope_dir,
        backend = backend,
        score_threshold = score_threshold,
        polygon_format = polygon_format
    )
    x$create_gobject(
        data_to_use = data_to_use,
        feat_type = feat_type,
        split_keyword = split_keyword,
        poly_z_indices = poly_z_indices,
        FOVs = FOVs,
        dropcols = dropcols,
        spat_unit = spat_unit,
        load_transcripts = load_transcripts,
        load_expression = load_expression,
        load_cellmeta = load_cellmeta,
        load_images = load_images,
        instructions = instructions,
        verbose = verbose
    )
}


# PATHS ####

.merscope_detect_paths <- function(merscope_dir) {
    f <- function(...) {
        p <- file.path(merscope_dir, ...)
        if (file.exists(p)) normalizePath(p) else NULL
    }
    parent <- dirname(merscope_dir)
    fp <- function(...) {
        p <- file.path(parent, ...)
        if (file.exists(p)) normalizePath(p) else NULL
    }

    # transcripts: prefer parquet over csv
    tx_path <- f("detected_transcripts.parquet") %||% f("detected_transcripts.csv")

    # boundaries: parquet file, else the per-FOV hdf5 directory
    cell_bound_path <- f("cell_boundaries.parquet") %||% f("cell_boundaries")

    paths <- list(
        tx_path          = tx_path,
        cell_bound_path  = cell_bound_path,
        expr_path        = f("cell_by_gene.csv"),
        cell_meta_path   = f("cell_metadata.csv"),
        cell_cat_path    = f("cell_categories.csv"),
        cell_numcat_path = f("cell_numeric_categories.csv"),
        experiment_path  = fp("experiment.json"),
        image_dir        = f("images"),
        transform_path   = f("images", "micron_to_mosaic_pixel_transform.csv")
    )

    if (is.null(paths$tx_path) && is.null(paths$expr_path)) {
        stop("[merscope] no `detected_transcripts.*` or `cell_by_gene.csv` ",
             "found in:\n  ", merscope_dir,
             "\nIs this a MERSCOPE region directory?", call. = FALSE)
    }

    lapply(paths, function(p) p %||% NA_character_)
}


# MODULAR ####


## transcripts ####

# Name unnamed columns in a MERSCOPE transcript schema.
.merscope_tx_schema <- function(path) {
    sc <- arrow::schema(arrow::open_dataset(path))
    nms <- names(sc)
    unnamed <- !nzchar(nms)
    if (any(unnamed)) {
        nms[unnamed] <- paste0(".merscope_unnamed_", seq_len(sum(unnamed)))
    }
    list(
        schema = do.call(arrow::schema,
            lapply(seq_along(nms), function(i) arrow::field(nms[i], sc[[i]]$type))),
        dropcols = nms[unnamed]
    )
}

# Lazy transcript query with the schema patched and index columns removed.
.merscope_tx_query <- function(path, score_threshold = NULL, FOVs = NULL,
                               dropcols = character(0L),
                               flip_vertical = FALSE) {
    sch <- .merscope_tx_schema(path)
    a <- arrow::open_dataset(sources = path, schema = sch$schema)
    drop <- unique(c(sch$dropcols, as.character(dropcols)))
    if (length(drop) > 0L) {
        a <- dplyr::select(a, -dplyr::any_of(drop))
    }
    if (!is.null(score_threshold)) {
        a <- dplyr::filter(a, transcript_score >= score_threshold)
    }
    if (!is.null(FOVs)) {
        fovs_ <- as.integer(FOVs)
        a <- dplyr::filter(a, fov %in% fovs_)
    }
    if ("cell_id" %in% names(a)) {
        a <- dplyr::mutate(a, cell_id = arrow::cast(cell_id, arrow::string()))
    }
    if (isTRUE(flip_vertical)) {
        a <- dplyr::mutate(a, global_y = -global_y)
    }
    a
}

# Disk-backed MERSCOPE transcript ingestion.
.merscope_transcript_disk <- function(
    path,
    gsource,
    feat_type = c("rna", "neg_probe"),
    split_keyword = list("Blank"),
    FOVs = NULL,
    dropcols,
    score_threshold = NULL,
    flip_vertical = TRUE,
    output = c("giottoPoints", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path) || is.na(path)) {
        stop("[merscope_transcript_disk] no transcripts path provided",
             call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("giottoPoints", "store"))
    GiottoUtils::package_check("arrow")
    GiottoUtils::package_check("dplyr")

    GiottoUtils::vmsg("[merscope_transcript_disk] loading transcripts",
                      .v = verbose)

    score_thr <- score_threshold
    fovs_     <- if (is.null(FOVs)) NULL else as.integer(FOVs)
    drop      <- as.character(dropcols)

    flip <- isTRUE(flip_vertical)
    read_fun <- function(x, schema = NULL) {
        .merscope_tx_query(x, score_threshold = score_thr, FOVs = fovs_,
                           dropcols = drop, flip_vertical = flip)
    }

    fs <- fileStore(path = path, read_fun = read_fun)
    qs <- as(fs, "queryableStore")

    tx_store <- sourceWrite(
        gsource, qs,
        store_type = "parquetGeomTile",
        type = "points",
        id_col = "gene",
        sdimx = "global_x",
        sdimy = "global_y",
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


## polygons ####

# Determine the z-plane architecture of a MERSCOPE boundary table.
.merscope_zplane_architecture <- function(path, poly_z_indices = NULL,
                                          n_sample = 200L, verbose = NULL) {
    ds <- arrow::open_dataset(path)
    if (!"ZIndex" %in% names(ds)) {
        return(list(kind = "2d", z_indices = NA_integer_, use_z = NA_integer_))
    }

    z_all <- sort(ds |> dplyr::distinct(ZIndex) |> dplyr::collect() |>
                  dplyr::pull(ZIndex))
    if (!is.null(poly_z_indices)) {
        keep <- intersect(as.integer(poly_z_indices), as.integer(z_all))
        if (!length(keep)) {
            stop("[merscope] none of `poly_z_indices` (",
                 paste(poly_z_indices, collapse = ", "),
                 ") are present. Available: ", paste(z_all, collapse = ", "),
                 call. = FALSE)
        }
        z_all <- keep
    }

    if (length(z_all) == 1L) {
        GiottoUtils::vmsg(.v = verbose,
            "[merscope] single z-plane; 2D architecture")
        return(list(kind = "2d", z_indices = z_all, use_z = z_all[[1L]]))
    }

    zA <- z_all[[1L]]; zB <- z_all[[2L]]
    ids_A <- ds |> dplyr::filter(ZIndex == zA) |> dplyr::select(EntityID) |>
             dplyr::collect() |> dplyr::pull(EntityID)
    ids_B <- ds |> dplyr::filter(ZIndex == zB) |> dplyr::select(EntityID) |>
             dplyr::collect() |> dplyr::pull(EntityID)

    if (!setequal(ids_A, ids_B)) {
        return(.merscope_stop_3d(z_all, reason = "cell ID sets differ between z-planes"))
    }

    # sample cells and compare serialized geometry across all planes
    set.seed(1L)
    samp <- if (length(ids_A) > n_sample) sample(ids_A, n_sample) else ids_A
    g <- ds |> dplyr::filter(EntityID %in% samp) |>
         dplyr::select(EntityID, ZIndex, Geometry) |> dplyr::collect()
    sp <- split(g, g$EntityID)
    same <- vapply(sp, function(gi) {
        gi <- gi[order(gi$ZIndex), ]
        ref <- as.raw(gi$Geometry[[1L]])
        all(vapply(seq_len(nrow(gi)),
                   function(i) identical(as.raw(gi$Geometry[[i]]), ref),
                   logical(1L)))
    }, logical(1L))

    if (!all(same)) {
        return(.merscope_stop_3d(z_all, reason = sprintf(
            "geometry differs across z-planes for %d of %d sampled cells",
            sum(!same), length(same))))
    }

    GiottoUtils::vmsg(.v = verbose, sprintf(
        "[merscope] [Auto-Detection] Identical IDs AND vertices across %d z-planes (Replicated 2D Data).\n[merscope] [Auto-Correction] Using z-plane %s, dropping %d redundant plane(s).",
        length(z_all), zA, length(z_all) - 1L))
    list(kind = "replicated_2d", z_indices = z_all, use_z = zA)
}

.merscope_stop_3d <- function(z_indices, reason) {
    stop(GiottoUtils::wrap_txt(sprintf(
        "[merscope] true 3D segmentation detected (%s).
Only 2D and replicated-2D boundaries are supported.
Available z-indices: %s
To load a single plane, pass e.g. poly_z_indices = %s",
        reason, paste(z_indices, collapse = ", "), z_indices[[1L]]
    )), call. = FALSE)
}

# Disk-backed MERSCOPE polygon ingestion.
.merscope_poly_disk <- function(
    path,
    gsource,
    name = "cell",
    poly_z_indices = NULL,
    polygon_format = c("parquet", "hdf5"),
    FOVs = NULL,
    metadata_path = NULL,
    threshold = 25000L,
    flip_vertical = TRUE,
    output = c("giottoPolygon", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path) || is.na(path)) {
        stop("[merscope_poly_disk] no boundary path provided", call. = FALSE)
    }
    checkmate::assert_class(gsource, "gsource")
    checkmate::assert_character(name, len = 1L)
    output <- match.arg(output, choices = c("giottoPolygon", "store"))
    polygon_format <- match.arg(polygon_format, c("parquet", "hdf5"))
    GiottoUtils::package_check("arrow")
    GiottoUtils::package_check("terra")

    if (polygon_format == "hdf5") {
        stop(GiottoUtils::wrap_txt(
            "[merscope] `polygon_format = 'hdf5'` is not implemented.
Use cell_boundaries.parquet."
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)

    GiottoUtils::vmsg(
        sprintf("[merscope_poly_disk] loading boundary '%s'", name),
        .v = verbose
    )

    arch <- .merscope_zplane_architecture(path,
        poly_z_indices = poly_z_indices, verbose = verbose)

    ds <- arrow::open_dataset(path)
    q <- ds
    if (!is.na(arch$use_z) && "ZIndex" %in% names(ds)) {
        use_z <- arch$use_z
        q <- dplyr::filter(q, ZIndex == use_z)
    }

    if (!is.null(FOVs)) {
        if (is.null(metadata_path) || is.na(metadata_path)) {
            stop("[merscope_poly_disk] `FOVs` requires cell_metadata.csv to ",
                 "map EntityID -> fov, but no metadata path was found.",
                 call. = FALSE)
        }
        cm <- data.table::fread(metadata_path, select = c("EntityID", "fov"))
        keep_ids <- cm[["EntityID"]][cm[["fov"]] %in% as.integer(FOVs)]
        if (!length(keep_ids)) {
            stop("[merscope_poly_disk] no cells found in FOV(s): ",
                 paste(FOVs, collapse = ", "), call. = FALSE)
        }
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "[merscope_poly_disk] FOV subset -> %s cells",
            format(length(keep_ids), big.mark = ",")))
        q <- dplyr::filter(q, EntityID %in% keep_ids)
    }

    x <- q |> dplyr::select(EntityID, Geometry) |> dplyr::collect()

    if (nrow(x) == 0L) {
        GiottoUtils::vmsg(.v = verbose, "[merscope_poly_disk] no polygons found")
        return(NULL)
    }
    GiottoUtils::vmsg(.v = verbose,
        sprintf("[merscope_poly_disk] decoding %s geometries",
                format(nrow(x), big.mark = ",")))

    sv <- terra::vect(as.list(x$Geometry))
    terra::values(sv) <- data.frame(
        poly_ID = as.character(x$EntityID),
        stringsAsFactors = FALSE
    )
    if (isTRUE(flip_vertical)) {
        gm <- terra::geom(sv)
        gm[, "y"] <- -gm[, "y"]
        vals <- terra::values(sv)
        sv <- terra::vect(gm, type = "polygons", crs = terra::crs(sv))
        terra::values(sv) <- vals
    }

    pg <- sourceWrite(gsource, sv, store_type = "parquetGeom")

    poly_store <- sourceWrite(gsource, pg,
        store_type = "parquetGeomTile", threshold = threshold,
        verbose = verbose, ...)

    if (output == "store") return(poly_store)

    GiottoClass::createGiottoPolygon(x = poly_store, name = name)
}


## expression ####

.merscope_feat_classes <- function(feat_ids, split_keyword = "Blank",
                                   feat_type = c("rna", "neg_probe")) {
    kw <- as.character(unlist(split_keyword))
    if (length(feat_type) != length(kw) + 1L) {
        stop("[merscope] `feat_type` needs one name more than `split_keyword`",
             call. = FALSE)
    }
    out <- rep(feat_type[[1L]], length(feat_ids))
    for (i in seq_along(kw)) {
        hit <- grepl(kw[[i]], feat_ids, fixed = TRUE)
        out[hit] <- feat_type[[i + 1L]]
    }
    out
}

# Disk-backed MERSCOPE expression ingestion.
.merscope_expression_disk <- function(
    path,
    gsource,
    spat_unit = "cell",
    split_by_type = TRUE,
    feat_type = c("rna", "neg_probe"),
    split_keyword = "Blank",
    remove_zero_rows = TRUE,
    output = c("exprObj", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path) || is.na(path)) {
        stop("[merscope_expression_disk] no expression path provided",
             call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("exprObj", "store"))

    GiottoUtils::vmsg("[merscope_expression_disk] loading cell x gene matrix",
                      .v = verbose)

    # `cell` is the MERSCOPE cell id column in cell_by_gene.csv
    inp <- csvWideInput(path, cell_id_col = "cell")
    feat_classes_vec <- .merscope_feat_classes(inp@feat_ids,
                                               split_keyword = split_keyword,
                                               feat_type = feat_type)

    pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                      verbose = verbose, ...)

    if (isTRUE(remove_zero_rows)) {
        col_id <- NULL # NSE binding
        ds <- storeRead(pe)
        present <- ds |> dplyr::distinct(col_id) |> dplyr::collect()
        keep_idx <- sort(as.integer(present$col_id))
        if (length(keep_idx) > 0L &&
            length(keep_idx) < length(feat_classes_vec)) {
            pe <- pe[keep_idx, , drop = FALSE]
            feat_classes_vec <- feat_classes_vec[keep_idx]
        }
    }

    uniq_classes <- unique(feat_classes_vec)
    if (length(uniq_classes) > 1L && isTRUE(split_by_type)) {
        store_list <- lapply(uniq_classes, function(fc) {
            idx <- which(feat_classes_vec == fc)
            pe[idx, , drop = FALSE]
        })
        names(store_list) <- uniq_classes
    } else {
        nm <- if (length(uniq_classes) == 1L) uniq_classes else "rna"
        store_list <- stats::setNames(list(pe), nm)
    }

    if (output == "store") return(store_list)

    lapply(seq_along(store_list), function(i) {
        methods::new("exprObj",
            name       = "raw",
            exprMat    = store_list[[i]],
            spat_unit  = spat_unit,
            feat_type  = names(store_list)[[i]],
            provenance = spat_unit
        )
    })
}


## cell metadata ####

.merscope_cellmeta <- function(path, spat_unit = "cell", FOVs = NULL,
                               verbose = NULL, ...) {
    if (missing(path) || is.na(path)) {
        stop("[merscope_cellmeta] no cell metadata path provided", call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    GiottoUtils::vmsg("[merscope_cellmeta] loading cell metadata", .v = verbose)

    dt <- data.table::fread(path)

    if (!is.null(FOVs) && "fov" %in% names(dt)) {
        n_before <- nrow(dt)
        dt <- dt[dt[["fov"]] %in% as.integer(FOVs)]
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "[merscope_cellmeta] FOV subset -> %s of %s cells",
            format(nrow(dt), big.mark = ","), format(n_before, big.mark = ",")))
    }

    if ("EntityID" %in% names(dt)) {
        data.table::setnames(dt, "EntityID", "cell_ID")
        dt[, cell_ID := as.character(cell_ID)]
    }
    GiottoClass::createCellMetaObj(
        metadata = dt,
        spat_unit = spat_unit,
        feat_type = "rna",
        provenance = spat_unit
    )
}


## spatial locations ####

.merscope_spatlocs <- function(path, spat_unit = "cell", FOVs = NULL,
                               flip_vertical = TRUE,
                               verbose = NULL, ...) {
    if (missing(path) || is.na(path)) {
        stop("[merscope_spatlocs] no cell metadata path provided", call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    GiottoUtils::vmsg("[merscope_spatlocs] loading cell centers", .v = verbose)
    sdimy <- NULL # NSE var

    dt <- data.table::fread(path)
    if (!is.null(FOVs) && "fov" %in% names(dt)) {
        dt <- dt[dt[["fov"]] %in% as.integer(FOVs)]
    }
    sl <- data.table::data.table(
        sdimx = dt[["center_x"]],
        sdimy = dt[["center_y"]],
        cell_ID = as.character(dt[["EntityID"]])
    )
    if (isTRUE(flip_vertical)) {
        sl[, sdimy := -sdimy]
    }
    GiottoClass::createSpatLocsObj(
        coordinates = sl,
        name = "raw",
        spat_unit = spat_unit,
        provenance = spat_unit,
        verbose = FALSE
    )
}


## images ####

.merscope_image_paths <- function(load_images, image_dir) {
    load_images <- as.list(load_images)
    nms <- names(load_images)
    if (is.null(nms)) nms <- rep("", length(load_images))
    paths <- character(length(load_images))
    for (i in seq_along(load_images)) {
        v <- as.character(load_images[[i]])
        if (file.exists(v)) {
            paths[[i]] <- normalizePath(v)
        } else {
            if (is.na(image_dir)) {
                stop("[merscope] no `images` directory found for image ", v,
                     call. = FALSE)
            }
            paths[[i]] <- file.path(image_dir, paste0("mosaic_", v, ".tif"))
            if (!file.exists(paths[[i]])) {
                stop("[merscope] image not found: ", paths[[i]], call. = FALSE)
            }
        }
        if (!nzchar(nms[[i]])) {
            nms[[i]] <- sub("^mosaic_", "", sub("(?i)\\.tiff?$", "",
                basename(paths[[i]]), perl = TRUE))
        }
    }
    stats::setNames(as.list(paths), nms)
}

.merscope_image <- function(path, name = "image", transform_path = NULL,
                            negative_y = TRUE, verbose = NULL) {
    if (missing(path)) {
        stop("[merscope_image] no image path provided", call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    if (is.null(transform_path) || is.na(transform_path)) {
        stop("[merscope_image] no micron_to_mosaic_pixel_transform.csv found ",
             "for image placement", call. = FALSE)
    }
    checkmate::assert_file_exists(transform_path)
    GiottoUtils::vmsg(.v = verbose,
        sprintf("[merscope_image] loading image as '%s'", name))

    tm <- as.matrix(data.table::fread(transform_path, header = FALSE))
    sx <- tm[1L, 1L]
    sy <- tm[2L, 2L]
    tx <- tm[1L, 3L]
    ty <- tm[2L, 3L]

    img <- GiottoClass::createGiottoLargeImage(
        path, name = name, negative_y = negative_y, verbose = FALSE
    )
    img <- GiottoClass::rescale(img, fx = 1 / sx, fy = 1 / sy, x0 = 0, y0 = 0)
    dy <- if (isTRUE(negative_y)) ty / sy else -ty / sy
    GiottoClass::spatShift(img, dx = -tx / sx, dy = dy)
}
