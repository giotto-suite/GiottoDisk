# MERSCOPE (Vizgen MERFISH) ingest pipeline for `gDirSource`-managed projects.
#
# Laid out like the Xenium disk reader: class and init (loaders, then
# create_gobject), the import function, then the modular loaders. Transcripts,
# polygons and expression are written to the project as stores; cell metadata,
# cell centers and images stay in memory. `create_gobject()` takes the same
# arguments as Giotto's in-memory MERSCOPE reader, so
# Giotto::createGiottoMerscopeObject(backend =) runs it unchanged.
#
# MERSCOPE specifics kept: `FOVs` (applied to every data type, the expression
# matrix included), replicated z-planes in the boundaries, and the two boundary
# formats (`cell_boundaries.parquet`, or per-FOV HDF5 files in
# `cell_boundaries/`). Transcripts may be parquet or csv; a csv is converted to
# parquet once and cached, since tiling scans its source once per tile.


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
        bound_path <- if (poly_fmt == "hdf5") cell_bound_dir else cell_bound_path
        # not `transform_path = transform_path` below: a default that names its
        # own argument cannot be evaluated
        tf_path <- transform_path

        ## transcripts ####
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

        ## polygons ####
        poly_fun <- function(
            path = bound_path,
            name = "cell",
            poly_z_indices = NULL,
            polygon_format = poly_fmt,
            FOVs = NULL,
            metadata_path = cell_meta_path,
            cells = NULL,
            flip_vertical = TRUE,
            output = c("giottoPolygon", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_poly_disk(
                path = path, gsource = gsrc, name = name,
                poly_z_indices = poly_z_indices,
                polygon_format = polygon_format,
                FOVs = FOVs, metadata_path = metadata_path, cells = cells,
                flip_vertical = flip_vertical,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_polys <- poly_fun

        ## cell metadata (in memory) ####
        cm_fun <- function(path = cell_meta_path, spat_unit = "cell",
                           FOVs = NULL, cells = NULL, verbose = NULL, ...) {
            .merscope_cellmeta(path = path, spat_unit = spat_unit, FOVs = FOVs,
                               cells = cells, verbose = verbose)
        }
        .Object@calls$load_cellmeta <- cm_fun

        ## spatial locations (cell centers, in memory) ####
        sl_fun <- function(path = cell_meta_path, spat_unit = "cell",
                           FOVs = NULL, cells = NULL, flip_vertical = TRUE,
                           verbose = NULL, ...) {
            .merscope_spatlocs(path = path, spat_unit = spat_unit, FOVs = FOVs,
                               cells = cells, flip_vertical = flip_vertical,
                               verbose = verbose)
        }
        .Object@calls$load_spatlocs <- sl_fun

        ## expression ####
        ex_fun <- function(
            path = expr_path,
            spat_unit = "cell",
            split_by_type = TRUE,
            feat_type = c("rna", "neg_probe"),
            split_keyword = "Blank",
            FOVs = NULL,
            metadata_path = cell_meta_path,
            cells = NULL,
            remove_zero_rows = TRUE,
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_expression_disk(
                path = path, gsource = gsrc, spat_unit = spat_unit,
                split_by_type = split_by_type,
                feat_type = feat_type, split_keyword = split_keyword,
                FOVs = FOVs, metadata_path = metadata_path, cells = cells,
                remove_zero_rows = remove_zero_rows,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_expression <- ex_fun

        ## images (in memory, file-backed) ####
        img_fun <- function(path, name = "image", transform_path = tf_path,
                            negative_y = TRUE, verbose = NULL, ...) {
            .merscope_image(path = path, name = name,
                            transform_path = transform_path,
                            negative_y = negative_y, verbose = verbose)
        }
        .Object@calls$load_image <- img_fun

        ## create_gobject ####
        # Same arguments as Giotto's in-memory MerscopeReader$create_gobject().
        gobject_fun <- function(
            transcript_path = tx_path,
            boundary_path = bound_path,
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
            image_transform_path = tf_path,
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

            # cell_metadata.csv maps cells to FOVs and holds the cell centers:
            # read at most once, and only when a step needs it
            cell_env <- new.env()
            get_cells <- function() {
                if (is.null(cell_env$cells)) {
                    cell_env$cells <- .merscope_read_cells(metadata_path,
                        FOVs = FOVs)
                }
                cell_env$cells
            }
            fov_cells <- function() if (!is.null(FOVs)) get_cells()

            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            # transcripts
            if (load_transcripts) {
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    FOVs = FOVs,
                    dropcols = dropcols,
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, tx_list, verbose = FALSE)
            }

            # polygons; their centroids are attached as spatial locations
            b <- NULL
            if (data_to_use == "subcellular" && !is.na(boundary_path)) {
                b <- funs$load_polys(
                    path = boundary_path,
                    name = spat_unit,
                    poly_z_indices = poly_z_indices,
                    FOVs = FOVs,
                    metadata_path = metadata_path,
                    cells = fov_cells(),
                    verbose = verbose
                )
                if (!is.null(b)) {
                    g <- GiottoClass::setGiotto(g, b,
                        centroids_to_spatlocs = TRUE, verbose = FALSE)
                }
            }

            # expression
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    spat_unit = spat_unit,
                    split_by_type = TRUE,
                    feat_type = feat_type,
                    split_keyword = unlist(split_keyword),
                    FOVs = FOVs,
                    metadata_path = metadata_path,
                    cells = fov_cells(),
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, ex)
            }

            # cell metadata
            if (load_cellmeta) {
                cx <- funs$load_cellmeta(path = metadata_path,
                    spat_unit = spat_unit, cells = get_cells(),
                    verbose = verbose)
                g <- GiottoClass::setGiotto(g, cx, verbose = FALSE)
            }

            # images
            if (!is.null(load_images)) {
                img_paths <- .merscope_image_paths(load_images, image_dir)
                imglist <- lapply(seq_along(img_paths), function(img_i) {
                    funs$load_image(
                        path = img_paths[[img_i]],
                        name = names(img_paths)[[img_i]],
                        transform_path = image_transform_path,
                        verbose = verbose
                    )
                })
                g <- GiottoClass::setGiotto(g, imglist, verbose = FALSE)
            }

            # no polygons loaded: cell centers from the cell metadata
            if (is.null(b) && !is.na(metadata_path)) {
                sl <- funs$load_spatlocs(path = metadata_path,
                    spat_unit = spat_unit, cells = get_cells(),
                    verbose = verbose)
                g <- GiottoClass::setGiotto(g, sl, verbose = FALSE)
            }

            GiottoUtils::vmsg(.v = verbose, "done")
            g
        }
        .Object@calls$create_gobject <- gobject_fun

        .Object
    }
)

# * access ####
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


# CREATE READER ####

#' @title Import a Vizgen MERSCOPE assay (disk-backed)
#' @name importMerscopeDisk
#' @description
#' Disk-backed MERSCOPE reader. Produces a `MerscopeDiskReader` whose
#' `load_transcripts()`, `load_polys()` and `load_expression()` calls write to a
#' `gDirSource`-managed project vault; cell metadata, cell centers and images
#' stay in memory. `$create_gobject()` takes the same arguments as
#' [Giotto::importMerscope()]'s, and [Giotto::createGiottoMerscopeObject()]
#' with `backend =` uses this reader.
#'
#' Transcripts may be `detected_transcripts.parquet` or `.csv`; a csv is
#' converted to parquet once, cached by path, size and modification time.
#' Cell boundaries may be `cell_boundaries.parquet` or (with
#' `polygon_format = "hdf5"`) the per-FOV HDF5 files in `cell_boundaries/`.
#' `FOVs` applies to transcripts, polygons, expression, cell metadata and cell
#' centers. The y coordinates are negated (`y` becomes `-y`).
#' @param merscope_dir MERSCOPE region directory
#' @param backend a `gsource` (typically `gDirSource`) project backend, or a
#'   directory path for one
#' @param score_threshold `numeric` or `NULL` (default). Minimum
#'   `transcript_score` retained; `NULL` keeps all transcripts. Exports without
#'   a `transcript_score` column are not filtered.
#' @param polygon_format `character`. `"parquet"` (default) or `"hdf5"`.
#' @returns `MerscopeDiskReader` object
#' @seealso [importXeniumDisk()]
#' @examples
#' \dontrun{
#' x <- importMerscopeDisk("path/to/region_0", backend = "path/to/project")
#' g <- x$create_gobject()
#' }
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


# MODULAR ####

## paths ####

.merscope_detect_paths <- function(merscope_dir) {
    f <- function(...) {
        p <- file.path(merscope_dir, ...)
        if (file.exists(p)) normalizePath(p) else NULL
    }
    bound_dir <- file.path(merscope_dir, "cell_boundaries")

    paths <- list(
        # transcripts: parquet over csv
        tx_path         = f("detected_transcripts.parquet") %||%
                          f("detected_transcripts.csv"),
        cell_bound_path = f("cell_boundaries.parquet"),
        cell_bound_dir  = if (dir.exists(bound_dir)) normalizePath(bound_dir),
        expr_path       = f("cell_by_gene.csv"),
        cell_meta_path  = f("cell_metadata.csv"),
        image_dir       = f("images"),
        transform_path  = f("images", "micron_to_mosaic_pixel_transform.csv")
    )

    if (is.null(paths$tx_path) && is.null(paths$expr_path)) {
        stop("[merscope] no `detected_transcripts.*` or `cell_by_gene.csv` ",
             "found in:\n  ", merscope_dir,
             "\nIs this a MERSCOPE region directory?", call. = FALSE)
    }

    lapply(paths, function(p) p %||% NA_character_)
}


## transcript ####

# Schema of a parquet transcript file with unnamed (pandas index) columns
# named, and `cell_id` read as string: MERSCOPE cell IDs are 64-bit.
.merscope_tx_schema <- function(path, format = "parquet") {
    sc <- arrow::open_dataset(path, format = format)$schema
    nms <- names(sc)
    unnamed <- !nzchar(nms)
    if (any(unnamed)) {
        nms[unnamed] <- paste0(".merscope_unnamed_", seq_len(sum(unnamed)))
    }
    list(
        schema = do.call(arrow::schema, lapply(seq_along(nms), function(i) {
            arrow::field(nms[[i]],
                if (nms[[i]] == "cell_id") arrow::string() else sc[[i]]$type)
        })),
        dropcols = nms[unnamed]
    )
}

# Lazy transcript query with the schema patched and index columns removed.
.merscope_tx_query <- function(path, score_threshold = NULL, FOVs = NULL,
                               dropcols = character(0L),
                               flip_vertical = FALSE) {
    transcript_score <- fov <- global_y <- NULL # NSE vars
    sch <- .merscope_tx_schema(path)
    a <- arrow::open_dataset(sources = path, schema = sch$schema)
    drop <- unique(c(sch$dropcols, as.character(dropcols)))
    if (length(drop) > 0L) {
        a <- dplyr::select(a, -dplyr::any_of(drop))
    }
    if (!is.null(score_threshold) && "transcript_score" %in% names(a)) {
        a <- dplyr::filter(a, transcript_score >= score_threshold)
    }
    if (!is.null(FOVs)) {
        fovs_ <- as.integer(FOVs)
        a <- dplyr::filter(a, fov %in% fovs_)
    }
    if (isTRUE(flip_vertical)) {
        a <- dplyr::mutate(a, global_y = -global_y)
    }
    a
}

# detected_transcripts.csv converted once to a parquet dataset with the
# vendor's columns (an unnamed index column dropped, `cell_id` as string),
# cached by path, size and modification time like the zarr conversions.
# Tiling scans its source once per tile; a csv would be parsed in full each
# time.
.merscope_csv_to_parquet <- function(path, verbose = NULL) {
    info <- file.info(path)
    key <- .hash(list(
        path = normalizePath(path, winslash = "/", mustWork = FALSE),
        size = as.numeric(info$size),
        mtime = as.numeric(info$mtime),
        what = "merscope-transcripts-csv-v1"
    ))
    cache_dir <- file.path(getArtifactDumpDir(), "merscope_cache", key)
    target <- file.path(cache_dir, "detected_transcripts")
    if (dir.exists(target)) {
        GiottoUtils::vmsg(.v = verbose,
            "[merscope] using cached parquet transcripts:", target)
        return(target)
    }
    GiottoUtils::vmsg(.v = verbose,
        "[merscope] converting csv transcripts to parquet (once):", path)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    sch <- .merscope_tx_schema(path, format = "csv")
    q <- arrow::open_dataset(path, format = "csv", schema = sch$schema,
                             skip = 1L) |>
        dplyr::select(-dplyr::any_of(sch$dropcols))
    # written aside and renamed, so an interrupted conversion is never reused
    tmp <- file.path(cache_dir, ".writing")
    unlink(tmp, recursive = TRUE)
    arrow::write_dataset(q, tmp, format = "parquet")
    file.rename(tmp, target)
    target
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

    if (!grepl("\\.parquet$", path, ignore.case = TRUE)) {
        path <- .merscope_csv_to_parquet(path, verbose = verbose)
    }

    # capture bare values only -- read_fun is shipped to workers
    score_thr <- score_threshold
    fovs_     <- if (is.null(FOVs)) NULL else as.integer(FOVs)
    drop      <- as.character(dropcols)
    flip      <- isTRUE(flip_vertical)
    read_fun <- function(x, schema = NULL) {
        .merscope_tx_query(x, score_threshold = score_thr, FOVs = fovs_,
                           dropcols = drop, flip_vertical = flip)
    }

    fs <- fileStore(path = path, read_fun = read_fun)
    qs <- as(fs, "queryableStore")
    # the same query for the Rust writer (GiottoIngest), if installed:
    # unnamed (pandas index) columns are named "" in the file
    tx_cols <- names(arrow::open_dataset(path)$schema)
    qs@params$points_spec <- .ingest_points_spec(path, flip_y = flip,
        filters = c(
            if (!is.null(score_thr) && "transcript_score" %in% tx_cols) {
                list(list(col = "transcript_score", op = ">=", value = score_thr))
            },
            if (!is.null(fovs_)) list(list(col = "fov", op = "in", value = fovs_))
        ),
        to_string = "cell_id", drop = unique(c("", drop)))

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

    .tx_points_out(tx_store, output, feat_type, split_keyword)
}


## polygon ####

# The z-plane to load from cell_boundaries.parquet: the lowest of
# `poly_z_indices` (all planes when NULL), after checking that the planes are
# replicates. NA when the file has no ZIndex.
.merscope_zplane_architecture <- function(path, poly_z_indices = NULL,
                                          n_sample = 200L, verbose = NULL) {
    ZIndex <- EntityID <- Geometry <- NULL # NSE vars
    ds <- arrow::open_dataset(path)
    if (!"ZIndex" %in% names(ds)) {
        return(list(kind = "2d", z_indices = NA_integer_, use_z = NA_integer_))
    }
    ds <- dplyr::mutate(ds, EntityID = arrow::cast(EntityID, arrow::string()))

    z_all <- sort(ds |> dplyr::distinct(ZIndex) |> dplyr::collect() |>
                  dplyr::pull(ZIndex))
    z_all <- .merscope_select_z(z_all, poly_z_indices)
    if (length(z_all) == 1L) {
        GiottoUtils::vmsg(.v = verbose,
            "[merscope] single z-plane; 2D architecture")
        return(list(kind = "2d", z_indices = z_all, use_z = z_all[[1L]]))
    }

    zA <- z_all[[1L]]
    zB <- z_all[[2L]]
    ids_A <- ds |> dplyr::filter(ZIndex == zA) |> dplyr::select(EntityID) |>
             dplyr::collect() |> dplyr::pull(EntityID)
    ids_B <- ds |> dplyr::filter(ZIndex == zB) |> dplyr::select(EntityID) |>
             dplyr::collect() |> dplyr::pull(EntityID)
    if (!setequal(ids_A, ids_B)) {
        .merscope_stop_3d(z_all, reason = "cell ID sets differ between z-planes")
    }

    # compare the serialized geometry of a sample of cells across the planes
    samp <- .merscope_sample_ids(ids_A, n_sample)
    g <- ds |> dplyr::filter(EntityID %in% samp, ZIndex %in% z_all) |>
         dplyr::select(EntityID, ZIndex, Geometry) |> dplyr::collect()
    same <- vapply(split(g, g$EntityID), function(gi) {
        gi <- gi[order(gi$ZIndex), ]
        ref <- as.raw(gi$Geometry[[1L]])
        all(vapply(seq_len(nrow(gi)),
                   function(i) identical(as.raw(gi$Geometry[[i]]), ref),
                   logical(1L)))
    }, logical(1L))
    if (!all(same)) {
        .merscope_stop_3d(z_all, reason = sprintf(
            "geometry differs across z-planes for %d of %d sampled cells",
            sum(!same), length(same)))
    }

    GiottoUtils::vmsg(.v = verbose, sprintf(
        "[merscope] identical polygons across %d z-planes; using z-plane %s",
        length(z_all), zA))
    list(kind = "replicated_2d", z_indices = z_all, use_z = zA)
}

# z-planes to consider: all of `z_all`, or the requested ones that exist
.merscope_select_z <- function(z_all, poly_z_indices = NULL) {
    if (is.null(poly_z_indices)) return(z_all)
    keep <- sort(intersect(as.integer(poly_z_indices), as.integer(z_all)))
    if (!length(keep)) {
        stop("[merscope] none of `poly_z_indices` (",
             paste(poly_z_indices, collapse = ", "),
             ") are present. Available: ", paste(z_all, collapse = ", "),
             call. = FALSE)
    }
    keep
}

# up to n IDs spread evenly over `ids`: deterministic, leaves the RNG alone
.merscope_sample_ids <- function(ids, n) {
    ids <- sort(unique(ids))
    if (length(ids) <= n) return(ids)
    ids[unique(round(seq(1, length(ids), length.out = n)))]
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

# cell_boundaries.parquet -> SpatVector (`poly_ID`) of one z-plane. With
# `flip_vertical`, y is negated in the WKB itself (wk), before terra decodes
# it: no vertex matrix is built.
.merscope_poly_parquet <- function(path, poly_z_indices = NULL, ids = NULL,
                                   flip_vertical = TRUE, verbose = NULL) {
    ZIndex <- EntityID <- Geometry <- NULL # NSE vars
    checkmate::assert_file_exists(path)
    arch <- .merscope_zplane_architecture(path,
        poly_z_indices = poly_z_indices, verbose = verbose)

    q <- arrow::open_dataset(path) |>
        dplyr::mutate(EntityID = arrow::cast(EntityID, arrow::string()))
    if (!is.na(arch$use_z)) {
        use_z <- arch$use_z
        q <- dplyr::filter(q, ZIndex == use_z)
    }
    if (!is.null(ids)) q <- dplyr::filter(q, EntityID %in% ids)
    x <- q |> dplyr::select(EntityID, Geometry) |> dplyr::collect()
    if (nrow(x) == 0L) return(NULL)
    GiottoUtils::vmsg(.v = verbose,
        sprintf("[merscope_poly_disk] decoding %s geometries",
                format(nrow(x), big.mark = ",")))

    wkb <- as.list(x$Geometry)
    if (isTRUE(flip_vertical)) {
        GiottoUtils::package_check("wk")
        wkb <- unclass(wk::wk_transform(wk::new_wk_wkb(wkb),
            wk::wk_affine_scale(1, -1)))
    }
    sv <- terra::vect(wkb)
    terra::values(sv) <- data.frame(poly_ID = x$EntityID)
    sv
}

# cell_boundaries/feature_data_<fov>.hdf5 -> SpatVector (`poly_ID`) of one
# z-plane, read with Giotto's Vizgen HDF5 reader (raw vertices, not smoothed).
# Planes are 0-based like ZIndex; the HDF5 reader counts from 1. The replicate
# check reads two planes of one FOV file; all FOVs are then read for one plane.
.merscope_poly_hdf5 <- function(path, poly_z_indices = NULL, FOVs = NULL,
                                flip_vertical = TRUE, n_sample = 200L,
                                verbose = NULL) {
    GiottoUtils::package_check("rhdf5", repository = "Bioc")
    checkmate::assert_directory_exists(path)
    y <- NULL # NSE var
    files <- list.files(path, pattern = "^feature_data_\\d+\\.hdf5$")
    all_fovs <- sort(as.integer(sub("^feature_data_(\\d+)\\.hdf5$", "\\1",
        files)))
    sel <- if (is.null(FOVs)) all_fovs else intersect(all_fovs, FOVs)
    if (!length(sel)) return(NULL)

    read_planes <- function(fovs, z) {
        Giotto::readPolygonFilesVizgenHDF5(
            boundaries_path = path,
            fovs = fovs,
            z_indices = z + 1L,
            flip_y_axis = FALSE,
            calc_centroids = FALSE,
            H5Fopen_flags = "H5F_ACC_RDONLY",
            verbose = FALSE,
            output = "data.table"
        )
    }

    z <- .merscope_select_z(0:6, poly_z_indices)
    if (length(z) > 1L) {
        planes <- read_planes(sel[[1L]], z[1:2])
        if (length(planes) > 1L) {
            a <- planes[[1L]]
            b <- planes[[2L]]
            if (!setequal(unique(a$cell_id), unique(b$cell_id))) {
                .merscope_stop_3d(z, reason = "cell ID sets differ between z-planes")
            }
            samp <- .merscope_sample_ids(a$cell_id, n_sample)
            va <- a[a$cell_id %in% samp]
            vb <- b[b$cell_id %in% samp]
            va <- va[order(va$cell_id)]
            vb <- vb[order(vb$cell_id)]
            if (nrow(va) != nrow(vb) || !isTRUE(all.equal(va$x, vb$x)) ||
                !isTRUE(all.equal(va$y, vb$y))) {
                .merscope_stop_3d(z, reason = "geometry differs across z-planes")
            }
        }
    }
    v <- read_planes(sel, z[[1L]])
    if (!length(v)) return(NULL)
    v <- v[[1L]]
    if (isTRUE(flip_vertical)) v[, y := -y]
    ids <- unique(v$cell_id)
    m <- cbind(geom = match(v$cell_id, ids), part = 1, x = v$x, y = v$y,
               hole = 0)
    sv <- terra::vect(m, type = "polygons")
    terra::values(sv) <- data.frame(poly_ID = as.character(ids))
    sv
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
    cells = NULL,
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

    GiottoUtils::vmsg(
        sprintf("[merscope_poly_disk] loading boundary '%s'", name),
        .v = verbose
    )

    sv <- switch(polygon_format,
        "parquet" = {
            ids <- NULL
            if (!is.null(FOVs)) {
                if (is.null(cells)) {
                    cells <- .merscope_read_cells(metadata_path, FOVs = FOVs)
                }
                ids <- .merscope_fov_cell_ids(cells, FOVs)
            }
            .merscope_poly_parquet(path, poly_z_indices = poly_z_indices,
                ids = ids, flip_vertical = flip_vertical, verbose = verbose)
        },
        "hdf5" = .merscope_poly_hdf5(path, poly_z_indices = poly_z_indices,
            FOVs = FOVs, flip_vertical = flip_vertical, verbose = verbose)
    )
    if (is.null(sv)) {
        GiottoUtils::vmsg(.v = verbose, "[merscope_poly_disk] no polygons found")
        return(NULL)
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

# Disk-backed MERSCOPE expression ingestion. With `FOVs`, the store keeps the
# cells of those FOVs only, like every other data type.
.merscope_expression_disk <- function(
    path,
    gsource,
    spat_unit = "cell",
    split_by_type = TRUE,
    feat_type = c("rna", "neg_probe"),
    split_keyword = "Blank",
    FOVs = NULL,
    metadata_path = NULL,
    cells = NULL,
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

    if (!is.null(FOVs)) {
        if (is.null(cells)) {
            cells <- .merscope_read_cells(metadata_path, FOVs = FOVs)
        }
        ids <- .merscope_fov_cell_ids(cells, FOVs)
        pe <- pe[, which(colnames(pe) %in% ids), drop = FALSE]
    }

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

# cell_metadata.csv as a data.table with `cell_ID` as character (MERSCOPE
# EntityIDs are 64-bit integers, beyond what a double holds exactly). Older
# exports have an unnamed first column holding the IDs.
.merscope_read_cells <- function(path, FOVs = NULL) {
    if (is.null(path) || is.na(path)) {
        stop("[merscope] cell_metadata.csv is needed here (to map cells to ",
             "FOVs, or for cell centers), but none was found.", call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    hdr <- names(data.table::fread(path, nrows = 0L))
    id_col <- if ("EntityID" %in% hdr) "EntityID" else hdr[[1L]]
    dt <- data.table::fread(path,
        colClasses = stats::setNames("character", id_col))
    data.table::setnames(dt, id_col, "cell_ID")
    if (!is.null(FOVs) && "fov" %in% names(dt)) {
        dt <- dt[dt[["fov"]] %in% as.integer(FOVs)]
    }
    dt
}

# cell IDs of the selected FOVs
.merscope_fov_cell_ids <- function(cells, FOVs) {
    ids <- cells[["cell_ID"]][cells[["fov"]] %in% as.integer(FOVs)]
    if (!length(ids)) {
        stop("[merscope] no cells found in FOV(s): ",
             paste(FOVs, collapse = ", "), call. = FALSE)
    }
    ids
}

.merscope_cellmeta <- function(path, spat_unit = "cell", FOVs = NULL,
                               cells = NULL, verbose = NULL) {
    GiottoUtils::vmsg("[merscope_cellmeta] loading cell metadata", .v = verbose)
    if (is.null(cells)) cells <- .merscope_read_cells(path, FOVs = FOVs)
    GiottoClass::createCellMetaObj(
        metadata = cells,
        spat_unit = spat_unit,
        feat_type = "rna",
        provenance = spat_unit
    )
}


## spatial locations ####

.merscope_spatlocs <- function(path, spat_unit = "cell", FOVs = NULL,
                               cells = NULL, flip_vertical = TRUE,
                               verbose = NULL) {
    GiottoUtils::vmsg("[merscope_spatlocs] loading cell centers", .v = verbose)
    if (is.null(cells)) cells <- .merscope_read_cells(path, FOVs = FOVs)
    sdimy <- NULL # NSE var
    sl <- data.table::data.table(
        sdimx = cells[["center_x"]],
        sdimy = cells[["center_y"]],
        cell_ID = cells[["cell_ID"]]
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
