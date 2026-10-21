# Atera ingest pipeline for `gDirSource`-managed projects.
#
# A standalone reader for Atera output (zarr only): its own path detection
# and object assembly, so it can diverge from Xenium. It follows the same
# strategy as the other disk readers -- each data type is written to the
# vault by `sourceWrite()` as its own store, one piece at a time -- and
# reuses lower-level loaders where the data type and its reading are the
# same as Xenium's:
#   polygons      .tenx_poly_disk() on the cells archive (same as Xenium zarr)
#   expression    .tenx_expression_disk() (shared 10x path; the Atera matrix
#                 streams through tenxZarrCsrInput)
#   cell metadata Giotto's 10x cells reader, on the converted cells archive
#   images        Giotto's 10x OME-TIFF image loader
# Transcripts are Atera-specific: the zarr grid is streamed straight into
# tile partitions (tenxZarrTxInput, R/zarr-tx-store.R).

# CLASS ####

#' @title Atera disk reader
#' @name AteraDiskReader-class
#' @description
#' Disk-backed reader for Atera output. Create with [importAteraDisk()];
#' `$create_gobject()` builds the object and `$load_*()` load one data type.
#' @keywords internal
setClass("AteraDiskReader",
    slots = list(
        atera_dir = "character",
        qv = "ANY",
        micron = "numeric",
        paths = "list",
        calls = "list",
        backend = "ANY"
    ),
    prototype = list(
        qv = 20,
        backend = NULL
    )
)

# * paths ####

# The Atera v1 output files the reader uses, by exact name: the folder also
# holds archives it does not read (binned_transcripts, csc_cell_feature_matrix)
# and an instrument panel_config.json that is not a gene panel.
.atera_detect_paths <- function(dir) {
    exp_path <- file.path(dir, "experiment.spatial")
    img_dir <- file.path(dir, "morphology_2d")
    list(
        tx_path = .zarr_find_product(dir, "transcripts"),
        cells_path = .zarr_find_product(dir, "cells"),
        expr_path = .zarr_find_product(dir, "cell_feature_matrix"),
        experiment_path = if (file.exists(exp_path)) exp_path else NULL,
        img_paths = .atera_image_paths(img_dir)
    )
}

# Morphology images named by their channel (`ch0000_dapi.ome.tif` ->
# "dapi"), DAPI first. Empty when the folder is absent.
.atera_image_paths <- function(img_dir) {
    if (!dir.exists(img_dir)) return(list())
    f <- list.files(img_dir, pattern = "\\.ome\\.tif$", full.names = TRUE,
        ignore.case = TRUE)
    if (!length(f)) return(list())
    nm <- sub("^ch\\d+_(.*?)\\.ome\\.tif$", "\\1", basename(f),
        ignore.case = TRUE)
    ord <- order(!grepl("^dapi$", nm, ignore.case = TRUE))
    stats::setNames(as.list(f[ord]), nm[ord])
}

# Pixel size from experiment.spatial. NULL when absent.
.atera_pixel_size <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    jsonlite::fromJSON(path)$pixel_size
}

# * init ####
setMethod(
    "initialize", signature("AteraDiskReader"),
    function(.Object, ..., atera_dir, qv, micron, backend) {
        obj <- callNextMethod(.Object, ...)
        if (!missing(atera_dir)) obj@atera_dir <- atera_dir
        if (!missing(qv)) obj@qv <- qv
        if (!missing(micron)) obj@micron <- micron
        if (!missing(backend)) {
            if (is.character(backend)) backend <- gDirSource(path = backend)
            checkmate::assert_class(backend, "gsource")
            obj@backend <- backend
        }
        if (is.null(obj@backend)) {
            stop("[AteraDiskReader] `backend` is required", call. = FALSE)
        }
        if (length(obj@atera_dir) == 0L) return(obj)
        checkmate::assert_directory_exists(obj@atera_dir)

        obj@paths <- .atera_detect_paths(obj@atera_dir)
        tx <- obj@paths$tx_path
        if (is.null(tx) || !identical(
            .zarr_detect_product(tx, "transcripts"), "atera-zarr-v6")) {
            stop("[AteraDiskReader] no Atera v1 transcripts.zarr.zip in ",
                obj@atera_dir, ". Output in the Xenium layout (e.g. Atera ",
                "preview data) reads with importXeniumDisk().", call. = FALSE)
        }
        if (length(obj@micron) == 0L) {
            obj@micron <- .atera_pixel_size(obj@paths$experiment_path) %||% {
                warning(GiottoUtils::wrap_txt("No experiment.spatial found.",
                    "Guessing 0.2125 as micron scaling"), call. = FALSE)
                0.2125
            }
        }
        # distinct names for the path defaults below
        .tx_path <- obj@paths$tx_path
        .cells_path <- obj@paths$cells_path
        .expr_path <- obj@paths$expr_path
        .img_paths <- obj@paths$img_paths
        gsrc <- obj@backend

        ## transcripts (Atera-specific) ####
        tx_fun <- function(
            path = .tx_path,
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
            .atera_transcript_disk(
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

        ## polygons (shared 10x cells archive) ####
        poly_fun <- function(
            path = .cells_path,
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
            .tenx_poly_disk(
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

        ## expression (shared 10x path) ####
        ex_fun <- function(
            path = .expr_path,
            gene_ids = "symbols",
            remove_zero_rows = TRUE,
            split_by_type = TRUE,
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            .tenx_expression_disk(
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

        ## cell metadata and centroids (10x cells schema, in memory) ####
        cmeta_fun <- function(
            path = .cells_path,
            dropcols = c("x_centroid", "y_centroid"),
            verbose = NULL
        ) {
            do.call(Giotto:::.xenium_cellmeta, list(
                path = .zarr_ensure_parquet(path, what = "cells"),
                dropcols = dropcols, verbose = verbose
            ))
        }
        obj@calls$load_cellmeta <- cmeta_fun

        sloc_fun <- function(
            path = .cells_path,
            flip_vertical = TRUE,
            verbose = NULL
        ) {
            do.call(Giotto:::.xenium_spatlocs, list(
                path = .zarr_ensure_parquet(path, what = "cells"),
                flip_vertical = flip_vertical, verbose = verbose
            ))
        }
        obj@calls$load_spatlocs <- sloc_fun

        ## images (10x OME-TIFF) ####
        img_fun <- function(
            path,
            name = "image",
            output_dir,
            micron = obj@micron,
            negative_y = TRUE,
            flip_vertical = FALSE,
            flip_horizontal = FALSE,
            verbose = NULL,
            ...
        ) {
            a <- list(path = path, name = name, micron = micron,
                negative_y = negative_y, flip_vertical = flip_vertical,
                flip_horizontal = flip_horizontal, verbose = verbose, ...)
            if (!missing(output_dir)) a$output_dir <- output_dir
            do.call(Giotto:::.xenium_image, a)
        }
        obj@calls$load_image <- img_fun

        ## create_gobject ####
        ## aligned images (10x image alignment) ####
        img_aff_fun <- function(
            path,
            imagealignment_path,
            name = "aligned_image",
            micron = obj@micron,
            verbose = NULL
        ) {
            Giotto::read10xAffineImage(
                file = path,
                imagealignment_path = imagealignment_path,
                name = name,
                micron = micron,
                verbose = verbose
            )
        }
        obj@calls$load_aligned_image <- img_aff_fun

        ## create_gobject ####
        # Takes the arguments of the Giotto 10x readers' create_gobject so
        # Giotto's createGiottoAteraObject(backend = ...) can call it.
        # `gene_panel_json_path` is accepted for that reason; Atera ships no
        # gene panel.
        gobject_fun <- function(
            transcript_path = .tx_path,
            load_bounds = list(
                cell = "cell",
                nucleus = "nucleus"
            ),
            gene_panel_json_path = NULL,
            expression_path = .expr_path,
            metadata_path = .cells_path,
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
            load_images = "morphology",
            load_aligned_images = NULL,
            load_transcripts = TRUE,
            load_expression = TRUE,
            load_cellmeta = FALSE,
            instructions = NULL,
            verbose = NULL
        ) {
            if (!load_transcripts && !load_expression) {
                warning(GiottoUtils::wrap_txt("One of either transcripts or",
                    "expression info should be loaded for a fully",
                    "functioning object"))
            }
            if (length(gene_panel_json_path) &&
                nzchar(gene_panel_json_path[[1L]])) {
                warning("[AteraDiskReader] `gene_panel_json_path` is ignored: ",
                    "Atera ships no gene panel", call. = FALSE)
            }
            funs <- obj@calls
            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            if (load_transcripts) {
                tx_list <- funs$load_transcripts(path = transcript_path,
                    feat_type = feat_type, split_keyword = split_keyword,
                    verbose = verbose)
                g <- GiottoClass::setGiotto(g, tx_list, verbose = FALSE)
            }

            # named list: polygon name -> path, where "cell" / "nucleus"
            # stand for the cells archive (a character vector of those
            # shortnames works too)
            if (length(load_bounds)) {
                if (is.null(names(load_bounds))) {
                    names(load_bounds) <- unlist(load_bounds)
                }
                load_bounds <- as.list(load_bounds)
                is_short <- unlist(load_bounds) %in% c("cell", "nucleus")
                load_bounds[is_short] <- list(.cells_path)
                blist <- lapply(names(load_bounds), function(b) {
                    if (b == "nucleus") {
                        funs$load_polys(path = load_bounds[[b]],
                            name = "nucleus", split_geom = TRUE,
                            split_geom_fmt = "nucleus_%d", verbose = verbose)
                    } else {
                        funs$load_polys(path = load_bounds[[b]], name = b,
                            verbose = verbose)
                    }
                })
                # centroids are attached as spatlocs at set time
                g <- GiottoClass::setGiotto(g, do.call(c, blist),
                    centroids_to_spatlocs = TRUE, verbose = FALSE)
            } else {
                sl <- funs$load_spatlocs(path = metadata_path,
                    verbose = verbose)
                if (!is.null(sl)) g <- GiottoClass::setGiotto(g, sl,
                    verbose = FALSE)
            }

            if (load_expression) {
                ex <- funs$load_expression(path = expression_path,
                    verbose = verbose)
                g <- GiottoClass::setGiotto(g, ex)
            }

            if (load_cellmeta) {
                cx <- funs$load_cellmeta(path = metadata_path,
                    verbose = verbose)
                g <- GiottoClass::setGiotto(g, cx, verbose = FALSE)
            }

            # "morphology" (or Xenium's "focus") = every morphology_2d
            # channel; else a named list of image paths
            if (length(load_images) == 1L && is.character(load_images) &&
                load_images %in% c("morphology", "focus")) {
                load_images <- .img_paths
            }
            if (length(load_images)) {
                load_images <- as.list(load_images)
                checkmate::assert_list(load_images, names = "named")
                imglist <- do.call(c, lapply(names(load_images), function(nm) {
                    funs$load_image(path = load_images[[nm]], name = nm,
                        verbose = verbose)
                }))
                g <- GiottoClass::setGiotto(g, imglist)
            }

            # list of c(image path, alignment matrix path), by name
            if (length(load_aligned_images)) {
                checkmate::assert_list(load_aligned_images, names = "named")
                aimglist <- lapply(names(load_aligned_images), function(nm) {
                    funs$load_aligned_image(
                        path = load_aligned_images[[nm]][1L],
                        imagealignment_path = load_aligned_images[[nm]][2L],
                        name = nm
                    )
                })
                g <- GiottoClass::setGiotto(g, aimglist)
            }
            g
        }
        obj@calls$create_gobject <- gobject_fun

        obj
    }
)

# transcripts ####

# Atera v1 transcripts -> parquetGeomTileStore, streamed from the zarr grid
# (no intermediate table, no re-tiling).
.atera_transcript_disk <- function(
    path,
    gsource,
    feat_type,
    split_keyword,
    flip_vertical = TRUE,
    dropcols = c(),
    qv_threshold = 20,
    cores = GiottoUtils::determine_cores(),
    output = c("giottoPoints", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path) || is.null(path)) {
        stop("[atera_transcript_disk] no transcripts.zarr.zip found",
            call. = FALSE)
    }
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("giottoPoints", "store"))
    GiottoUtils::vmsg("[atera_transcript_disk] loading transcripts",
        .v = verbose)
    tx_store <- sourceWrite(
        gsource, tenxZarrTxInput(path),
        store_type = "parquetGeomTile",
        qv_threshold = qv_threshold,
        flip_vertical = flip_vertical,
        dropcols = dropcols,
        workers = cores,
        verbose = verbose,
        ...
    )
    .tx_points_out(tx_store, output, feat_type, split_keyword)
}

# access ####

#' @export
setMethod("$", signature("AteraDiskReader"), function(x, name) {
    if (name %in% c("atera_dir", "qv", "micron")) {
        return(methods::slot(x, name))
    }
    x@calls[[name]]
})

#' @export
setMethod("$<-", signature("AteraDiskReader"), function(x, name, value) {
    if (name %in% c("atera_dir", "qv")) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }
    stop("Only 'atera_dir' and 'qv' can be set", call. = FALSE)
})

#' @export
`.DollarNames.AteraDiskReader` <- function(x, pattern) {
    dn <- c("atera_dir", "qv", "micron")
    if (length(x@calls)) dn <- c(dn, paste0(names(x@calls), "()"))
    utils::.DollarNames(dn, pattern)
}

setMethod("show", signature("AteraDiskReader"), function(object) {
    cat(sprintf("GiottoDisk <%s>\n", class(object)[1L]))
    d <- object@atera_dir
    cat("dir       :", if (length(d)) GiottoUtils::str_abbreviate(d) else "",
        "\n")
    cat("qv_cutoff :", paste(object@qv, collapse = ", "), "\n")
    cat("micron    :", paste(object@micron, collapse = ", "), "\n")
    if (length(object@calls)) {
        cat("funs      :", paste0(names(object@calls), "()",
            collapse = ", "), "\n")
    }
})

#' @title Import an Atera assay (disk-backed)
#' @name importAteraDisk
#' @description
#' Disk-backed reader for Atera output (zarr archives). Each data type is
#' written to a `gDirSource`-managed project vault as its own store:
#' transcripts stream from the zarr grid straight into a tiled store, the
#' expression matrix streams in cell order, and cell / nucleus polygons go
#' through the same tiled polygon ingest as Xenium. The pixel size is read
#' from `experiment.spatial`.
#' @param atera_dir Atera output directory
#' @param backend a `gsource` (typically `gDirSource`) project backend, or a
#'   path to one
#' @param qv_threshold minimum Phred-scaled quality score retained for
#'   transcripts
#' @returns `AteraDiskReader` object
#' @examples
#' \dontrun{
#' rdr <- importAteraDisk("/path/to/atera_outs", backend = "my_project")
#' g <- rdr$create_gobject()
#' }
#' @export
importAteraDisk <- function(atera_dir = NULL, backend, qv_threshold = 20) {
    if (missing(backend)) {
        stop("[importAteraDisk] `backend` is required", call. = FALSE)
    }
    a <- list(Class = "AteraDiskReader", backend = backend, qv = qv_threshold)
    if (!is.null(atera_dir)) a$atera_dir <- atera_dir
    do.call(methods::new, args = a)
}
