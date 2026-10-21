#' @include convenience-xenium.R reader-shared.R
NULL

# VisiumHD ingest pipeline for `gDirSource`-managed projects.
#
# Disk-backed counterpart to Giotto's `VisiumHDReader`. Expression goes
# through the shared 10x path (`.tenx_expression_disk`). Everything else is
# the inherited closures: tissue positions, cell metadata, scalefactors and
# images read `tissue_positions.parquet` lazily through arrow and apply every
# spatial filter, micron scaling and the y flip; the bin points
# (`load_transcripts`, the 2 um bins by default) load as a giottoBinPoints
# (points and counts in the vault); polygons and tessellations
# are built in memory and written to the vault by the backed object's setter
# when attached, as in StereoSeqDiskReader.
#
# The object takes the same shape as the disk Stereo-seq one -- bin
# spat_units named as in memory (`bin002`/`bin008`/`bin016`, `cell`),
# `rna`/`raw` expression, `raw` spatlocs, bin points -- and
# several units can be loaded into one object: `create_gobject()` takes a
# vector of bins and an existing `gobject` to add into, and
# `load_bin_mapping` records each bin's parent units as cell metadata.



# CLASS ####



setClass(
    "VisiumHDDiskReader",
    contains = "VisiumHDReader",
    slots = list(
        backend = "ANY"
    ),
    prototype = list(
        backend = NULL
    )
)

# * init ####
setMethod(
    "initialize", signature("VisiumHDDiskReader"),
    function(.Object, ..., backend) {
        obj <- callNextMethod(.Object, ...)

        if (!missing(backend)) {
            if (is.character(backend)) {
                backend <- gDirSource(path = backend)
            }
            checkmate::assert_class(backend, "gsource")
            obj@backend <- backend
        }
        if (is.null(obj@backend)) {
            stop("[VisiumHDDiskReader] `backend` is required", call. = FALSE)
        }
        # parent returns before path detection when no directory is set
        if (length(obj@visiumhd_dir) == 0L) return(obj)

        list2env(obj@paths, envir = environment())
        # the parent's own call: it builds polygon loaders only for segmented
        # outputs. (`binpath2`, the 2 um output, may be absent for either.)
        is_seg <- !is.null(obj@calls$load_polygon)

        # Distinct names for detected paths used as `gobject_fun` defaults;
        # `binpath = binpath` in a formal would be a recursive default
        # (see StereoSeqDiskReader).
        .def_binpath <- binpath
        .def_binpath2 <- if (is_seg) NULL else obj@paths$binpath2
        .def_map_path <- if (is_seg) NULL
                         else .visiumhd_mapping_path(obj@visiumhd_dir, binpath)

        gsrc <- obj@backend
        parent <- obj@calls

        ## expression (disk) ####
        ex_fun <- function(
            path = .def_binpath,
            feature_id_type = obj@feature_id_type,
            remove_zero_rows = TRUE,
            split_by_type = TRUE,
            expression_source = obj@expression_source,
            bin = obj@bin,
            barcodes = obj@barcodes,
            outdir = obj@outdir,
            force_untar = FALSE,
            untar_params = list(),
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            .visiumhd_expression_disk(
                path = path,
                gsource = gsrc,
                feature_id_type = feature_id_type,
                remove_zero_rows = remove_zero_rows,
                split_by_type = split_by_type,
                expression_source = expression_source,
                bin = bin,
                barcodes = barcodes,
                outdir = outdir,
                force_untar = force_untar,
                untar_params = untar_params,
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_expression <- ex_fun

        ## bin mapping ####
        map_fun <- function(
            path = .def_map_path,
            spat_unit = sprintf("bin%03d", obj@bin),
            barcodes = NULL,
            verbose = NULL
        ) {
            .visiumhd_bin_mapping(path = path, spat_unit = spat_unit,
                                  barcodes = barcodes, verbose = verbose)
        }
        if (!is_seg) obj@calls$load_bin_mapping <- map_fun

        ## bin points (disk) ####
        # The bins as points: their counts as a parquetExprStore and one point
        # per bin with counts, as a store-backed giottoBinPoints per feature
        # type. Takes the
        # arguments and defaults (the 2 um output) of the in-memory
        # load_transcripts(), which builds giottoBinPoints.
        tx_fun <- function(
            expr_path = .def_binpath2,
            feature_id_type = obj@feature_id_type,
            remove_zero_rows = TRUE,
            split_by_type = TRUE,
            expression_source = obj@expression_source,
            tissue_positions_path = .def_binpath2,
            scalefactors_path = .def_binpath2,
            micron = obj@micron,
            outdir = obj@outdir,
            bin = 2L,
            tissue_only = obj@tissue_only,
            barcodes = obj@barcodes,
            array_subset_row = obj@array_subset_row,
            array_subset_col = obj@array_subset_col,
            pxl_subset_row = obj@pxl_subset_row,
            pxl_subset_col = obj@pxl_subset_col,
            filter = obj@filter,
            filter_coverage_cutoff = obj@filter_coverage,
            force_untar = FALSE,
            untar_params = list(),
            verbose = NULL,
            ...
        ) {
            if (is.null(expr_path)) {
                stop("[VisiumHDDiskReader] no binned output found for the ",
                     "bin points", call. = FALSE)
            }
            sl <- parent$load_tissue_position(
                path = tissue_positions_path,
                bin = bin,
                micron = micron,
                scalefactors_path = scalefactors_path,
                barcodes = barcodes,
                tissue_only = tissue_only,
                array_subset_row = array_subset_row,
                array_subset_col = array_subset_col,
                pxl_subset_row = pxl_subset_row,
                pxl_subset_col = pxl_subset_col,
                filter = filter,
                filter_coverage_cutoff = filter_coverage_cutoff,
                outdir = outdir,
                force_untar = force_untar,
                output = "spatLocsObj",
                verbose = verbose
            )
            stores <- .visiumhd_expression_disk(
                path = expr_path,
                gsource = gsrc,
                feature_id_type = feature_id_type,
                remove_zero_rows = remove_zero_rows,
                split_by_type = split_by_type,
                expression_source = expression_source,
                bin = bin,
                barcodes = GiottoClass::spatIDs(sl),
                outdir = outdir,
                untar_params = untar_params,
                output = "store",
                verbose = verbose
            )
            # one bin points store, shared by the feature types
            bins <- NULL
            lapply(names(stores), function(ft) {
                bp <- .binpoints_disk(stores[[ft]], sl, gsrc,
                    feat_type = ft, bins = bins, verbose = verbose)
                bins <<- bp@spatial
                bp
            })
        }
        if (!is_seg) obj@calls$load_transcripts <- tx_fun


        ## create_gobject (binned) ####
        # Mirrors the parent's binned gobject_fun: same arguments, same load
        # order, the same barcode reconciliation between spatlocs and
        # expression. Differences: it initializes a backed giotto, adds into
        # `gobject` when one is given, accepts several `bin`s, and can attach
        # the bin hierarchy.
        bin_gobject_fun <- function(
            load_expression = TRUE,
            load_spatlocs = TRUE,
            load_metadata = TRUE,
            load_transcripts = FALSE,
            load_image = TRUE,
            load_bin_mapping = FALSE,
            create_tessellated_polys = FALSE,
            bin = obj@bin,
            micron = obj@micron,
            tissue_only = obj@tissue_only,
            barcodes = obj@barcodes,
            array_subset_row = obj@array_subset_row,
            array_subset_col = obj@array_subset_col,
            pxl_subset_row = obj@pxl_subset_row,
            pxl_subset_col = obj@pxl_subset_col,
            filter = obj@filter,
            filter_coverage_cutoff = obj@filter_coverage,
            expression_source = obj@expression_source,
            feature_id_type = obj@feature_id_type,
            expression_remove_zero_rows = TRUE,
            expression_split_by_type = TRUE,
            image_type = NULL,
            tessellate_shape = "hexagon",
            tessellate_shape_size = 400,
            tessellate_name = sprintf("%s%d",
                tessellate_shape, as.integer(tessellate_shape_size)),
            tissue_positions_path = .def_binpath,
            scalefactors_path = .def_binpath,
            expression_path = .def_binpath,
            image_path = .def_binpath,
            outdir = obj@outdir,
            force_untar = FALSE,
            untar_params = list(),
            gobject = NULL,
            instructions = NULL,
            verbose = NULL
        ) {
            a <- as.list(environment())

            # several units: one reader per bin, each adding into the same
            # object. Paths are re-detected per bin, and the image, the 2 um
            # bin points and any tessellation are loaded once.
            if (length(bin) > 1L) {
                a <- a[setdiff(names(a), c("tissue_positions_path",
                    "scalefactors_path", "expression_path", "image_path"))]
                return(.create_units(a, "bin", as.integer(bin),
                    once = c("load_image", "load_transcripts",
                             "create_tessellated_polys"),
                    create_for = function(b) {
                        rb <- obj
                        rb$bin <- b
                        rb@calls$create_gobject
                    }))
            }

            load_expression <- as.logical(load_expression)
            load_spatlocs <- as.logical(load_spatlocs)
            load_metadata <- as.logical(load_metadata)
            if (load_spatlocs && !load_expression) {
                stop("[VisiumHDDiskReader] load_spatlocs = TRUE requires ",
                     "load_expression = TRUE", call. = FALSE)
            }
            if (load_metadata && !load_expression) {
                stop("[VisiumHDDiskReader] load_metadata = TRUE requires ",
                     "load_expression = TRUE", call. = FALSE)
            }
            su <- sprintf("bin%03d", as.integer(bin))
            funs <- obj@calls
            # a forced re-extract happens once, before anything is read
            if (isTRUE(force_untar)) {
                .visiumhd_untar(expression_path, bin = bin, outdir = outdir,
                    force = TRUE, untar_params = untar_params,
                    verbose = verbose)
            }

            spat_filter <- list(
                array_subset_row = array_subset_row,
                array_subset_col = array_subset_col,
                pxl_subset_row = pxl_subset_row,
                pxl_subset_col = pxl_subset_col,
                filter = filter,
                filter_coverage_cutoff = filter_coverage_cutoff
            )

            sl <- NULL
            if (load_spatlocs) {
                sl <- do.call(funs$load_tissue_position, c(list(
                    path = tissue_positions_path,
                    bin = bin,
                    micron = micron,
                    scalefactors_path = scalefactors_path,
                    barcodes = barcodes,
                    tissue_only = tissue_only,
                    output = "spatLocsObj",
                    verbose = verbose
                ), spat_filter))
                barcodes <- GiottoClass::spatIDs(sl) %||% barcodes
            }

            expr_list <- NULL
            if (load_expression) {
                expr_list <- funs$load_expression(
                    path = expression_path,
                    bin = bin,
                    barcodes = barcodes,
                    feature_id_type = feature_id_type,
                    remove_zero_rows = expression_remove_zero_rows,
                    split_by_type = expression_split_by_type,
                    expression_source = expression_source,
                    outdir = outdir,
                    untar_params = untar_params,
                    verbose = verbose
                )
                barcodes <- GiottoClass::spatIDs(expr_list[[1L]]) %||%
                    barcodes
                if (!is.null(sl)) sl <- sl[barcodes]
            }

            # bin points on disk, with the arguments the parent's
            # gobject_fun gives its loader. Unlike the parent, `barcodes` is
            # not replaced by the 2 um barcodes afterwards: for bin != 2 that
            # filtered the cell metadata to nothing.
            tx_list <- NULL
            if (load_transcripts) {
                tx_list <- funs$load_transcripts(
                    bin = 2L,
                    micron = micron,
                    scalefactors_path = scalefactors_path,
                    tissue_only = FALSE,
                    barcodes = NULL,
                    pxl_subset_row = pxl_subset_row,
                    pxl_subset_col = pxl_subset_col,
                    filter = filter,
                    filter_coverage_cutoff = filter_coverage_cutoff,
                    array_subset_row = if (!is.null(array_subset_row)) {
                        array_subset_row * (bin / 2)
                    },
                    array_subset_col = if (!is.null(array_subset_col)) {
                        array_subset_col * (bin / 2)
                    },
                    feature_id_type = feature_id_type,
                    remove_zero_rows = expression_remove_zero_rows,
                    split_by_type = expression_split_by_type,
                    expression_source = expression_source,
                    verbose = verbose
                )
            }

            cmeta <- NULL
            if (load_metadata) {
                cmeta <- funs$load_cellmeta(
                    tissue_positions_path = tissue_positions_path,
                    bin = bin,
                    barcodes = barcodes,
                    tissue_only = tissue_only,
                    verbose = verbose
                )
                # Rows arrive in arrow scan order, which is not stable across
                # runs; put them in expression order, as `sl[barcodes]` does
                # for spatlocs above.
                if (!is.null(barcodes)) {
                    cdt <- cmeta[]
                    cmeta[] <- cdt[match(barcodes, cdt$cell_ID, nomatch = 0L)]
                }
            }

            gimg <- NULL
            if (load_image) {
                gimg <- funs$load_image(
                    path = image_path,
                    bin = bin,
                    image_type = image_type,
                    micron = micron,
                    scalefactors_path = scalefactors_path,
                    verbose = verbose
                )
            }

            tess_poly <- NULL
            if (create_tessellated_polys) {
                tess_poly <- funs$tessellate_polygon(
                    tissue_positions_path = tissue_positions_path,
                    shape = tessellate_shape,
                    shape_size = tessellate_shape_size,
                    name = tessellate_name,
                    bin = bin,
                    micron = micron,
                    scalefactors_path = scalefactors_path,
                    verbose = verbose
                )
            }

            g <- .backed_gobject_init(gobject, gsrc, instructions)
            if (!is.null(expr_list)) {
                g <- GiottoClass::setGiotto(g, expr_list, verbose = verbose)
            }
            if (!is.null(sl)) g <- GiottoClass::setGiotto(g, sl, verbose = verbose)
            if (!is.null(cmeta)) {
                g <- GiottoClass::setGiotto(g, cmeta, verbose = verbose)
            }
            if (!is.null(gimg)) {
                g <- GiottoClass::setGiotto(g, gimg, verbose = verbose)
            }
            if (!is.null(tess_poly)) {
                g <- GiottoClass::setGiotto(g, tess_poly, verbose = verbose)
            }
            if (!is.null(tx_list)) {
                g <- GiottoClass::setGiotto(g, tx_list, verbose = verbose)
            }
            if (isTRUE(load_bin_mapping) && load_expression) {
                g <- .add_parent_units(g,
                    spat_unit = su,
                    mapping = funs$load_bin_mapping(spat_unit = su,
                        barcodes = barcodes, verbose = verbose),
                    verbose = verbose
                )
            }
            g
        }

        ## create_gobject (segmented) ####
        # Mirrors the parent's segmented gobject_fun, plus `gobject`. The 2 um
        # points for a cell object come from createGiottoVisiumHDObjectCell(),
        # which loads them from the binned outputs as it does in memory.
        seg_gobject_fun <- function(
            load_expression = TRUE,
            load_polygons = c("cell", "nucleus"),
            graphclust_annotated = FALSE,
            load_image = TRUE,
            micron = obj@micron,
            barcodes = obj@barcodes,
            expression_source = obj@expression_source,
            feature_id_type = obj@feature_id_type,
            expression_remove_zero_rows = TRUE,
            expression_split_by_type = TRUE,
            image_type = NULL,
            scalefactors_path = .def_binpath,
            expression_path = .def_binpath,
            image_path = .def_binpath,
            geojson_path = .def_binpath,
            gobject = NULL,
            instructions = NULL,
            verbose = NULL,
            ...
        ) {
            not_used <- names(list(...))
            if (length(not_used) > 0L) {
                GiottoUtils::vmsg(.v = verbose, "[visiumHD] params:",
                    toString(not_used), "not used with segmentation outputs")
            }
            funs <- obj@calls

            poly_list <- NULL
            if (length(load_polygons) > 0L) {
                load_polygons <- match.arg(load_polygons,
                    choices = c("cell", "nucleus"), several.ok = TRUE)
                poly_list <- lapply(load_polygons, function(ptype) {
                    p <- parent$load_polygon(
                        path = geojson_path,
                        type = ptype,
                        graphclust_annotated = graphclust_annotated,
                        scalefactors_path = scalefactors_path,
                        micron = micron,
                        verbose = verbose
                    )
                    if (!is.null(barcodes)) p <- p[barcodes]
                    p
                })
            }

            expr_list <- NULL
            if (as.logical(load_expression)) {
                expr_list <- funs$load_expression(
                    path = expression_path,
                    barcodes = barcodes,
                    feature_id_type = feature_id_type,
                    remove_zero_rows = expression_remove_zero_rows,
                    split_by_type = expression_split_by_type,
                    expression_source = expression_source,
                    verbose = verbose
                )
            }

            gimg <- NULL
            if (as.logical(load_image)) {
                gimg <- funs$load_image(
                    path = image_path,
                    image_type = image_type,
                    micron = micron,
                    scalefactors_path = scalefactors_path,
                    verbose = verbose
                )
            }

            g <- .backed_gobject_init(gobject, gsrc, instructions)
            for (p in poly_list) {
                g <- GiottoClass::setGiotto(g, p, centroids_to_spatlocs = TRUE,
                                            verbose = verbose)
            }
            if (!is.null(expr_list)) {
                g <- GiottoClass::setGiotto(g, expr_list, verbose = verbose)
            }
            if (!is.null(gimg)) {
                g <- GiottoClass::setGiotto(g, gimg, verbose = verbose)
            }
            g
        }

        obj@calls$create_gobject <- if (is_seg) seg_gobject_fun
                                    else bin_gobject_fun
        obj
    }
)



# CREATE READER ####

#' @title Import a 10x Visium HD assay (disk-backed)
#' @name importVisiumHDDisk
#' @description
#' Disk-backed counterpart to [Giotto::importVisiumHD()]. Produces a
#' `VisiumHDDiskReader` whose loaders write into the `gDirSource`-managed
#' project vault:
#'
#' * expression (binned or segmented) streams from the 10x `.h5` (or mtx
#'   directory) into a `parquetExprStore`;
#' * cell / nucleus polygons (segmented outputs) and tessellations are
#'   written to the vault as parquet geometry when attached to the object;
#' * the bin points (`load_transcripts`, the 2 um bins by default) are a
#'   `giottoBinPoints` backed by two stores: one point per bin with counts,
#'   and the bins' counts.
#'
#' Tissue positions, cell metadata and images come from the
#' inherited `VisiumHDReader` loaders, so spat_unit names, barcodes,
#' coordinates and filters match the in-memory reader. `create_gobject()` additionally takes
#' `bin` as a vector (e.g. `c(2, 8, 16)`), a `gobject` to add into, and
#' `load_bin_mapping = TRUE` to record each unit's parent units (from
#' `barcode_mappings.parquet`) as cell metadata. `NNNum_outputs.tar` bin
#' outputs are extracted next to the tar (or into `outdir`) on first use.
#' @param visiumhd_dir Visium HD output directory: the `outs` directory, a
#'   `binned_outputs` / `square_XXXum` directory, or `segmented_outputs`.
#' @param backend a `gsource` (typically `gDirSource`) project backend, or a
#'   directory path.
#' @param bin,micron,outdir,expression_source,feature_id_type,tissue_only,barcodes,array_subset_row,array_subset_col,pxl_subset_row,pxl_subset_col,filter,filter_coverage_cutoff
#'   passed through to the parent `VisiumHDReader`; see
#'   [Giotto::importVisiumHD()].
#' @returns `VisiumHDDiskReader` object
#' @seealso [Giotto::importVisiumHD()] for the in-memory variant
#' @export
importVisiumHDDisk <- function(
    visiumhd_dir = NULL,
    backend,
    bin = 8,
    micron = FALSE,
    outdir = NULL,
    expression_source = "raw",
    feature_id_type = c("symbol", "ensembl"),
    tissue_only = FALSE,
    barcodes = NULL,
    array_subset_row = NULL,
    array_subset_col = NULL,
    pxl_subset_row = NULL,
    pxl_subset_col = NULL,
    filter = NULL,
    filter_coverage_cutoff = 0.5
) {
    if (missing(backend)) {
        stop("[importVisiumHDDisk] `backend` is required", call. = FALSE)
    }
    # Same argument handling as Giotto::importVisiumHD(), so adding
    # `backend =` to a call changes nothing about what is read.
    a <- list(Class = "VisiumHDDiskReader", backend = backend)
    if (!is.null(visiumhd_dir)) a$visiumhd_dir <- visiumhd_dir
    a$bin <- as.integer(bin)
    a$micron <- as.logical(micron)
    if (!is.null(outdir)) a$outdir <- outdir
    a$expression_source <- match.arg(expression_source, c("raw", "filtered"))
    a$feature_id_type <- match.arg(feature_id_type, c("symbol", "ensembl"))
    a$tissue_only <- as.logical(tissue_only)
    if (!is.null(barcodes)) a$barcodes <- barcodes
    if (!is.null(array_subset_row)) a$array_subset_row <- array_subset_row
    if (!is.null(array_subset_col)) a$array_subset_col <- array_subset_col
    if (!is.null(pxl_subset_row)) a$pxl_subset_row <- pxl_subset_row
    if (!is.null(pxl_subset_col)) a$pxl_subset_col <- pxl_subset_col
    if (!is.null(filter)) a$filter <- filter
    if (!is.null(filter_coverage_cutoff) && !is.na(filter_coverage_cutoff)) {
        a$filter_coverage <- filter_coverage_cutoff
    }
    do.call(methods::new, args = a)
}



# MODULAR ####


## expression ####

# Resolve the 10x matrix under a bin or segmented output directory and ingest
# it through the shared 10x path. The unit is read off the matrix name --
# `*_feature_cell_matrix` is segmented ("cell"), `*_feature_bc_matrix` a bin
# ("bin%03d") -- so a direct file path works as well as a directory.
.visiumhd_expression_disk <- function(
    path,
    gsource,
    feature_id_type = c("symbol", "ensembl"),
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    expression_source = c("raw", "filtered"),
    bin = 8L,
    barcodes = NULL,
    outdir = NULL,
    force_untar = FALSE,
    untar_params = list(),
    output = c("exprObj", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path) || is.null(path) || !file.exists(path)) {
        stop("[visiumhd_expression_disk] no expression path found",
             call. = FALSE)
    }
    feature_id_type <- match.arg(feature_id_type, c("symbol", "ensembl"))
    expression_source <- match.arg(expression_source)
    output <- match.arg(output)
    checkmate::assert_character(barcodes, null.ok = TRUE)

    path <- .visiumhd_untar(path, bin = bin, outdir = outdir,
        force = force_untar, untar_params = untar_params, verbose = verbose)
    mat <- .visiumhd_matrix_path(path, expression_source)
    spat_unit <- if (grepl("_feature_cell_matrix", basename(mat))) {
        "cell"
    } else {
        sprintf("bin%03d", as.integer(bin))
    }
    GiottoUtils::vmsg(.v = verbose, "[visiumhd_expression_disk]", spat_unit,
                      ":", mat)

    stores <- .tenx_expression_disk(
        path = mat,
        gsource = gsource,
        gene_ids = if (feature_id_type == "symbol") "symbols" else "ensembl",
        remove_zero_rows = remove_zero_rows,
        split_by_type = split_by_type,
        spat_unit = spat_unit,
        output = "store",
        verbose = verbose,
        ...
    )

    # barcode filter, applied after zero-feature removal as in memory
    if (!is.null(barcodes)) {
        stores <- lapply(stores, function(pe) {
            keep <- which(pe@cell_ids %in% barcodes)
            if (length(keep) < length(pe@cell_ids)) pe <- pe[, keep, drop = FALSE]
            pe
        })
    }
    if (output == "store") return(stores)

    lapply(names(stores), function(ft) {
        methods::new("exprObj",
            name       = "raw",
            exprMat    = stores[[ft]],
            spat_unit  = spat_unit,
            feat_type  = ft,
            provenance = spat_unit
        )
    })
}

# `path` is a matrix file / mtx directory, or an output directory holding
# one. The .h5 is preferred: tenxH5Input streams it by cell.
.visiumhd_matrix_path <- function(path, expression_source) {
    if (!dir.exists(path) || file.exists(file.path(path, "matrix.mtx.gz"))) {
        return(path)
    }
    stems <- paste0(expression_source,
                    c("_feature_bc_matrix", "_feature_cell_matrix"))
    cands <- file.path(path, c(paste0(stems, ".h5"), stems))
    hit <- cands[file.exists(cands)]
    if (!length(hit)) {
        stop("[visiumhd_expression_disk] no ", expression_source,
             " feature matrix under ", path, call. = FALSE)
    }
    hit[[1L]]
}


## tar outputs ####

# A `NNNum_outputs.tar` path becomes its extracted `square_NNNum` directory
# (next to the tar, or into `outdir`); any other path is returned as is.
.visiumhd_untar <- function(path, bin, outdir = NULL, force = FALSE,
                            untar_params = list(), verbose = NULL) {
    if (!any(GiottoUtils::file_extension(path) %in% "tar")) return(path)
    do.call(.visiumhd_untar_bin, c(list(path, outdir = outdir,
        bin = bin, force = force, verbose = verbose), untar_params))
}

# Extract a binned-output tar to `outdir/square_NNNum` unless that folder
# already holds more than the spatial info. GiottoDisk's own copy of the
# in-memory reader's helper of the same name, so no Giotto internal is called.
.visiumhd_untar_bin <- function(f, outdir = NULL, bin = 8, pattern = NULL,
                                force = FALSE, verbose = NULL, ...) {
    checkmate::assert_file_exists(f, extension = c("tar", "tar.gz"))
    checkmate::assert_character(pattern, null.ok = TRUE)
    bin <- as.integer(bin)
    if (is.null(outdir)) outdir <- dirname(f)
    # the tar holds one top folder, `square_NNNum`
    data_subdir <- file.path(outdir, sprintf("square_%03dum", bin))
    # spatial/ may have been extracted on its own before
    data_info <- list.files(data_subdir, include.dirs = TRUE, recursive = FALSE)
    if (!length(data_info) || setequal(data_info, "spatial")) force <- TRUE
    if (!force) return(data_subdir)

    files <- NULL
    if (!is.null(pattern)) {
        tar_manifest <- utils::untar(f, list = TRUE, ...)
        files <- tar_manifest[grepl(pattern = pattern, x = tar_manifest)]
    }
    GiottoUtils::vmsg(.v = verbose, "untarring...\n", f)
    utils::untar(f, exdir = outdir, files = files, ...)
    data_subdir
}


## bin mapping ####

# `barcode_mappings.parquet` sits in the `outs` directory, above
# binned_outputs/ and segmented_outputs/.
.visiumhd_mapping_path <- function(visiumhd_dir, binpath) {
    dirs <- unique(c(visiumhd_dir, dirname(binpath), dirname(dirname(binpath))))
    cands <- file.path(dirs, "barcode_mappings.parquet")
    hit <- cands[file.exists(cands)]
    if (length(hit)) hit[[1L]] else NULL
}

# Read the rows of `barcode_mappings.parquet` for one unit, with columns
# renamed to Giotto spat_units (`square_008um` -> `bin008`, `cell_id` ->
# `cell`). `cell_ID` holds the unit's own ids. The in_nucleus / in_cell flags
# describe 2 um bins, so they travel only with `bin002`.
.visiumhd_bin_mapping <- function(path, spat_unit, barcodes = NULL,
                                  verbose = NULL) {
    if (is.null(path) || !file.exists(path)) {
        stop("[visiumhd_bin_mapping] barcode_mappings.parquet not found",
             call. = FALSE)
    }
    ds <- arrow::open_dataset(path)
    cols <- names(ds)
    unit_cols <- grep("^square_\\d+um$", cols, value = TRUE)
    su_of <- function(x) {
        ifelse(x == "cell_id", "cell",
               sprintf("bin%03d", as.integer(gsub("\\D", "", x))))
    }
    key <- cols[su_of(cols) == spat_unit & cols %in% c(unit_cols, "cell_id")]
    if (length(key) != 1L) {
        stop("[visiumhd_bin_mapping] no mapping column for spat_unit '",
             spat_unit, "'", call. = FALSE)
    }
    parents <- setdiff(c(unit_cols, "cell_id"), key)
    flags <- if (identical(spat_unit, "bin002")) {
        intersect(c("in_nucleus", "in_cell"), cols)
    }
    GiottoUtils::vmsg(.v = verbose, "[visiumhd_bin_mapping]", spat_unit,
                      "<-", toString(su_of(parents)))

    q <- dplyr::select(ds, dplyr::all_of(c(key, parents, flags)))
    key_sym <- rlang::sym(key)
    q <- dplyr::filter(q, !is.na(!!key_sym))
    if (!is.null(barcodes)) {
        q <- dplyr::filter(q, !!key_sym %in% barcodes)
    }
    dt <- data.table::as.data.table(dplyr::collect(q))
    data.table::setnames(dt, c(key, parents), c("cell_ID", su_of(parents)))
    dt
}
