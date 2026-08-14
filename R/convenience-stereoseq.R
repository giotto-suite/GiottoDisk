# Stereo-seq ingest pipeline for `gDirSource`-managed projects.
#
# Disk-backed counterpart to Giotto's `StereoSeqReader`. Routes the
# expression matrix through GiottoDisk (`parquetExprStore` written into the
# project vault via `sourceWrite`). Bin spatial locations come out of that
# same streaming pass rather than a second read. The remaining modalities
# (images, masks, binpoints, polygons) stay on the inherited in-mem closures
# and can be ported following the same pattern when needed.



# CLASS ####



setClass(
    "StereoSeqDiskReader",
    contains = "StereoSeqReader",
    slots = list(
        backend = "ANY"
    ),
    prototype = list(
        backend = NULL
    )
)

# * init ####
setMethod(
    "initialize", signature("StereoSeqDiskReader"),
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
            stop("[StereoSeqDiskReader] `backend` is required",
                 call. = FALSE)
        }

        # Mirror @paths into the init frame so closures can reference
        # parent-detected paths (gef_path, bin1_gef_path, image_dir)
        # directly via default-arg expressions, same convention as
        # XeniumDiskReader.
        list2env(obj@paths, envir = environment())

        # Distinct names for the detected paths that `gobject_fun` takes as
        # parameters. `gef_path = gef_path` in its formals is a recursive
        # default argument reference -- the promise resolves to the parameter
        # itself, not to the binding above -- so the path never arrives.
        # Giotto's own reader sidesteps this the same way, with
        # `.default_mask_path2`.
        .def_gef_path  <- gef_path
        .def_bin1_path <- bin1_gef_path
        .def_image_dir <- image_dir
        .def_mask_path <- mask_path

        gsrc <- obj@backend
        type_       <- obj@type
        bin_size_   <- obj@bin_size
        gene_col_   <- obj@gene_column

        # expression (disk override)
        ex_fun <- function(
            path        = gef_path,
            type        = type_,
            bin_size    = bin_size_,
            gene_column = gene_col_,
            spat_unit   = if (type_ == "bin") bin_size_ else "cell",
            output      = c("exprObj", "store"),
            verbose     = NULL,
            ...
        ) {
            .stereoseq_expression_disk(
                path        = path,
                gsource     = gsrc,
                type        = type,
                bin_size    = bin_size,
                gene_column = gene_column,
                spat_unit   = spat_unit,
                output      = output,
                verbose     = verbose,
                ...
            )
        }
        obj@calls$load_expression <- ex_fun

        # create_gobject (disk variant). Mirrors the parent's gobject_fun
        # body but: (a) initializes the giotto via createGiottoObject(
        # backend = gsrc, ...) instead of bare giotto(), and (b) the
        # expression branch goes through funs$load_expression -- the disk
        # override -- so the parquetExprStore lands in the vault.
        # Spatlocs / image / binpoints / polygons / mask remain inherited
        # (in-mem) closures and can be ported following the same pattern.
        gobject_fun <- function(
            load_expression = TRUE,
            load_spatlocs   = TRUE,
            load_binpoints  = FALSE,
            load_image      = TRUE,
            load_mask       = TRUE,
            load_polygons   = (obj@type == "cell"),
            type            = obj@type,
            bin_size        = obj@bin_size,
            gene_column     = obj@gene_column,
            negative_y      = obj@negative_y,
            gef_path        = .def_gef_path,
            bin1_path       = .def_bin1_path,
            image_path      = .def_image_dir,
            mask_path       = .def_mask_path,
            instructions    = NULL,
            verbose         = NULL) {

            spat_unit <- if (type == "bin") bin_size else "cell"

            # Match parent's guard: spatlocs needs expression coords
            if (load_spatlocs && !load_expression) {
                GiottoUtils::vmsg(.v = verbose,
                    "[StereoSeqDiskReader] load_spatlocs = TRUE ignored:",
                    "requires load_expression = TRUE.")
                load_spatlocs <- FALSE
            }

            funs <- obj@calls

            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            # expression (disk; overridden closure)
            bin_coords <- NULL
            if (load_expression) {
                ex <- funs$load_expression(
                    path        = gef_path,
                    type        = type,
                    bin_size    = bin_size,
                    gene_column = gene_column,
                    spat_unit   = spat_unit,
                    verbose     = verbose
                )
                bin_coords <- attr(ex, "bin_coords")
                g <- GiottoClass::setGiotto(g, ex, verbose = verbose)
            }

            # spatlocs. For bins, the ingest stream already produced the
            # (x, y) -> bin_ID map, so build from that. The inherited closure
            # would otherwise read the whole `geneExp/<bin>/expression`
            # dataset into memory a second time -- bin coordinates only exist
            # inside the expression records -- which is exactly the read the
            # disk backend exists to avoid. Cellbin keeps the inherited path:
            # its coordinates come from the small `cellBin/cell` table.
            if (load_spatlocs) {
                sl <- if (!is.null(bin_coords)) {
                    .stereoseq_spatlocs_from_coords(
                        bin_coords = bin_coords,
                        negative_y = negative_y,
                        spat_unit  = spat_unit,
                        verbose    = verbose
                    )
                } else {
                    funs$load_spatlocs(
                        path        = gef_path,
                        type        = type,
                        bin_size    = bin_size,
                        gene_column = gene_column,
                        negative_y  = negative_y,
                        spat_unit   = spat_unit,
                        verbose     = verbose
                    )
                }
                g <- GiottoClass::setGiotto(g, sl, verbose = verbose)
            }

            # binpoints (inherited; in-mem)
            if (load_binpoints) {
                if (is.null(bin1_path)) {
                    stop("[StereoSeqDiskReader] no bin GEF for bin1 binpoints. ",
                         "Provide `bin1_path` explicitly or ensure a *.tissue.gef",
                         " or *.gef is present.", call. = FALSE)
                }
                gbp <- funs$load_binpoints(
                    path        = bin1_path,
                    gene_column = gene_column,
                    negative_y  = negative_y,
                    spat_unit   = spat_unit,
                    verbose     = verbose
                )
                g <- GiottoClass::setGiotto(g, gbp, verbose = verbose)
            }

            # image (inherited; in-mem)
            if (load_image) {
                gimg <- funs$load_image(
                    path = image_path, negative_y = negative_y,
                    verbose = verbose
                )
                if (!is.null(gimg)) {
                    g <- GiottoClass::setGiotto(g, gimg, verbose = verbose)
                }
            }

            # polygons (cell-type only) / mask
            gpoly <- NULL
            if (type == "cell" && load_polygons) {
                gpoly <- funs$load_polygons(
                    path = gef_path, negative_y = negative_y,
                    verbose = verbose
                )
            } else if (load_mask) {
                if (is.null(mask_path)) {
                    warning("[StereoSeqDiskReader] no mask file provided; ",
                            "skipping mask polygons.", call. = FALSE)
                } else {
                    gpoly <- funs$load_mask(
                        path = mask_path, negative_y = negative_y,
                        verbose = verbose
                    )
                }
            }
            if (!is.null(gpoly)) {
                # Drop the terra centroids the in-memory loaders pre-compute,
                # so they are derived from the `parquetGeomBase` at attach
                # time instead. This is the same choice the Xenium disk
                # reader makes (see `.xenium_polygons_disk`). Keeping them
                # leaves a centroid SpatVector with no `poly_ID` attribute
                # next to a store-backed `@spatVector`, and subsetting the
                # polygon -- `subset(gpolygon@spatVectorCentroids, poly_ID
                # %in% cell_ids)` in GiottoClass -- then fails, taking
                # `filterGiotto()` / `subsetGiotto()` down with it.
                gpoly@spatVectorCentroids <- NULL
                g <- GiottoClass::setGiotto(g, gpoly, verbose = verbose)
            }

            g
        }
        obj@calls$create_gobject <- gobject_fun

        obj
    }
)



# CREATE READER ####

#' @title Import a BGI Stereo-seq assay (disk-backed)
#' @name importStereoSeqDisk
#' @description
#' Disk-backed counterpart to [Giotto::importStereoSeq()]. Produces a
#' `StereoSeqDiskReader` whose `load_expression()` call streams the
#' source `.gef` file into a `parquetExprStore` written to the
#' `gDirSource`-managed project vault. For `type = "bin"`, spatial
#' locations are derived from that same streaming pass. The remaining
#' modalities (images, masks, binpoints, polygons, and cellbin spatial
#' locations) come from the inherited `StereoSeqReader` closures.
#' @param stereoseq_dir Stereo-seq output directory
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#'   Naming matches [GiottoClass::createGiottoObject()]'s `backend` param.
#' @param type,bin_size,gene_column,negative_y,gef_type passed through to
#'   the parent `StereoSeqReader` initializer.
#' @returns `StereoSeqDiskReader` object
#' @seealso [Giotto::importStereoSeq()] for the in-memory variant
#' @export
importStereoSeqDisk <- function(
    stereoseq_dir = NULL,
    backend,
    type        = c("bin", "cell"),
    bin_size    = "bin100",
    gene_column = c("geneName", "geneID"),
    negative_y  = TRUE,
    gef_type
) {
    if (missing(backend)) {
        stop("[importStereoSeqDisk] `backend` is required", call. = FALSE)
    }
    type <- match.arg(type)
    gene_column <- match.arg(gene_column)
    # Must match Giotto::importStereoSeq()'s defaults. Switching a call from
    # the in-mem reader to this one by adding `backend =` should not change
    # which .gef gets read.
    if (missing(gef_type)) {
        gef_type <- if (type == "bin") "tissue" else "adjusted_cellbin"
    }
    a <- list(
        Class       = "StereoSeqDiskReader",
        backend     = backend,
        type        = type,
        bin_size    = bin_size,
        gene_column = gene_column,
        negative_y  = as.logical(negative_y),
        gef_type    = gef_type
    )
    if (!is.null(stereoseq_dir)) a$stereoseq_dir <- stereoseq_dir
    do.call(new, args = a)
}



# MODULAR ####


## expression ####

# Disk-backed Stereo-seq expression ingestion. Picks the appropriate
# exprInput marker based on `type`:
#   - "cell"  -> cellbinGefInput  (reads cellBin/cell + cellBin/gene)
#   - "bin"   -> binGefInput      (reads geneExp/<bin>/gene + ../expression)
# Routes through sourceWrite(gsource, inp, store_type = "parquetExpr").
# The .gef file is never fully materialized -- the iterator streams
# gene-chunks via rhdf5 hyperslabs.
.stereoseq_expression_disk <- function(
    path,
    gsource,
    type        = c("bin", "cell"),
    bin_size    = "bin100",
    gene_column = c("geneName", "geneID"),
    spat_unit   = NULL,
    output      = c("exprObj", "store"),
    verbose     = NULL,
    ...
) {
    if (missing(path) || is.null(path) || length(path) == 0L ||
        !file.exists(path)) {
        stop("[stereoseq_expression_disk] no .gef path provided",
             call. = FALSE)
    }
    checkmate::assert_class(gsource, "gsource")
    type <- match.arg(type)
    gene_column <- match.arg(gene_column)
    output <- match.arg(output, choices = c("exprObj", "store"))

    if (is.null(spat_unit)) {
        spat_unit <- if (type == "bin") bin_size else "cell"
    }

    GiottoUtils::vmsg("[stereoseq_expression_disk] type:", type,
                       " gef:", path, .v = verbose)

    if (type == "cell") {
        inp <- cellbinGefInput(path, gene_column = gene_column)
    } else {
        # Passed through verbatim: the GEF group key carries the "bin"
        # prefix (`geneExp/bin100/...`), matching StereoSeqReader's
        # `$bin_size` and the hardcoded "geneExp/bin1/expression" that
        # Giotto's binpoints reader uses.
        inp <- binGefInput(path,
                            bin_size    = as.character(bin_size),
                            gene_column = gene_column)
    }

    pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                       verbose = verbose, ...)

    # Bin coordinates live inside the expression records, so the ingest pass
    # is the only one that sees them for free. binGefInput's iterator
    # accumulates the (x, y) -> bin_ID map and publishes it here; carried out
    # as an attribute so `gobject_fun` can build spatlocs without a second
    # full read of the gef. NULL for cellbin, whose coordinates come from the
    # small cellBin/cell table instead.
    bin_coords <- inp@params$coord_env$bin_coords

    out <- if (output == "store") {
        pe
    } else {
        # Wrapped in a list to match what the in-memory
        # `StereoSeqReader$load_expression()` hands back. `setGiotto()` takes
        # either, but a reader used piecewise -- as the importer vignette
        # does -- must be substitutable with `backend =` set or unset.
        list(methods::new("exprObj",
            name       = "raw",
            exprMat    = pe,
            spat_unit  = spat_unit,
            feat_type  = "rna",
            provenance = spat_unit
        ))
    }

    attr(out, "bin_coords") <- bin_coords
    out
}


# Build a spatLocsObj from the (x, y) -> bin_ID map captured during ingest.
# Mirrors `.stereoseq_build_spatlocs()`'s bin branch in Giotto exactly --
# same `bin_<id>` naming, same y flip -- so IDs line up with the expression
# store column names whichever path produced them.
.stereoseq_spatlocs_from_coords <- function(
    bin_coords, negative_y = TRUE, spat_unit = NULL, verbose = NULL
) {
    x <- y <- bin_ID <- cell_ID <- NULL  # data.table vars

    spat_locs <- data.table::copy(bin_coords)
    data.table::setorder(spat_locs, bin_ID)
    spat_locs[, cell_ID := paste0("bin_", bin_ID)]
    spat_locs <- spat_locs[, .(cell_ID, x, y)]
    spat_locs[, x := as.integer(x)]
    spat_locs[, y := as.integer(y)]
    if (isTRUE(negative_y)) {
        spat_locs[, y := 0L - y]
    }

    GiottoUtils::vmsg(.v = verbose,
        "[stereoseq_spatlocs_disk]", nrow(spat_locs),
        "bins from the ingest stream")

    GiottoClass::createSpatLocsObj(
        coordinates = spat_locs,
        name        = "raw",
        spat_unit   = spat_unit,
        provenance  = spat_unit,
        verbose     = FALSE
    )
}
