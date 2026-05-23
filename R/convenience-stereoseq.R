# Stereo-seq ingest pipeline for `gDirSource`-managed projects.
#
# Disk-backed counterpart to Giotto's `StereoSeqReader`. Currently routes
# only the expression matrix through GiottoDisk (`parquetExprStore`
# written into the project vault via `sourceWrite`). Other modalities
# (spatlocs, images, masks, binpoints, polygons) remain on the inherited
# in-mem closures and can be ported following the same pattern when
# needed.



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
            gef_path        = gef_path,
            bin1_path       = bin1_gef_path,
            image_path      = image_dir,
            mask_path       = mask_path,
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
            if (load_expression) {
                ex <- funs$load_expression(
                    path        = gef_path,
                    type        = type,
                    bin_size    = bin_size,
                    gene_column = gene_column,
                    spat_unit   = spat_unit,
                    verbose     = verbose
                )
                g <- GiottoClass::setGiotto(g, ex, verbose = verbose)
            }

            # spatlocs (inherited; in-mem; reads gef separately from the
            # disk expression path, so two reads when both are requested)
            if (load_spatlocs) {
                sl <- funs$load_spatlocs(
                    path        = gef_path,
                    type        = type,
                    bin_size    = bin_size,
                    gene_column = gene_column,
                    negative_y  = negative_y,
                    spat_unit   = spat_unit,
                    verbose     = verbose
                )
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
#' `gDirSource`-managed project vault. Other modalities (spatlocs,
#' images, masks, binpoints, polygons) remain in-memory via the
#' inherited `StereoSeqReader` closures.
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
    bin_size    = "bin50",
    gene_column = c("geneName", "geneID"),
    negative_y  = TRUE,
    gef_type
) {
    if (missing(backend)) {
        stop("[importStereoSeqDisk] `backend` is required", call. = FALSE)
    }
    type <- match.arg(type)
    gene_column <- match.arg(gene_column)
    if (missing(gef_type)) {
        gef_type <- if (type == "bin") "tissue" else "cellbin"
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
    bin_size    = "bin50",
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
        # binGefInput's bin_size slot is just the numeric key under
        # geneExp/, e.g. "50" for bin50. The StereoSeqReader-side
        # convention is "bin50"/"bin100"/... -- strip the "bin" prefix.
        bin_key <- sub("^bin", "", as.character(bin_size))
        inp <- binGefInput(path,
                            bin_size    = bin_key,
                            gene_column = gene_column)
    }

    pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                       verbose = verbose, ...)

    if (output == "store") return(pe)

    methods::new("exprObj",
        name       = "raw",
        exprMat    = pe,
        spat_unit  = spat_unit,
        feat_type  = "rna",
        provenance = spat_unit
    )
}
