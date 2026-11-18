#' @include convenience-xenium.R reader-shared.R
NULL

# Standard Visium ingest for `gDirSource`-managed projects.
#
# Disk-backed counterpart to Giotto's `VisiumReader`. Expression goes through
# the shared 10x path (`.tenx_expression_disk`) into a `parquetExprStore`.
# Spot locations, spot metadata, scalefactors, the micron scale, the image and
# the circular spot polygons are the inherited in-memory loaders: a Visium run has a few
# thousand spots. Polygons are written to the vault by the backed object's
# setter when attached.



# CLASS ####



setClass(
    "VisiumDiskReader",
    contains = "VisiumReader",
    slots = list(
        backend = "ANY"
    ),
    prototype = list(
        backend = NULL
    )
)

# * init ####
setMethod(
    "initialize", signature("VisiumDiskReader"),
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
            stop("[VisiumDiskReader] `backend` is required", call. = FALSE)
        }

        list2env(obj@paths, envir = environment())
        gsrc <- obj@backend
        feat_id <- obj@feature_id_type

        ## expression (disk) ####
        ex_fun <- function(
            path = expr_path,
            feature_id_type = feat_id,
            output = c("exprObj", "store"),
            verbose = NULL
        ) {
            .visium_expression_disk(
                path = path,
                gsource = gsrc,
                feature_id_type = feature_id_type,
                output = output,
                verbose = verbose
            )
        }
        obj@calls$load_expression <- ex_fun

        ## giotto object ####
        # Same arguments and load order as the parent's create_gobject(), plus
        # a backed `gobject` to add into.
        gobject_fun <- function(
            load_expression = TRUE,
            load_spatlocs = TRUE,
            load_metadata = TRUE,
            load_image = TRUE,
            load_polygons = TRUE,
            expression_path = expr_path,
            tissue_positions_path = positions_path,
            scalefactors_path = json_path,
            image_path = img_path,
            feature_id_type = feat_id,
            gobject = NULL,
            instructions = NULL,
            verbose = NULL
        ) {
            funs <- obj@calls
            g <- .backed_gobject_init(gobject, gsrc, instructions)

            barcodes <- NULL
            feat_types <- "rna"
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    feature_id_type = feature_id_type,
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, ex, verbose = FALSE)
                barcodes <- colnames(ex[[1L]][])
                feat_types <- vapply(ex, GiottoClass::featType, character(1L))
            }

            sl <- funs$load_tissue_position(
                path = tissue_positions_path, barcodes = barcodes
            )
            if (load_spatlocs) g <- GiottoClass::setGiotto(g, sl, verbose = FALSE)

            if (load_metadata) {
                cx <- funs$load_cellmeta(
                    path = tissue_positions_path,
                    barcodes = barcodes,
                    feat_type = feat_types
                )
                g <- GiottoClass::setGiotto(g, cx, verbose = FALSE)
            }

            json <- funs$load_scalefactors(path = scalefactors_path)

            if (load_image && !is.null(image_path)) {
                img <- funs$load_image(
                    path = image_path, scalefactors = json, verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, img, verbose = FALSE)
            }

            ms <- funs$load_micron_scale(
                path = tissue_positions_path, scalefactors = json
            )
            if (load_polygons) {
                poly <- funs$load_polygon(spatlocs = sl, micron_scale = ms)
                g <- GiottoClass::setPolygonInfo(
                    gobject = g,
                    x = poly,
                    centroids_to_spatlocs = FALSE,
                    verbose = FALSE,
                    initialize = TRUE
                )
            }
            GiottoClass::instructions(g, "micron_scale") <- ms

            g
        }
        obj@calls$create_gobject <- gobject_fun

        obj
    }
)



# CREATE READER ####

#' @title Import a 10x Visium assay (disk-backed)
#' @name importVisiumDisk
#' @description
#' Disk-backed counterpart to [Giotto::importVisium()]. Produces a
#' `VisiumDiskReader` whose expression loader streams the 10x feature-barcode
#' matrix (mtx directory or `.h5`) into a `parquetExprStore` in the
#' `gDirSource`-managed project vault. Spot locations, spot metadata, the
#' image and the circular spot polygons come from the inherited
#' `VisiumReader` loaders, so barcodes, coordinates and feat types match the
#' in-memory reader. `create_gobject()` takes the in-memory reader's
#' arguments, so `Giotto::createGiottoVisiumObject(backend =)` uses it, and
#' additionally a backed `gobject` to add into.
#' @param visium_dir Space Ranger `outs` directory for a standard Visium run
#' @param backend a `gsource` (typically `gDirSource`) project backend, or a
#'   directory path.
#' @param expression_source,feature_id_type passed through to the parent
#'   `VisiumReader`; see [Giotto::importVisium()].
#' @returns `VisiumDiskReader` object
#' @seealso [Giotto::importVisium()] for the in-memory variant
#' @examples
#' \dontrun{
#' g <- importVisiumDisk("path/to/outs", backend = "path/to/project")$
#'     create_gobject()
#' }
#' @export
importVisiumDisk <- function(
    visium_dir = NULL,
    backend,
    expression_source = c("raw", "filtered"),
    feature_id_type = c("symbol", "ensembl")
) {
    if (missing(backend)) {
        stop("[importVisiumDisk] `backend` is required", call. = FALSE)
    }
    # Same argument handling as Giotto::importVisium(), so adding
    # `backend =` to a call changes nothing about what is read.
    a <- list(
        Class = "VisiumDiskReader",
        backend = backend,
        expression_source = match.arg(expression_source, c("raw", "filtered")),
        feature_id_type = match.arg(feature_id_type, c("symbol", "ensembl"))
    )
    if (!is.null(visium_dir)) a$visium_dir <- visium_dir
    do.call(methods::new, args = a)
}



# MODULAR ####


## expression ####

# Shared 10x ingest, then the feat types the in-memory reader keeps:
# "Gene Expression" as rna and "Antibody Capture" as protein.
.visium_expression_disk <- function(
    path,
    gsource,
    feature_id_type = c("symbol", "ensembl"),
    output = c("exprObj", "store"),
    verbose = NULL
) {
    if (missing(path) || is.null(path) || !file.exists(path)) {
        stop("[visium_expression_disk] no expression path found",
             call. = FALSE)
    }
    feature_id_type <- match.arg(feature_id_type, c("symbol", "ensembl"))
    output <- match.arg(output)
    GiottoUtils::vmsg(.v = verbose, "[visium_expression_disk]", path)

    stores <- .tenx_expression_disk(
        path = path,
        gsource = gsource,
        gene_ids = if (feature_id_type == "symbol") "symbols" else "ensembl",
        spat_unit = "cell",
        output = "store",
        verbose = verbose
    )
    stores <- list(rna = stores[["rna"]], protein = stores[["Antibody_Capture"]])
    stores <- stores[!vapply(stores, is.null, logical(1L))]
    if (output == "store") return(stores)

    lapply(names(stores), function(ft) {
        methods::new("exprObj",
            name       = "raw",
            exprMat    = stores[[ft]],
            spat_unit  = "cell",
            feat_type  = ft,
            provenance = "cell"
        )
    })
}
