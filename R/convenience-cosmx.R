# CosMx ingest pipeline for `gDirSource`-managed projects.
#
# Disk-backed counterpart to Giotto's `CosmxReader`. Currently routes only
# the expression matrix through GiottoDisk (`parquetExprStore` written
# into the project vault via `sourceWrite`). Transcripts / polys / images
# / cellmeta remain on the inherited in-mem closures from `CosmxReader`;
# they can be ported following the same pattern when needed.



# CLASS ####



setClass(
    "CosMxDiskReader",
    contains = "CosmxReader",
    slots = list(
        backend = "ANY"
    ),
    prototype = list(
        backend = NULL
    )
)

# * init ####
setMethod(
    "initialize", signature("CosMxDiskReader"),
    function(.Object, ..., backend) {
        obj <- callNextMethod(.Object, ...)

        if (!missing(backend)) {
            # Match createGiottoObject(backend = ...): character -> coerce
            # to gDirSource (path).
            if (is.character(backend)) {
                backend <- gDirSource(path = backend)
            }
            checkmate::assert_class(backend, "gsource")
            obj@backend <- backend
        }
        if (is.null(obj@backend)) {
            stop("[CosMxDiskReader] `backend` is required", call. = FALSE)
        }

        # Mirror @paths into the init frame so subclass closures can
        # reference path names directly via default-arg expressions
        # (same convention as XeniumDiskReader). Each name resolves to a
        # bare character string at call time.
        list2env(obj@paths, envir = environment())
        gsrc  <- obj@backend
        slide <- obj@slide
        fovs_ <- obj@fovs

        # expression (disk override)
        ex_fun <- function(
            path = expr_path,
            feat_type = c("rna", "negprobes"),
            split_keyword = list("NegPrb"),
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            .cosmx_expression_disk(
                path = path,
                gsource = gsrc,
                slide = slide,
                fovs = fovs_ %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_expression <- ex_fun

        # create_gobject (disk variant). Mirrors parent's gobject_fun but
        # initializes the giotto object with backend = gsrc so any further
        # artifacts are vault-resident.
        gobject_fun <- function(
            transcript_path = tx_path,
            expression_path = expr_path,
            metadata_path = meta_path,
            cell_labels_dir = cell_labels_dir,
            composite_img_dir = composite_img_dir,
            overlay_img_dir = overlay_img_dir,
            feat_type = c("rna", "negprobes"),
            split_keyword = list("NegPrb"),
            load_images = list(
                composite = "composite",
                overlay = "overlay"
            ),
            image_negative_y = NULL,
            load_expression = FALSE,
            load_cellmeta = TRUE,
            load_transcripts = TRUE,
            instructions = NULL,
            cores = GiottoUtils::determine_cores(),
            verbose = NULL) {
            load_expression  <- as.logical(load_expression)
            load_cellmeta    <- as.logical(load_cellmeta)
            load_transcripts <- as.logical(load_transcripts)

            funs <- obj@calls

            # init gobject with disk backend
            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            # transcripts (inherited; in-mem)
            if (isTRUE(load_transcripts)) {
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    cores = cores,
                    verbose = verbose
                )
                for (tx in tx_list) {
                    g <- GiottoClass::setGiotto(g, tx)
                }
            }

            # polys (inherited; in-mem)
            poly_args <- list(path = cell_labels_dir, verbose = FALSE)
            if (!is.null(image_negative_y)) {
                poly_args$shift_vertical_step <- if (isTRUE(image_negative_y)) FALSE else 1
            }
            polys <- do.call(funs$load_polys, poly_args)
            g <- GiottoClass::setGiotto(g, polys, verbose = verbose)

            # images (inherited; in-mem)
            if (!is.null(load_images)) {
                load_images[load_images == "composite"] <- composite_img_dir
                load_images[load_images == "overlay"]   <- overlay_img_dir
                imglist  <- list()
                dirnames <- names(load_images)
                for (imdir_i in seq_along(load_images)) {
                    img_args <- list(
                        path = load_images[[imdir_i]],
                        img_type = dirnames[[imdir_i]],
                        verbose = verbose
                    )
                    if (!is.null(image_negative_y)) {
                        img_args$negative_y <- image_negative_y
                    }
                    dir_imgs <- do.call(funs$load_images, img_args)
                    imglist <- c(imglist, dir_imgs)
                }
                g <- GiottoClass::addGiottoLargeImage(
                    g, largeImages = imglist, verbose = FALSE
                )
            }

            # expression & meta cell-ID intersection (mirrors parent)
            allowed_ids <- GiottoClass::spatIDs(polys)

            # expression (disk; overridden closure)
            if (load_expression) {
                exlist <- funs$load_expression(
                    path = expression_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    verbose = verbose
                )
                for (ex in exlist) {
                    # Subset on cell IDs via lazy parquetExprStore slicing
                    bool <- colnames(ex[]) %in% allowed_ids
                    ex[] <- ex[][, bool]
                    g <- GiottoClass::setGiotto(g, ex, verbose = verbose)
                }
            }

            # cellmeta (inherited; in-mem)
            if (load_cellmeta) {
                cx <- funs$load_cellmeta(
                    path = metadata_path, cores = cores, verbose = verbose
                )
                cx[] <- cx[][cell_ID %in% allowed_ids, ]
                g <- GiottoClass::setGiotto(g, cx, verbose = verbose)
            }

            # spatlocs from polygon centroids (mirrors parent)
            g <- GiottoClass::addSpatialCentroidLocations(g, verbose = FALSE)

            # add fovs metadata column
            g$fov <- gsub("^c_\\d+_(\\d+)_\\d+$", "\\1",
                          GiottoClass::pDataDT(g)$cell_ID)

            g
        }
        obj@calls$create_gobject <- gobject_fun

        obj
    }
)



# CREATE READER ####

#' @title Import a NanoString CosMx assay (disk-backed)
#' @name importCosMxDisk
#' @description
#' Disk-backed counterpart to [Giotto::importCosMx()]. Produces a
#' `CosMxDiskReader` whose `load_expression()` call writes a
#' `parquetExprStore` into a `gDirSource`-managed project vault.
#' Transcripts / polys / images / cellmeta remain in-memory via the
#' inherited `CosmxReader` closures.
#' @param cosmx_dir CosMx output directory
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#'   Naming matches [GiottoClass::createGiottoObject()]'s `backend` param.
#' @param slide,fovs,version,micron,px2um,poly_pref passed through to the
#'   parent `CosmxReader` initializer.
#' @returns `CosMxDiskReader` object
#' @seealso [Giotto::importCosMx()] for the in-memory variant
#' @export
importCosMxDisk <- function(cosmx_dir = NULL,
                              backend,
                              slide = 1,
                              fovs = NULL,
                              version = "default",
                              micron = FALSE,
                              px2um = 0.12028,
                              poly_pref = c("mask", "csv")) {
    if (missing(backend)) {
        stop("[importCosMxDisk] `backend` is required", call. = FALSE)
    }
    a <- list(
        Class = "CosMxDiskReader",
        backend = backend,
        slide = slide,
        version = version,
        micron = micron,
        px2um = px2um,
        poly_pref = match.arg(poly_pref)
    )
    if (!is.null(cosmx_dir)) a$cosmx_dir <- cosmx_dir
    if (!is.null(fovs)) a$fovs <- fovs
    do.call(new, args = a)
}



# MODULAR ####


## expression ####

# Disk-backed CosMx expression ingestion. Wraps the wide-format
# exprMat_file CSV in a csvWideInput marker that:
#   - drops cell_ID == 0 (background) and (optionally) restricts FOVs
#   - skips the non-feature `fov` column
# Routes through sourceWrite(gsource, inp, store_type = "parquetExpr").
# CosMx-specific cell ID reconstruction (`c_<slide>_<fov>_<cell_ID>`)
# is applied post-write; split_keyword feat_type splits are applied
# lazily via pe[i, ] gene-row slicing -- no parquet rewrite.
.cosmx_expression_disk <- function(
    path,
    gsource,
    slide = 1,
    fovs = NULL,
    feat_type = c("rna", "negprobes"),
    split_keyword = list("NegPrb"),
    output = c("exprObj", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path) || length(path) == 0L || !file.exists(path)) {
        stop("[cosmx_expression_disk] no exprMat_file path provided",
             call. = FALSE)
    }
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("exprObj", "store"))

    GiottoUtils::vmsg("[cosmx_expression_disk] streaming CSV ->",
                       "parquetExprStore", .v = verbose)

    .fovs <- if (!is.null(fovs)) as.integer(fovs) else NULL
    row_filter <- function(chunk) {
        keep <- chunk[["cell_ID"]] != 0L
        if (!is.null(.fovs)) {
            keep <- keep & chunk[["fov"]] %in% .fovs
        }
        keep
    }

    inp <- csvWideInput(
        csv_path        = path,
        cell_id_col     = "cell_ID",
        skip_cols       = "fov",
        row_filter_fun  = row_filter
    )

    pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                       verbose = verbose, ...)

    # Reconstruct globally-unique cell IDs `c_<slide>_<fov>_<cell_ID>`
    # by re-reading just the (fov, cell_ID) columns from the source CSV.
    # Mirrors dev_stream's .cosmx_expression_parquet logic.
    cell_ID <- NULL  # NSE binding
    id_dt <- data.table::fread(path, select = c("fov", "cell_ID"))
    id_dt <- id_dt[cell_ID != 0L, ]
    if (!is.null(.fovs)) id_dt <- id_dt[fov %in% .fovs, ]
    if (nrow(id_dt) != length(pe@cell_ids)) {
        stop("[cosmx_expression_disk] cell-row mismatch: filtered CSV ",
             "rows = ", nrow(id_dt), ", parquetExprStore cells = ",
             length(pe@cell_ids), ". Filter logic disagrees with the ",
             "streamed write.", call. = FALSE)
    }
    pe@cell_ids <- sprintf("c_%d_%d_%d", slide, id_dt$fov, id_dt$cell_ID)

    # split_keyword feat_type splits via lazy gene-row slicing.
    feat_ids  <- pe@feat_ids
    expr_list <- vector("list", length(feat_type))
    names(expr_list) <- feat_type
    if (length(split_keyword) == 0L) {
        expr_list <- list(pe)
        names(expr_list) <- feat_type[[1L]]
    } else {
        remaining <- rep(TRUE, length(feat_ids))
        for (key_i in seq_along(split_keyword)) {
            bool <- grepl(pattern = split_keyword[[key_i]], x = feat_ids) &
                    remaining
            if (!any(bool)) {
                expr_list[[key_i + 1L]] <- NULL
                next
            }
            expr_list[[key_i + 1L]] <- pe[which(bool), , drop = FALSE]
            remaining <- remaining & !bool
        }
        expr_list[[1L]] <- pe[which(remaining), , drop = FALSE]
        expr_list <- Filter(Negate(is.null), expr_list)
    }

    if (output == "store") return(expr_list)

    # Wrap each split in an exprObj. Use new() to bypass
    # .evaluate_expr_matrix when needed (older GiottoClass installs).
    lapply(seq_along(expr_list), function(i) {
        methods::new("exprObj",
            name       = "raw",
            exprMat    = expr_list[[i]],
            spat_unit  = "cell",
            feat_type  = names(expr_list)[[i]],
            provenance = "cell"
        )
    })
}
