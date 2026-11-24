# CosMx ingest pipeline for `gDirSource`-managed projects.
#
# Disk-backed counterpart to Giotto's `CosmxReader`. Expression, protein,
# polygons and transcripts are written to the project vault through
# `sourceWrite()`; images and cell metadata use the inherited in-mem closures.



# CLASS ####



setClass(
    "CosMxDiskReader",
    contains = "CosmxReader",
    slots = list(
        backend = "ANY",
        protein_dir = "character"
    ),
    prototype = list(
        backend = NULL,
        protein_dir = character(0L)
    )
)

# * init ####
setMethod(
    "initialize", signature("CosMxDiskReader"),
    function(.Object, ..., backend, protein_dir) {
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

        if (!missing(protein_dir)) {
            if (length(protein_dir)) {
                checkmate::assert_directory_exists(protein_dir)
            }
            obj@protein_dir <- as.character(protein_dir)
        }
        # CosMx protein parsing lives on the parent reader (@helpers), so the
        # two readers share one copy; only the storage step is ours
        if (!methods::.hasSlot(obj, "helpers") || !length(obj@helpers)) {
            stop("[CosMxDiskReader] the CosmxReader class in use has no ",
                 "protein helpers: GiottoDisk was installed against an older ",
                 "Giotto. Reinstall GiottoDisk.", call. = FALSE)
        }
        h <- obj@helpers
        prot_paths <- h$protein_paths(obj@protein_dir)
        obj@paths[names(prot_paths)] <- prot_paths

        # Mirror @paths into the init frame so subclass closures can
        # reference path names directly via default-arg expressions
        # (same convention as XeniumDiskReader). Each name resolves to a
        # bare character string at call time.
        list2env(obj@paths, envir = environment())
        # Dotted copies: a closure default written as `meta_path = meta_path`
        # resolves against the closure's OWN formals first and self-references,
        # so the closures below bind these instead.
        .poly_path <- obj@paths$poly_path
        .meta_path <- obj@paths$meta_path
        .expr_path <- obj@paths$expr_path
        .tx_path <- obj@paths$tx_path
        .composite_dir <- obj@paths$composite_img_dir
        .overlay_dir <- obj@paths$overlay_img_dir
        .prot_expr_path <- obj@paths$protein_expr_path
        .prot_meta_path <- obj@paths$protein_meta_path
        .prot_plex_path <- obj@paths$protein_plex_path
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

        # protein values are decimals: never the integer cosmxscan scanner
        prot_fun <- function(
            path = .prot_expr_path,
            feat_type = c("protein", "protein_neg"),
            split_keyword = NULL,
            drop_cols = "^Channel-",
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            split_keyword <- h$protein_spec(path, feat_type, split_keyword,
                                            plex_path = .prot_plex_path)
            .cosmx_expression_disk(
                path = path,
                gsource = gsrc,
                slide = slide,
                fovs = fovs_ %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                drop_cols = drop_cols,
                scanner = "csv",
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_protein <- prot_fun

        # polygons (disk override)
        poly_fun <- function(
            path = .poly_path,
            name = "cell",
            part_col = NULL,
            flip_vertical = FALSE,
            calc_centroids = TRUE,
            keep_ids = NULL,
            expr_ids = keep_ids,
            output = c("giottoPolygon", "store"),
            verbose = NULL,
            ...
        ) {
            .cosmx_poly_disk(
                path = path,
                gsource = gsrc,
                name = name,
                part_col = part_col,
                flip_vertical = flip_vertical,
                fovs = fovs_ %none% NULL,
                keep_ids = keep_ids,
                expr_ids = expr_ids,
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_polys <- poly_fun

        # transcripts (disk override); `col_classes`, `micron` and `px2um`
        # have no disk counterpart, as for the polygons
        tx_fun <- function(
            path = .tx_path,
            feat_type = c("rna", "negprobes"),
            split_keyword = list("NegPrb"),
            dropcols = c(),
            cores = GiottoUtils::determine_cores(),
            output = c("giottoPoints", "store"),
            verbose = NULL,
            ...
        ) {
            .cosmx_transcript_disk(
                path = path,
                gsource = gsrc,
                version = obj@version,
                fovs = fovs_ %none% NULL,
                feat_type = feat_type,
                split_keyword = split_keyword,
                dropcols = dropcols,
                cores = cores,
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_transcripts <- tx_fun

        # create_gobject (disk variant). Mirrors parent's gobject_fun but
        # initializes the giotto object with backend = gsrc; defaults differ
        # where the disk path reads differently (matrix on, transcripts off).
        gobject_fun <- function(
            transcript_path = .tx_path,
            cell_labels_dir = .poly_path,
            expression_path = .expr_path,
            metadata_path = .meta_path,
            feat_type = c("rna", "negprobes", "falsecode"),
            split_keyword = list("^Negative", "^SystemControl"),
            load_images = list(
                composite = "composite",
                overlay = "overlay"
            ),
            image_negative_y = NULL,
            load_expression = TRUE,
            load_cellmeta = TRUE,
            load_transcripts = FALSE,
            instructions = NULL,
            cores = GiottoUtils::determine_cores(),
            verbose = NULL,
            # protein (multiomic exports)
            load_protein = !is.null(protein_expression_path),
            protein_expression_path = .prot_expr_path,
            protein_metadata_path = .prot_meta_path,
            protein_feat_type = c("protein", "protein_neg"),
            protein_split_keyword = NULL,
            protein_drop = "^Channel-",
            protein_meta_cols = c("nCount_RNA", "nFeature_RNA",
                                  "nCount_negprobes", "nFeature_negprobes",
                                  "Panel", "assay_type"),
            keep_cells = "intersect",
            # disk only
            polygon_path = cell_labels_dir,
            poly_read_fun = NULL,
            meta_cols = c("fov", "Area.um2", "nCount_RNA", "nFeature_RNA"),
            load_polygons = TRUE
        ) {
            keep_cells <- match.arg(keep_cells, choices = "intersect")
            load_protein <- isTRUE(load_protein)
            if (load_protein && (!length(protein_expression_path) ||
                                 !file.exists(protein_expression_path))) {
                stop("[CosMxDiskReader] `load_protein = TRUE` but no protein ",
                     "exprMat_file was found or supplied: set `protein_dir` ",
                     "or `protein_expression_path`.", call. = FALSE)
            }
            # polygons are read from the vertex CSV; the parent's mask
            # (CellLabels) route has no disk counterpart
            if (isTRUE(load_polygons) && length(polygon_path) &&
                dir.exists(polygon_path)) {
                stop("[CosMxDiskReader] polygons are read from the ",
                     "*-polygons.csv file; mask images (CellLabels) are not ",
                     "supported with a backend. Pass the CSV as ",
                     "`cell_labels_dir`/`polygon_path`.", call. = FALSE)
            }
            if (!is.null(load_images)) {
                checkmate::assert_list(load_images)
                if (is.null(names(load_images))) {
                    stop("[CosMxDiskReader] 'load_images' must be a named ",
                         "list of image directories", call. = FALSE)
                }
            }

            funs <- obj@calls

            # init gobject with disk backend
            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            # expression before polygons: a cell with a polygon but no
            # expression fails the gobject init check, so its polygon is dropped
            exlist <- NULL
            if (isTRUE(load_expression)) {
                exlist <- funs$load_expression(
                    path = expression_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    verbose = verbose
                )
            }
            plist <- NULL
            if (load_protein) {
                plist <- funs$load_protein(
                    path = protein_expression_path,
                    feat_type = protein_feat_type,
                    split_keyword = protein_split_keyword,
                    drop_cols = protein_drop,
                    verbose = verbose
                )
            }
            rna_ids <- if (length(exlist)) {
                colnames(exlist[[1L]][])
            } else if (load_protein && isTRUE(load_transcripts)) {
                .cosmx_read_cellmeta(metadata_path, character(0L),
                    fovs = fovs_ %none% NULL, cores = cores,
                    verbose = verbose)$cell_ID
            }
            prot_ids <- if (length(plist)) colnames(plist[[1L]][])
            keep <- h$match_cells(list(rna = rna_ids, protein = prot_ids),
                                  keep_cells = keep_cells)

            # polygons (disk; overridden closure)
            allowed_ids <- NULL
            poly_ids <- keep
            if (isTRUE(load_polygons)) {
                polys <- funs$load_polys(
                    path = polygon_path, name = "cell",
                    read_fun = poly_read_fun,
                    keep_ids = keep,
                    expr_ids = unique(c(rna_ids, prot_ids)),
                    verbose = verbose
                )
                allowed_ids <- GiottoClass::spatIDs(polys)
                g <- GiottoClass::setGiotto(
                    g, polys, centroids_to_spatlocs = TRUE, verbose = FALSE
                )
                if (!is.null(keep) && !any(keep %in% allowed_ids)) {
                    stop("[CosMxDiskReader] expression and polygon ",
                         "cell IDs do not intersect; both should be ",
                         "'c_<slide>_<fov>_<cell>'.", call. = FALSE)
                }
                poly_ids <- allowed_ids
            }
            h$report_cells(keep, poly_ids)
            # also drops the "dropped" attribute, so `keep` compares to ids
            if (!is.null(keep)) keep <- keep[keep %in% poly_ids]

            # images (inherited in-mem closure); folders not found are skipped,
            # since flat-file exports ship none
            if (!is.null(load_images)) {
                found <- list(composite = .composite_dir,
                              overlay = .overlay_dir)
                for (i in seq_along(load_images)) {
                    im <- load_images[[i]]
                    if (is.character(im) && length(im) == 1L &&
                        im %in% names(found)) {
                        load_images[i] <- list(found[[im]])
                    }
                }
                absent <- vapply(load_images, function(p) {
                    !length(p) || !dir.exists(p[[1L]])
                }, logical(1L))
                if (any(absent)) {
                    warning(GiottoUtils::wrap_txt(sprintf(
                        "[CosMxDiskReader] no %s image folder found;
                        skipping. Pass `load_images` paths to load images
                        from elsewhere.",
                        paste(names(load_images)[absent], collapse = ", ")
                    )), call. = FALSE)
                    load_images <- load_images[!absent]
                }
                imglist <- list()
                for (i in seq_along(load_images)) {
                    img_args <- list(
                        path = load_images[[i]],
                        img_type = names(load_images)[[i]],
                        verbose = verbose
                    )
                    if (!is.null(image_negative_y)) {
                        img_args$negative_y <- image_negative_y
                    }
                    imglist <- c(imglist, do.call(funs$load_images, img_args))
                }
                if (length(imglist)) {
                    g <- GiottoClass::setGiotto(g, imglist, verbose = FALSE)
                }
            }

            # GiottoClass does not check the cells of later expression sets
            for (ex in c(exlist, plist)) {
                if (!is.null(keep)) {
                    ids <- colnames(ex[])
                    if (!identical(ids, keep)) {
                        ex[] <- ex[][, match(keep, ids)]
                    }
                    if (!identical(colnames(ex[]), keep)) {
                        stop("[CosMxDiskReader] internal: cells of feat_type '",
                             ex@feat_type, "' could not be matched to the ",
                             "selected cells.", call. = FALSE)
                    }
                }
                g <- GiottoClass::setGiotto(g, ex, verbose = FALSE)
            }

            # transcripts (disk): off by default, the largest object here
            if (isTRUE(load_transcripts)) {
                if (is.null(transcript_path) ||
                    !file.exists(transcript_path)) {
                    stop("[CosMxDiskReader] `load_transcripts = TRUE` but no ",
                         "tx_file was found or supplied.", call. = FALSE)
                }
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    cores = cores,
                    verbose = verbose
                )
                g <- GiottoClass::setGiotto(g, tx_list, verbose = FALSE)
            }

            # cellmeta, read from the vendor CSV rather than through the
            # parent closure: `cell` is the vendor's GLOBAL id and the only
            # column that joins to the polygons, while `cell_ID` in that same
            # file is the FOV-local integer. Read once, set on both feat_types.
            rna_ft <- if (length(exlist)) feat_type[[1L]]
            prot_ft <- if (length(plist)) protein_feat_type[[1L]]
            meta_ids <- keep %null% allowed_ids
            fovs_use <- fovs_ %none% NULL
            if (isTRUE(load_cellmeta)) {
                cx <- NULL
                if (length(exlist) || !length(plist)) {
                    cx <- .cosmx_read_cellmeta(metadata_path, meta_cols,
                        ids = meta_ids, fovs = fovs_use, cores = cores,
                        verbose = verbose)
                    if (!is.null(cx)) {
                        g <- GiottoClass::addCellMetadata(
                            g, feat_type = rna_ft, new_metadata = cx,
                            by_column = TRUE, column_cell_ID = "cell_ID"
                        )
                    }
                }
                if (length(plist)) {
                    shared <- setdiff(meta_cols, protein_meta_cols)
                    pm <- .cosmx_read_cellmeta(protein_metadata_path,
                        cols = if (is.null(cx)) {
                            c(shared, protein_meta_cols)
                        } else {
                            protein_meta_cols
                        },
                        ids = meta_ids, fovs = fovs_use, cores = cores,
                        verbose = verbose)
                    if (!is.null(pm)) {
                        own <- intersect(protein_meta_cols, names(pm))
                        data.table::setnames(pm, own,
                            h$protein_meta_names(own, protein_feat_type))
                    }
                    if (!is.null(cx)) {
                        sx <- cx[, intersect(c("cell_ID", shared), names(cx)),
                                 with = FALSE]
                        pm <- if (is.null(pm)) sx else
                            merge(sx, pm, by = "cell_ID", all.x = TRUE,
                                  sort = FALSE)
                    }
                    if (!is.null(pm) && ncol(pm) > 1L) {
                        g <- GiottoClass::addCellMetadata(
                            g, feat_type = prot_ft, new_metadata = pm,
                            by_column = TRUE, column_cell_ID = "cell_ID"
                        )
                    }
                }
            }

            # add fovs metadata column. Never `g$fov <- ...` here -- on this
            # object that is a silent no-op, no error and no column.
            fts <- c(rna_ft, prot_ft)
            if (!length(fts)) {
                g <- .cosmx_add_fov(g)
            }
            for (ft in fts) {
                g <- .cosmx_add_fov(g, feat_type = ft)
            }

            g
        }
        obj@calls$create_gobject <- gobject_fun

        obj
    }
)

# * access ####
#' @export
setMethod("$", signature("CosMxDiskReader"), function(x, name) {
    if (identical(name, "protein_dir")) return(x@protein_dir)
    callNextMethod()
})

#' @export
setMethod("$<-", signature("CosMxDiskReader"), function(x, name, value) {
    if (identical(name, "protein_dir")) {
        return(initialize(x, protein_dir = value))
    }
    callNextMethod()
})



# CREATE READER ####

#' @title Import a NanoString CosMx assay (disk-backed)
#' @name importCosMxDisk
#' @description
#' Disk-backed counterpart to [Giotto::importCosMx()]. Produces a
#' `CosMxDiskReader` whose `load_expression()`, `load_polys()` and
#' `load_transcripts()` calls write to a `gDirSource`-managed project vault:
#' a `parquetExprStore` for the expression matrix and `parquetGeomTile` stores
#' for the polygons and transcripts. Images and cell metadata remain in-memory
#' via the inherited `CosmxReader` closures.
#'
#' The transcript flatfile is converted to parquet once, in a cache under
#' [getArtifactDumpDir()] (first import converts, later imports reuse).
#' Polygons of cells that have no expression (cells removed by the vendor's
#' QC) are not loaded, with a warning.
#'
#' `fovs` applies to every input while it is read: the expression and protein
#' matrices (both the `cosmxscan` and the CSV route), the polygons, the cell
#' metadata and the transcripts. A gzipped file is still decompressed in full,
#' since gzip cannot seek.
#' @section Protein (multiomic exports):
#' The protein export of a CosMx multiomic run is a second folder whose files
#' have the same names as the RNA export, given as `protein_dir`. Its
#' `exprMat_file`, `metadata_file` and `plex-*.txt` are found by pattern.
#' `$load_protein()` writes the protein matrix to a `parquetExprStore` as raw
#' mean intensities (decimals, so never through `cosmxscan`), without the
#' `Channel-*` morphology columns (they repeat the metadata `Mean.*`), split
#' into the targets (`"protein"`) and the isotype controls (`"protein_neg"`:
#' `CodeClass` `"Negative"` in the plex file, else names matching `"IgG"`).
#'
#' `$create_gobject()` then also takes:
#' * `load_protein` (default: `TRUE` when a protein matrix was found),
#'   `protein_expression_path`, `protein_metadata_path`;
#' * `protein_feat_type` (`c("protein", "protein_neg")`; one name keeps the
#'   isotypes with the targets), `protein_split_keyword` (`NULL`: from the
#'   plex file), `protein_drop` (`"^Channel-"`; `NULL` keeps those columns);
#' * `protein_meta_cols`: protein metadata columns to add to the protein
#'   feat_type, renamed because the vendor file reuses the RNA names for
#'   protein totals: `nCount_RNA` -> `nCount_protein` (the sum of the targets
#'   **and** the `Channel-*` intensities, so not the sum of the stored
#'   matrix), `nFeature_RNA` -> `nFeature_protein`, `nCount_negprobes` ->
#'   `nCount_protein_neg`, `nFeature_negprobes` -> `nFeature_protein_neg`,
#'   any other name gets a `_protein` suffix (`Panel_protein`,
#'   `assay_type_protein`). The `meta_cols` not in `protein_meta_cols`
#'   (`fov`, `Area.um2`) describe the cell and are set on both feat_types;
#' * `keep_cells`: `"intersect"` (the only mode) keeps the cells in the RNA
#'   export, protein expression and a polygon, warning with the number
#'   dropped for each reason. The RNA cells are those of the RNA matrix, or,
#'   when RNA is loaded only as transcripts, those of the RNA cell metadata;
#'   with neither, every protein cell with a polygon is kept. Every feature
#'   type is put in the same cell order (RNA file order). With
#'   `load_expression = FALSE` the object holds no RNA expression.
#' @param cosmx_dir CosMx output directory
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#'   Naming matches [GiottoClass::createGiottoObject()]'s `backend` param.
#' @param slide,fovs,version,micron,px2um,poly_pref passed through to the
#'   parent `CosmxReader` initializer. `poly_pref` defaults to `"csv"`
#'   because the disk polygon loader reads the vertex CSV, not masks.
#' @param protein_dir protein export folder of a multiomic run (optional).
#'   Kept separate from `cosmx_dir` because the file names are the same.
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
                              poly_pref = c("csv", "mask"),
                              protein_dir = NULL) {
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
    if (!is.null(protein_dir)) a$protein_dir <- protein_dir
    do.call(new, args = a)
}


# MODULAR ####


## expression ####

# Disk-backed CosMx expression ingestion of a wide-format exprMat_file:
#   - drops cell_ID == 0 (background) and (optionally) restricts FOVs
#   - skips the non-feature `fov` column and columns matching `drop_cols`
# Routes through sourceWrite(gsource, inp, store_type = "parquetExpr").
# split_keyword feat_type splits are applied lazily via pe[i, ] gene-row
# slicing -- no parquet rewrite.
.cosmx_expression_disk <- function(
    path,
    gsource,
    slide = 1,
    fovs = NULL,
    feat_type = c("rna", "negprobes"),
    split_keyword = list("NegPrb"),
    drop_cols = NULL,
    scanner = c("auto", "cosmxscan", "csv"),
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
    scanner <- match.arg(scanner, choices = c("auto", "cosmxscan", "csv"))

    .fovs <- if (!is.null(fovs)) as.integer(fovs) else NULL

    drop <- character(0L)
    if (length(drop_cols)) {
        hdr <- .cosmx_csv_header(path)
        drop <- grep(paste(drop_cols, collapse = "|"), hdr, value = TRUE)
        drop <- setdiff(drop, c("fov", "cell_ID"))
    }

    # cosmxscan emits only nonzeros: ~100 s against ~974 s for csvWideInput on
    # a 20,378-feature panel. It reads integer counts only and cannot skip
    # interior columns, so it is used only where both hold.
    if (scanner == "auto") {
        use_scanner <- requireNamespace("cosmxscan", quietly = TRUE) &&
            grepl("[.]gz$", path, ignore.case = TRUE) &&
            !length(drop) &&
            .cosmx_values_are_integer(path)
        scanner <- if (use_scanner) "cosmxscan" else "csv"
    }
    if (scanner == "cosmxscan" && length(drop)) {
        stop("[cosmx_expression_disk] cosmxscan cannot skip the columns ",
             "matched by `drop_cols`; use scanner = \"csv\"", call. = FALSE)
    }
    na_ids <- character(0L)
    row_filter <- function(chunk) {
        keep <- chunk[["cell_ID"]] != 0L
        if (!is.null(.fovs)) {
            keep <- keep & chunk[["fov"]] %in% .fovs
        }
        feats <- setdiff(names(chunk), c("fov", "cell_ID", drop))
        na <- keep & !stats::complete.cases(chunk[, feats, with = FALSE])
        if (any(na)) {
            na_ids <<- c(na_ids, sprintf("c_%d_%d_%d", slide,
                                         chunk[["fov"]][na],
                                         chunk[["cell_ID"]][na]))
        }
        keep & !na
    }

    if (scanner == "cosmxscan") {
        GiottoUtils::vmsg("[cosmx_expression_disk] scanning .gz ->",
                          "parquetExprStore", .v = verbose)
        inp <- .cosmx_scan_input(path = path, slide = slide, fovs = .fovs)
        pe <- tryCatch(
            sourceWrite(gsource, inp, store_type = "parquetExpr",
                        verbose = verbose, ...),
            error = function(e) {
                if (!grepl("assumes integer counts", conditionMessage(e),
                           fixed = TRUE)) {
                    stop(e)
                }
                stop("[cosmx_expression_disk] ", basename(path), " has ",
                     "non-integer values after its first row; read it with ",
                     "scanner = \"csv\"", call. = FALSE)
            }
        )
        # fov and cell_ID arrive with the counts, so the global IDs are
        # composed during the stream and need no second pass
    }
    if (scanner == "csv") {
        GiottoUtils::vmsg("[cosmx_expression_disk] streaming CSV ->",
                          "parquetExprStore", .v = verbose)
        inp <- csvWideInput(
            csv_path        = path,
            cell_id_col     = "cell_ID",
            skip_cols       = c("fov", drop),
            row_filter_fun  = row_filter
        )
        pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                          verbose = verbose, ...)

        # Reconstruct globally-unique cell IDs `c_<slide>_<fov>_<cell_ID>`
        # by re-reading just the (fov, cell_ID) columns from the source CSV.
        cell_ID <- NULL  # NSE binding
        id_dt <- data.table::fread(path, select = c("fov", "cell_ID"))
        id_dt <- id_dt[cell_ID != 0L, ]
        if (!is.null(.fovs)) id_dt <- id_dt[fov %in% .fovs, ]
        ids <- sprintf("c_%d_%d_%d", slide, id_dt$fov, id_dt$cell_ID)
        ids <- ids[!ids %in% na_ids]
        if (length(ids) != length(pe@cell_ids)) {
            stop("[cosmx_expression_disk] cell-row mismatch: filtered CSV ",
                 "rows = ", length(ids), ", parquetExprStore cells = ",
                 length(pe@cell_ids), ". Filter logic disagrees with the ",
                 "streamed write.", call. = FALSE)
        }
        pe@cell_ids <- ids
        if (length(na_ids)) {
            warning(GiottoUtils::wrap_txt(sprintf(
                "[cosmx_expression_disk] %d %s in %s %s NA values and %s left
                out, as cells missing from it (e.g. %s)",
                length(na_ids), if (length(na_ids) == 1L) "cell" else "cells",
                basename(path), if (length(na_ids) == 1L) "has" else "have",
                if (length(na_ids) == 1L) "is" else "are",
                paste(utils::head(na_ids, 2L), collapse = ", ")
            )), call. = FALSE)
        }
    }
    if (!length(pe@cell_ids)) {
        stop("[cosmx_expression_disk] no cells in ", basename(path),
             if (!is.null(.fovs)) {
                 paste0(" for fovs ", paste(.fovs, collapse = ", "))
             }, call. = FALSE)
    }

    .cosmx_expr_split(pe, feat_type, split_keyword, output,
                      verbose = verbose)
}

.cosmx_csv_header <- function(path) {
    con <- if (.is_real_gz(path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    hdr <- strsplit(readLines(con, n = 1L, warn = FALSE), ",",
                    fixed = TRUE)[[1L]]
    gsub("^\"|\"$", "", hdr)
}

# cosmxscan parses integers only: it stops at `.`, `e`, a quote or `NA`
.cosmx_values_are_integer <- function(path) {
    con <- if (.is_real_gz(path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    l <- readLines(con, n = 2L, warn = FALSE)
    length(l) < 2L || !grepl("[^0-9,\r-]", l[[2L]])
}


.cosmx_expr_split <- function(pe, feat_type, split_keyword, output,
                              verbose = NULL) {
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
            if (any(bool)) {
                expr_list[[key_i + 1L]] <- pe[which(bool), , drop = FALSE]
                remaining <- remaining & !bool
            }
        }
        expr_list[[1L]] <- pe[which(remaining), , drop = FALSE]
        empty <- vapply(expr_list, is.null, logical(1L))
        if (any(empty)) {
            GiottoUtils::vmsg("  no features for feat_type ",
                              paste(names(expr_list)[empty], collapse = ", "),
                              .v = verbose)
        }
        expr_list <- expr_list[!empty]
    }

    GiottoUtils::vmsg(
        sprintf("  feature split: %s",
                paste(sprintf("%s %s", names(expr_list),
                              vapply(expr_list,
                                     function(x) length(x@feat_ids), 1L)),
                      collapse = " | ")),
        .v = verbose
    )

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


## cells and metadata ####

.cosmx_read_cellmeta <- function(path, cols, ids = NULL, fovs = NULL,
                                 cores = 1L, verbose = NULL) {
    if (!length(path) || !file.exists(path)) return(NULL)
    cell_ID <- fov <- NULL # NSE bindings
    hdr <- names(data.table::fread(path, nrows = 0L))
    if (!"cell" %in% hdr) {
        GiottoUtils::vmsg("[CosMxDiskReader] metadata has no",
                          "global `cell` column -- skipping",
                          .v = verbose)
        return(NULL)
    }
    cx <- data.table::fread(
        path,
        select = intersect(unique(c("cell", "fov", cols)), hdr),
        nThread = cores
    )
    if (length(fovs) && "fov" %in% names(cx)) {
        cx <- cx[fov %in% as.integer(fovs), ]
    }
    if (!"fov" %in% cols && "fov" %in% names(cx)) cx[, fov := NULL]
    data.table::setnames(cx, "cell", "cell_ID")
    if (!is.null(ids)) cx <- cx[cell_ID %in% ids, ]
    cx
}

.cosmx_add_fov <- function(g, feat_type = NULL) {
    pd <- GiottoClass::pDataDT(g, feat_type = feat_type)
    if ("fov" %in% names(pd) || !length(pd$cell_ID)) return(g)
    fv <- sub("^c_\\d+_(\\d+)_\\d+$", "\\1", pd$cell_ID)
    # sub() leaves an unparsed id unchanged, which reads as a constant column
    if (any(fv == pd$cell_ID)) {
        stop("[CosMxDiskReader] cell_IDs are not of the form ",
             "c_<slide>_<fov>_<cell>", call. = FALSE)
    }
    GiottoClass::addCellMetadata(
        g, feat_type = feat_type,
        new_metadata = data.table::data.table(
            cell_ID = pd$cell_ID, fov = as.integer(fv)
        ),
        by_column = TRUE, column_cell_ID = "cell_ID"
    )
}


## polygon ####

# Disk-backed CosMx polygon ingestion, following `.tenx_poly_disk()`:
# a `parquetStore` intermediate first (polygons need `row_index` for stable
# vertex ordering), then a `parquetGeomTileStore`.
#
# CosMx specifics:
#   - polygons are keyed on `cell` (`c_<slide>_<fov>_<cellID>`); `cellID`
#     restarts per FOV.
#   - coordinates are read as double, set in the CSV reader: the tile writer
#     finds the x/y columns by type, and AtoMx writes round values as `1e+05`.
#   - `fovs` and `keep_ids` are subset() ops after `row_index` is fixed, so
#     vertex order is kept.
.cosmx_poly_disk <- function(
    path,
    gsource,
    name = "cell",
    part_col = NULL,
    flip_vertical = FALSE,
    fovs = NULL,
    keep_ids = NULL,
    expr_ids = keep_ids,
    output = c("giottoPolygon", "store"),
    read_fun = NULL,
    verbose = NULL,
    ...
) {
    checkmate::assert_file_exists(path)
    checkmate::assert_class(gsource, "gsource")
    checkmate::assert_character(name, len = 1L)
    output <- match.arg(output, choices = c("giottoPolygon", "store"))
    GiottoUtils::package_check("arrow")
    GiottoUtils::package_check("dplyr")

    GiottoUtils::vmsg(
        sprintf("[cosmx_poly_disk] loading boundary '%s'", name),
        .v = verbose
    )

    fmt <- if (any(grepl("[.]parquet$", path, ignore.case = TRUE))) {
        "parquet"
    } else {
        "csv"
    }
    fov_done <- FALSE
    if (length(fovs) && fmt == "csv") {
        sub_path <- .cosmx_csv_fov_subset(path, as.integer(fovs))
        if (!is.null(sub_path)) {
            on.exit(unlink(sub_path), add = TRUE)
            path <- sub_path
            fov_done <- TRUE
        }
    }

    # No dplyr filters or row-wise mutations beyond the type casts: they push
    # arrow into a threaded scan that delivers batches out of source order,
    # which scatters a cell's vertices and breaks ring construction.
    if (is.null(read_fun)) {
        # the two global columns are the ones that matter; fov / cellID /
        # x_local / y_local are typed only to keep parquetStore from warning
        # that int64 is poorly supported
        types <- list(
            x_global_px = arrow::float64(),
            y_global_px = arrow::float64(),
            x_local_px  = arrow::float64(),
            y_local_px  = arrow::float64(),
            fov         = arrow::int32(),
            cellID      = arrow::int32()
        )
        if (fmt == "csv") {
            hdr <- names(data.table::fread(path, nrows = 0L))
            col_types <- do.call(arrow::schema, types[names(types) %in% hdr])
            read_fun <- function(x, ...) {
                arrow::open_dataset(sources = x, format = "csv",
                                    col_types = col_types)
            }
        } else {
            # NSE bindings
            x_global_px <- y_global_px <- x_local_px <- y_local_px <- NULL
            fov <- cellID <- NULL
            read_fun <- function(x, ...) {
                a <- arrow::open_dataset(sources = x, format = fmt)
                dplyr::mutate(a,
                    x_global_px = as.double(x_global_px),
                    y_global_px = as.double(y_global_px),
                    x_local_px  = as.double(x_local_px),
                    y_local_px  = as.double(y_local_px),
                    fov         = as.integer(fov),
                    cellID      = as.integer(cellID)
                )
            }
        }
    }

    fs <- fileStore(path = path, read_fun = read_fun)
    qs <- methods::as(fs, "queryableStore")

    poly_intermediate <- sourceWrite(gsource, qs, store_type = "parquet",
                                     verbose = verbose)

    if (length(fovs) && !fov_done) {
        if ("fov" %in% colnames(poly_intermediate)) {
            keep_fovs <- as.integer(fovs)
            poly_intermediate <- subset(poly_intermediate, fov %in% keep_fovs)
        } else if (is.null(keep_ids)) {
            warning("[cosmx_poly_disk] the polygon file has no `fov` column; ",
                    "`fovs` is not applied to the polygons", call. = FALSE)
        }
    }
    if (!is.null(keep_ids)) {
        poly_intermediate <- .cosmx_poly_keep(poly_intermediate,
                                              keep_ids = keep_ids,
                                              expr_ids = expr_ids)
    }

    # part_col = NULL: CosMx emits one ring per cell, so there is nothing to
    # detect -- and a detector's regex could latch onto `cellID`.
    poly_store <- sourceWrite(
        gsource, poly_intermediate,
        store_type = "parquetGeomTile",
        type = "polygons",
        id_col = "cell",
        sdimx = "x_global_px",
        sdimy = "y_global_px",
        group_col = "cell",
        part_col = part_col,
        flip_vertical = isTRUE(flip_vertical),
        verbose = verbose,
        ...
    )

    if (output == "store") return(poly_store)

    # Centroids are computed lazily from the parquetGeomBase by
    # setGiotto(centroids_to_spatlocs = TRUE) at attach time.
    GiottoClass::createGiottoPolygon(x = poly_store, name = name)
}

# lines, not an arrow filter: arrow scans in threads and loses the row order
.cosmx_csv_fov_subset <- function(path, fovs, chunk_lines = 1e6) {
    hdr <- .cosmx_csv_header(path)
    j <- match("fov", hdr)
    if (is.na(j)) return(NULL)
    con <- if (.is_real_gz(path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    out <- tempfile("cosmx_fovs_", fileext = ".csv")
    con_out <- file(out, "wt")
    on.exit(close(con_out), add = TRUE)
    writeLines(readLines(con, n = 1L, warn = FALSE), con_out)
    repeat {
        lines <- readLines(con, n = chunk_lines, warn = FALSE)
        if (!length(lines)) break
        lines <- lines[nzchar(lines)]
        if (!length(lines)) next
        fov <- data.table::fread(text = lines, header = FALSE, select = j,
                                 sep = ",")[[1L]]
        writeLines(lines[fov %in% fovs], con_out)
    }
    out
}

# The filter op is safe for ring construction: `row_index` is already baked
# in and every tile is re-sorted on it.
.cosmx_poly_keep <- function(x, keep_ids, expr_ids = keep_ids) {
    cell <- NULL # NSE binding
    keep_ids <- as.character(keep_ids)
    all_ids <- dplyr::collect(dplyr::distinct(storeRead(x), cell))$cell
    drop_ids <- setdiff(all_ids, keep_ids)
    if (!length(drop_ids)) return(x)

    orphan_ids <- setdiff(drop_ids, as.character(expr_ids))
    if (length(orphan_ids)) {
        warning(GiottoUtils::wrap_txt(sprintf(
            "[cosmx_poly_disk] %d of %d polygons belong to cells with no
            expression (for example %s) and are not loaded.",
            length(orphan_ids), length(all_ids),
            paste(utils::head(sort(orphan_ids), 3L), collapse = ", ")
        )), call. = FALSE)
    }

    # a filter op through subset(): on the dropped ids, usually far fewer
    # than the kept ones
    subset(x, cell %in% drop_ids, negate = TRUE)
}


## transcript ####

# Disk-backed CosMx transcript ingestion, following `.atera_transcript_disk()`.
# No columns are dropped by default: `cell_ID` / `cell`, the vendor's cell of
# each detection, allow checking the overlaps.
.cosmx_transcript_disk <- function(
    path,
    gsource,
    version = "v6",
    fovs = NULL,
    feat_type = c("rna", "negprobes", "falsecode"),
    split_keyword = list("^Negative", "^SystemControl"),
    dropcols = c(),
    cores = GiottoUtils::determine_cores(),
    output = c("giottoPoints", "store"),
    verbose = NULL,
    ...
) {
    if (missing(path) || length(path) == 0L) {
        stop("[cosmx_transcript_disk] no tx_file path provided", call. = FALSE)
    }
    checkmate::assert(
        checkmate::check_file_exists(path),
        checkmate::check_directory_exists(path)
    )
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("giottoPoints", "store"))
    GiottoUtils::package_check("arrow")

    path <- .cosmx_ensure_parquet(path, version = version, verbose = verbose)

    GiottoUtils::vmsg("[cosmx_transcript_disk] loading transcripts",
        .v = verbose)
    tx_store <- sourceWrite(
        gsource, cosmxTxInput(path),
        store_type = "parquetGeomTile",
        fovs = if (!is.null(fovs)) as.integer(fovs),
        dropcols = as.character(dropcols),
        workers = cores,
        verbose = verbose,
        ...
    )

    .tx_points_out(tx_store, output, feat_type, split_keyword)
}

# Pinned schema: the leading rows (`cell_ID == 0`) have a blank `CellComp`,
# which arrow infers as `null`; AtoMx writes round coordinates as `1e+05`.
.cosmx_tx_schema <- function(version = c("v6", "legacy")) {
    version <- match.arg(tolower(version), choices = c("v6", "legacy"))
    i <- arrow::int32()
    d <- arrow::float64()
    s <- arrow::utf8()
    switch(version,
        "v6" = arrow::schema(
            fov = i, cell_ID = i, cell = s,
            x_local_px = i, y_local_px = i,
            x_global_px = d, y_global_px = d,
            z = i, target = s, CellComp = s
        ),
        # legacy has no `cell` column and orders global before local
        "legacy" = arrow::schema(
            fov = i, cell_ID = i,
            x_global_px = d, y_global_px = d,
            x_local_px = d, y_local_px = d,
            z = i, target = s, CellComp = s
        )
    )
}

# bump when converter output changes; part of the cache fingerprint so
# stale conversions are redone rather than reused
.cosmx_converter_version <- "0.2.0"

# Convert-once wrapper for the tx_file, mirroring `.zarr_ensure_parquet()`:
# streams it to parquet partitioned by FOV, cached under
# <getArtifactDumpDir()>/cosmx_cache/<fingerprint>/ (path, size, mtime and
# versions). manifest.json is written last, so a half-finished conversion is
# redone.
.cosmx_ensure_parquet <- function(path,
    version = "v6",
    compression = "zstd",
    verbose = NULL) {
    if (.cosmx_fov_partitioned(path)) return(path)
    is_parquet <- dir.exists(path) ||
        grepl("[.]parquet$", path, ignore.case = TRUE)
    info <- file.info(path)
    fp <- .hash(list(
        path = normalizePath(path, winslash = "/", mustWork = FALSE),
        size = as.numeric(info$size),
        mtime = as.numeric(info$mtime),
        version = tolower(version),
        converter = .cosmx_converter_version
    ))
    cache_dir <- file.path(getArtifactDumpDir(), "cosmx_cache", fp)
    manifest_path <- file.path(cache_dir, "manifest.json")
    target <- file.path(cache_dir, "transcripts")

    if (file.exists(manifest_path) && dir.exists(target)) {
        GiottoUtils::vmsg("[cosmx_transcript_disk] using cached conversion:",
            target, .v = verbose)
        return(target)
    }

    GiottoUtils::vmsg("[cosmx_transcript_disk] converting tx_file to parquet:",
        path, .v = verbose)
    unlink(target, recursive = TRUE)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    ds <- if (is_parquet) {
        arrow::open_dataset(path)
    } else {
        arrow::open_dataset(path,
            format = "csv",
            schema = .cosmx_tx_schema(version),
            skip = 1L # the pinned schema replaces the header
        )
    }
    arrow::write_dataset(ds, target,
        format = "parquet", compression = compression, partitioning = "fov")

    jsonlite::write_json(list(
        source = normalizePath(path, winslash = "/", mustWork = FALSE),
        version = tolower(version),
        converter = .cosmx_converter_version,
        products = list(transcripts = "transcripts")
    ), manifest_path, auto_unbox = TRUE, pretty = TRUE)
    target
}
