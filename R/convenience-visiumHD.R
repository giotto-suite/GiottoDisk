# VisiumHD ingest pipeline for `gDirSource`-managed projects.
#
# Disk-backed counterpart to Giotto's `VisiumHDReader`. Currently routes
# only the expression matrix through GiottoDisk (`parquetExprStore`
# written into the project vault via `sourceWrite`). Other modalities
# (tissue positions, scalefactors, images, polygons, transcripts,
# cellmeta) remain on the inherited in-mem closures and can be ported
# following the same pattern when needed.



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

        gsrc <- obj@backend
        bin_ <- obj@bin
        barcodes_ <- obj@barcodes
        feature_id_type_ <- obj@feature_id_type

        # expression (disk override). Mirrors the parent closure's signature
        # so the inherited gobject_fun can plumb identical args through.
        ex_fun <- function(
            path,
            outdir = obj@outdir,
            bin = bin_,
            expression_source = obj@expression_source,
            barcodes = barcodes_,
            feature_id_type = feature_id_type_,
            remove_zero_rows = TRUE,
            split_by_type = TRUE,
            output = c("exprObj", "store"),
            agg_type = c("bin", "cell"),
            verbose = NULL,
            ...
        ) {
            .visiumhd_expression_disk(
                path = path,
                gsource = gsrc,
                bin = bin,
                barcodes = barcodes,
                feature_id_type = feature_id_type,
                remove_zero_rows = remove_zero_rows,
                split_by_type = split_by_type,
                output = output,
                agg_type = match.arg(agg_type),
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_expression <- ex_fun

        obj
    }
)



# CREATE READER ####

#' @title Import a 10x VisiumHD assay (disk-backed)
#' @name importVisiumHDDisk
#' @description
#' Disk-backed counterpart to [Giotto::importVisiumHD()]. Produces a
#' `VisiumHDDiskReader` whose `load_expression()` call writes a
#' `parquetExprStore` into a `gDirSource`-managed project vault.
#' Other modalities (tissue positions, scalefactors, images, polygons,
#' transcripts, cellmeta) remain in-memory via the inherited
#' `VisiumHDReader` closures.
#' @param visiumhd_dir VisiumHD output directory
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#'   Naming matches [GiottoClass::createGiottoObject()]'s `backend` param.
#' @param bin,micron,outdir,expression_source,feature_id_type,tissue_only,barcodes,array_subset_row,array_subset_col,pxl_subset_row,pxl_subset_col,filter,filter_coverage_cutoff
#'   passed through to the parent `VisiumHDReader` initializer.
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
    a <- list(
        Class = "VisiumHDDiskReader",
        backend = backend,
        bin = as.integer(bin),
        micron = as.logical(micron),
        expression_source = match.arg(expression_source, c("raw", "filtered")),
        feature_id_type = match.arg(feature_id_type, c("symbol", "ensembl")),
        tissue_only = as.logical(tissue_only),
        filter_coverage = filter_coverage_cutoff
    )
    if (!is.null(visiumhd_dir))     a$visiumhd_dir     <- visiumhd_dir
    if (!is.null(outdir))           a$outdir           <- outdir
    if (!is.null(barcodes))         a$barcodes         <- barcodes
    if (!is.null(array_subset_row)) a$array_subset_row <- array_subset_row
    if (!is.null(array_subset_col)) a$array_subset_col <- array_subset_col
    if (!is.null(pxl_subset_row))   a$pxl_subset_row   <- pxl_subset_row
    if (!is.null(pxl_subset_col))   a$pxl_subset_col   <- pxl_subset_col
    if (!is.null(filter))           a$filter           <- filter
    do.call(new, args = a)
}



# MODULAR ####


## expression ####

# Disk-backed VisiumHD expression ingestion.
#
# Same 10x mtx-triple / .h5 format detection as the Xenium disk helper,
# plus VisiumHD-specific spat_unit naming ("bin016" / "bin008" / "cell")
# and an optional `barcodes` filter (tissue-only subset). Routes through
# sourceWrite(gsource, inp, store_type = "parquetExpr") with mtxInput /
# tenxH5Input. Slicing by feat_type and barcode is applied lazily via
# parquetExprStore `[i, j]` -- no parquet rewrite.
.visiumhd_expression_disk <- function(
    path,
    gsource,
    bin = 8L,
    barcodes = NULL,
    feature_id_type = c("symbol", "ensembl"),
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    output = c("exprObj", "store"),
    agg_type = c("bin", "cell"),
    verbose = NULL,
    ...
) {
    if (missing(path) || length(path) == 0L || !file.exists(path)) {
        stop("[visiumhd_expression_disk] no expression path provided",
             call. = FALSE)
    }
    checkmate::assert_class(gsource, "gsource")
    output <- match.arg(output, choices = c("exprObj", "store"))
    agg_type <- match.arg(agg_type)
    feature_id_type <- match.arg(feature_id_type, c("symbol", "ensembl"))

    feature_id_col <- switch(feature_id_type,
        "ensembl" = 1L,
        "symbol"  = 2L
    )

    # Detect format
    if (dir.exists(path)) {
        fmt <- "mtx"
    } else {
        ext <- tolower(tools::file_ext(path))
        fmt <- switch(ext,
            "h5"  = "h5",
            "mtx" = "mtx",
            "gz"  = "mtx",
            stop("[visiumhd_expression_disk] unsupported format: ", path,
                 call. = FALSE)
        )
    }
    GiottoUtils::vmsg("[visiumhd_expression_disk] format:", fmt,
                       " agg:", agg_type, .v = verbose)

    if (fmt == "mtx") {
        inp <- mtxInput(path, feature_id_col = feature_id_col)
        feat_classes_vec <- .tenx_feat_classes_mtx(path)
    } else {
        inp <- tenxH5Input(path, feature_id_col = feature_id_col)
        feat_classes_vec <- .tenx_feat_classes_h5(path)
    }

    pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                       verbose = verbose, ...)

    # remove_zero_rows via Arrow distinct(col_id) (single streaming pass)
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

    # barcodes filter: lazy column subset on pe
    if (!is.null(barcodes)) {
        bool <- pe@cell_ids %in% barcodes
        pe <- pe[, bool, drop = FALSE]
    }

    # spat_unit / provenance naming: VisiumHD convention is "bin<NNN>"
    # for binned reads, "cell" for segmented reads.
    su <- if (agg_type == "bin") sprintf("bin%03d", as.integer(bin)) else "cell"

    # 10x feature_type -> Giotto feat_type rename (mirrors the in-mem path)
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
        nm <- if (length(uniq_classes) == 1L) feat_class_to_name(uniq_classes)
              else "rna"
        store_list <- stats::setNames(list(pe), nm)
    }

    if (output == "store") return(store_list)

    lapply(seq_along(store_list), function(i) {
        methods::new("exprObj",
            name       = "raw",
            exprMat    = store_list[[i]],
            spat_unit  = su,
            feat_type  = names(store_list)[[i]],
            provenance = su
        )
    })
}
