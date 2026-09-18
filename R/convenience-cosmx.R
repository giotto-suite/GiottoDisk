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
        # Dotted copies: a closure default written as `meta_path = meta_path`
        # resolves against the closure's OWN formals first and self-references,
        # so the closures below bind these instead.
        .poly_path <- obj@paths$poly_path
        .meta_path <- obj@paths$meta_path
        .expr_path <- obj@paths$expr_path
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

        # polygons (disk override)
        poly_fun <- function(
            path = .poly_path,
            name = "cell",
            part_col = NULL,
            flip_vertical = FALSE,
            calc_centroids = TRUE,
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
                output = output,
                verbose = verbose,
                ...
            )
        }
        obj@calls$load_polys <- poly_fun

        # create_gobject (disk variant). Mirrors parent's gobject_fun but
        # initializes the giotto object with backend = gsrc so any further
        # artifacts are vault-resident.
        gobject_fun <- function(
            polygon_path = .poly_path,
            expression_path = .expr_path,
            metadata_path = .meta_path,
            expr_store = NULL,
            poly_read_fun = NULL,
            feat_type = c("rna", "negprobes", "falsecode"),
            split_keyword = list("^Negative", "^SystemControl"),
            meta_cols = c("fov", "Area.um2", "nCount_RNA", "nFeature_RNA"),
            load_polygons = TRUE,
            load_expression = TRUE,
            load_cellmeta = TRUE,
            instructions = NULL,
            cores = GiottoUtils::determine_cores(),
            verbose = NULL
        ) {
            funs <- obj@calls

            # init gobject with disk backend
            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            # polygons (disk; overridden closure). Also the source of the cell
            # IDs and of the centroids that become spatlocs.
            allowed_ids <- NULL
            if (isTRUE(load_polygons)) {
                polys <- funs$load_polys(
                    path = polygon_path, name = "cell",
                    read_fun = poly_read_fun, verbose = verbose
                )
                allowed_ids <- GiottoClass::spatIDs(polys)
                g <- GiottoClass::setGiotto(
                    g, polys, centroids_to_spatlocs = TRUE, verbose = FALSE
                )
            }

            # expression, sliced to the cells that actually carry a boundary.
            # `expr_store` is a parquetExprStore already written by
            # sourceWrite(); without one the inherited closure streams the
            # wide CSV instead.
            if (isTRUE(load_expression)) {
                exlist <- if (!is.null(expr_store)) {
                    .cosmx_expr_split(expr_store, feat_type, split_keyword,
                                      "exprObj", verbose = verbose)
                } else {
                    funs$load_expression(
                        path = expression_path,
                        feat_type = feat_type,
                        split_keyword = split_keyword,
                        verbose = verbose
                    )
                }
                for (ex in exlist) {
                    if (!is.null(allowed_ids)) {
                        bool <- colnames(ex[]) %in% allowed_ids
                        if (!any(bool)) {
                            stop("[CosMxDiskReader] expression and polygon ",
                                 "cell IDs do not intersect; both should be ",
                                 "'c_<slide>_<fov>_<cell>'.", call. = FALSE)
                        }
                        ex[] <- ex[][, bool]
                    }
                    g <- GiottoClass::setGiotto(g, ex, verbose = FALSE)
                }
            }

            # cellmeta, read from the vendor CSV rather than through the
            # parent closure: `cell` is the vendor's GLOBAL id and the only
            # column that joins to the polygons, while `cell_ID` in that same
            # file is the FOV-local integer.
            if (isTRUE(load_cellmeta) && !is.null(metadata_path) &&
                file.exists(metadata_path)) {
                cell_ID <- NULL # NSE binding
                hdr <- names(data.table::fread(metadata_path, nrows = 0L))
                if (!"cell" %in% hdr) {
                    GiottoUtils::vmsg("[CosMxDiskReader] metadata has no",
                                      "global `cell` column -- skipping",
                                      .v = verbose)
                } else {
                    cx <- data.table::fread(
                        metadata_path,
                        select = intersect(c("cell", meta_cols), hdr),
                        nThread = cores
                    )
                    data.table::setnames(cx, "cell", "cell_ID")
                    if (!is.null(allowed_ids)) {
                        cx <- cx[cell_ID %in% allowed_ids, ]
                    }
                    g <- GiottoClass::addCellMetadata(
                        g, new_metadata = cx,
                        by_column = TRUE, column_cell_ID = "cell_ID"
                    )
                }
            }

            # add fovs metadata column. Never `g$fov <- ...` here -- on this
            # object that is a silent no-op, no error and no column.
            pd <- GiottoClass::pDataDT(g)
            if (!"fov" %in% names(pd) && length(pd$cell_ID)) {
                fv <- sub("^c_\\d+_(\\d+)_\\d+$", "\\1", pd$cell_ID)
                # sub() returns its input unchanged when it does not match, so
                # an unparsed id would read as a constant column. Fail loudly.
                if (any(fv == pd$cell_ID)) {
                    stop("[CosMxDiskReader] cell_IDs are not of the form ",
                         "c_<slide>_<fov>_<cell>", call. = FALSE)
                }
                g <- GiottoClass::addCellMetadata(
                    g,
                    new_metadata = data.table::data.table(
                        cell_ID = pd$cell_ID, fov = as.integer(fv)
                    ),
                    by_column = TRUE, column_cell_ID = "cell_ID"
                )
            }

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
#'   parent `CosmxReader` initializer. `poly_pref` defaults to `"csv"`
#'   because the disk polygon loader reads the vertex CSV, not masks.
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
                              poly_pref = c("csv", "mask")) {
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

    .cosmx_expr_split(pe, feat_type, split_keyword, output,
                      verbose = verbose)
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


## polygon ####

# Disk-backed CosMx polygon ingestion, following `.xenium_poly_disk()`:
# a `parquetStore` intermediate first (polygons need `row_index` for stable
# vertex ordering), then a `parquetGeomTileStore`.
#
# Two CosMx-specific points:
#   - `id_col` / `group_col` are `cell` (the global `c_<slide>_<fov>_<cellID>`
#     form). The file also carries `cellID`, which restarts per FOV; keying on
#     it silently mis-groups vertices.
#   - the double cast in read_fun IS required, and is the one mutation this
#     loader performs. CosMx ships pixel coordinates as int64 and the geometry
#     tile writer identifies the x/y columns BY TYPE: with integer coords it
#     picks the wrong pair, coerces the string ids ("NAs introduced by
#     coercion"), never creates `poly_ID`, and dies inside tileApply with
#     "setcolorder: non-existing column(s): cols[6]='poly_ID'".
.cosmx_poly_disk <- function(
    path,
    gsource,
    name = "cell",
    part_col = NULL,
    flip_vertical = FALSE,
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

    # No dplyr filters or row-wise mutations beyond the type casts: they push
    # arrow into a threaded scan that delivers batches out of source order,
    # which scatters a cell's vertices and breaks ring construction.
    if (is.null(read_fun)) {
        read_fun <- function(x, ...) {
            a <- arrow::open_dataset(sources = x, format = fmt)
            # the two global columns are the ones that matter; fov / cellID /
            # x_local / y_local are cast only to keep parquetStore from
            # warning that int64 is poorly supported
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

    fs <- fileStore(path = path, read_fun = read_fun)
    qs <- methods::as(fs, "queryableStore")

    poly_intermediate <- sourceWrite(gsource, qs, store_type = "parquet",
                                     verbose = verbose)

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
