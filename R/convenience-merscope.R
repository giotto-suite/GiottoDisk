# MERSCOPE (Vizgen MERFISH) ingest pipeline for `gDirSource`-managed projects.
#
# Parallels `convenience-xenium.R`, but standalone: the other disk readers
# subclass an in-memory `*Reader` from {Giotto} and inherit its path detection,
# while {Giotto}'s MERSCOPE support is the older function-based
# `createGiottoMerscopeObject()`. This reader owns its `@paths` detection and
# `@calls` closures.
#
# Modality routing (matching the Xenium disk reader):
#   transcripts -> parquetGeomTile store (points)
#   polygons    -> parquetGeom -> parquetGeomTile store
#   expression  -> parquetExpr store
#   cell meta   -> in memory
#
# Z-planes: MERSCOPE writes one boundary row per (cell, z-plane) and exports are
# usually replicated 2D. `.merscope_zplane_architecture()` detects that and
# drops the redundant planes; genuine 3D is refused. Rationale: adr/0010.



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

# * accessors ####
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
            # mirror createGiottoObject(backend = ...) handling: a character
            # input is treated as a directory path and coerced to a gDirSource.
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

        # Mirror @paths into the init frame so the closures below can reference
        # path names directly via default-arg expressions. Each default
        # evaluates to a bare character string at call time, so nothing
        # scope-laddering back to this frame travels into parallel workers.
        # (Same convention as XeniumDiskReader.)
        list2env(.Object@paths, envir = environment())
        gsrc     <- .Object@backend
        score_th <- .Object@score_threshold
        poly_fmt <- .Object@polygon_format

        # transcripts
        tx_fun <- function(
            path = tx_path,
            feat_type = c("rna", "Blank"),
            split_keyword = list("Blank"),
            FOVs = NULL,
            dropcols = .MERSCOPE_TX_DROPCOLS,
            score_threshold = score_th,
            flip_axis = "none",
            flip_anchor = NULL,
            output = c("giottoPoints", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_transcript_disk(
                path = path, gsource = gsrc, feat_type = feat_type,
                split_keyword = split_keyword, FOVs = FOVs,
                dropcols = dropcols, score_threshold = score_threshold,
                flip_axis = flip_axis, flip_anchor = flip_anchor,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_transcripts <- tx_fun

        # polygons
        poly_fun <- function(
            path = cell_bound_path,
            name = "cell",
            poly_z_indices = NULL,
            polygon_format = poly_fmt,
            FOVs = NULL,
            metadata_path = cell_meta_path,
            flip_axis = "none",
            flip_anchor = NULL,
            output = c("giottoPolygon", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_poly_disk(
                path = path, gsource = gsrc, name = name,
                poly_z_indices = poly_z_indices,
                polygon_format = polygon_format,
                FOVs = FOVs, metadata_path = metadata_path,
                flip_axis = flip_axis, flip_anchor = flip_anchor,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_polys <- poly_fun

        # expression
        ex_fun <- function(
            path = expr_path,
            split_by_type = TRUE,
            split_keyword = "Blank",
            remove_zero_rows = TRUE,
            output = c("exprObj", "store"),
            verbose = NULL,
            ...
        ) {
            .merscope_expression_disk(
                path = path, gsource = gsrc, split_by_type = split_by_type,
                split_keyword = split_keyword,
                remove_zero_rows = remove_zero_rows,
                output = output, verbose = verbose, ...
            )
        }
        .Object@calls$load_expression <- ex_fun

        # cell metadata (in memory)
        cm_fun <- function(path = cell_meta_path, FOVs = NULL, verbose = NULL, ...) {
            .merscope_cellmeta(path = path, FOVs = FOVs, verbose = verbose, ...)
        }
        .Object@calls$load_cellmeta <- cm_fun

        # orchestrator
        gobject_fun <- function(
            transcript_path = tx_path,
            boundary_path = cell_bound_path,
            expression_path = expr_path,
            metadata_path = cell_meta_path,
            data_to_use = c("subcellular", "aggregate"),
            feat_type = c("rna", "Blank"),
            split_keyword = list("Blank"),
            poly_z_indices = NULL,
            FOVs = NULL,
            flip_axis = c("none", "y", "x", "both"),
            spat_unit = "cell",
            load_transcripts = TRUE,
            load_expression = TRUE,
            load_cellmeta = FALSE,
            instructions = NULL,
            verbose = NULL) {

            data_to_use <- match.arg(data_to_use, c("subcellular", "aggregate"))
            flip_axis <- match.arg(flip_axis, c("none", "y", "x", "both"))
            load_transcripts <- as.logical(load_transcripts)
            load_expression  <- as.logical(load_expression)
            load_cellmeta    <- as.logical(load_cellmeta)

            # `aggregate` uses the vendor cell x gene matrix only; there is no
            # subcellular geometry to place, so transcripts/polygons are skipped
            # entirely. Mirrors `data_to_use` in createGiottoMerscopeObject().
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

            # Reflection anchor. Giotto's flip_axis reflects about the midpoint
            # of the *transcript* coordinate range so that the polygon and
            # transcript layers stay registered with one another. That anchor
            # therefore has to be computed once, up front, from the transcript
            # extent and handed to both loaders -- neither layer may pick its
            # own.
            flip_anchor <- NULL
            if (flip_axis != "none") {
                flip_anchor <- .merscope_flip_anchor(
                    transcript_path, FOVs = FOVs, verbose = verbose
                )
                GiottoUtils::vmsg(.v = verbose, sprintf(
                    "[merscope] flip_axis '%s' about (x0 = %.2f, y0 = %.2f)",
                    flip_axis, flip_anchor[["x0"]], flip_anchor[["y0"]]
                ))
            }

            # ---- load first, construct the giotto object last ----------------
            #
            # ORDERING IS LOAD-BEARING. Xenium's order -- createGiottoObject(
            # backend =) first, then load modalities -- segfaults here (Windows
            # 0xC0000005) once a large Arrow scan runs afterwards. Reproducible
            # on the 151M-row detected_transcripts.parquet at any output size,
            # so not OOM; unaffected by giottodisk.duckdb_memory_limit,
            # arrow::set_cpu_count(1) or set_io_thread_count(1). Only
            # constructing the gobject AFTER the scans avoids it. Small
            # synthetic sources do not trigger it, which is presumably why
            # Xenium has not hit it.
            #
            # Likely an upstream GiottoDisk issue, not a MERSCOPE one; revert to
            # the Xenium ordering once fixed.

            tx_list <- NULL
            if (load_transcripts) {
                tx_list <- funs$load_transcripts(
                    path = transcript_path,
                    feat_type = feat_type,
                    split_keyword = split_keyword,
                    FOVs = FOVs,
                    flip_axis = flip_axis,
                    flip_anchor = flip_anchor,
                    verbose = verbose
                )
            }

            b <- NULL
            if (data_to_use == "subcellular" && !is.na(boundary_path)) {
                b <- funs$load_polys(
                    path = boundary_path,
                    name = spat_unit,
                    poly_z_indices = poly_z_indices,
                    FOVs = FOVs,
                    metadata_path = metadata_path,
                    flip_axis = flip_axis,
                    flip_anchor = flip_anchor,
                    verbose = verbose
                )
                if (!is.null(b) && !inherits(b, "list")) b <- list(b)
            }

            ex <- NULL
            if (load_expression) {
                ex <- funs$load_expression(
                    path = expression_path,
                    split_by_type = TRUE,
                    split_keyword = if (is.list(split_keyword)) {
                        unlist(split_keyword)
                    } else split_keyword,
                    verbose = verbose
                )
            }

            cx <- NULL
            if (load_cellmeta) {
                cx <- funs$load_cellmeta(path = metadata_path, FOVs = FOVs,
                                         verbose = verbose)
            }

            # init gobject with disk backend, then attach
            g <- GiottoClass::createGiottoObject(
                backend = gsrc,
                instructions = instructions
            )

            if (!is.null(tx_list)) {
                g <- GiottoClass::setGiotto(g, tx_list, verbose = FALSE)
            }
            if (!is.null(b)) {
                # centroids_to_spatlocs = TRUE attaches centroids as spatlocs at
                # set time, computed lazily from the parquetGeomBase.
                g <- GiottoClass::setGiotto(g, b,
                    centroids_to_spatlocs = TRUE, verbose = FALSE)
            }
            if (!is.null(ex)) g <- GiottoClass::setGiotto(g, ex)
            if (!is.null(cx)) g <- GiottoClass::setGiotto(g, cx, verbose = FALSE)

            GiottoUtils::vmsg(.v = verbose, "done")
            g
        }
        .Object@calls$create_gobject <- gobject_fun

        .Object
    }
)



# CREATE READER ####

#' @title Import a Vizgen MERSCOPE assay (disk-backed)
#' @name importMerscopeDisk
#' @description
#' Disk-backed MERSCOPE reader. Produces a `MerscopeDiskReader` whose
#' `load_transcripts()`, `load_polys()` and `load_expression()` calls write to a
#' `gDirSource`-managed project vault. Cell metadata remains in memory.
#'
#' Unlike [importXeniumDisk()], this reader does not subclass an in-memory
#' {Giotto} reader -- no `MerscopeReader` S4 class exists -- so it performs its
#' own path detection.
#' @param merscope_dir MERSCOPE region directory (the one containing
#'   `detected_transcripts.parquet`). `experiment.json` is looked up one level
#'   above, where the vendor writes it.
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#' @param score_threshold `numeric` or `NULL` (default). Minimum
#'   `transcript_score` retained. Vendor exports are already score-filtered, so
#'   unlike Xenium's `qv_threshold` this defaults to no additional filtering.
#' @param polygon_format `character`. `"parquet"` (default) or `"hdf5"`.
#' @returns `MerscopeDiskReader` object
#' @seealso [importXeniumDisk()]
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


#' @title Create a disk-backed MERSCOPE Giotto object
#' @name createGiottoMerscopeObjectDisk
#' @description
#' One-call entry point mirroring the shape of
#' `Giotto::createGiottoXeniumObject(backend = )`, for Vizgen MERSCOPE output.
#' @param merscope_dir MERSCOPE region directory
#' @param backend a `gsource` project backend
#' @param score_threshold minimum `transcript_score`, or `NULL` (default)
#' @param polygon_format `"parquet"` (default) or `"hdf5"`
#' @param data_to_use `"subcellular"` (default) builds from transcripts and
#'   polygons; `"aggregate"` uses only the vendor cell x gene matrix.
#' @param poly_z_indices which z-indices to use for polygons. `NULL` (default)
#'   auto-detects the architecture -- see the z-plane note in this file.
#' @param FOVs which FOVs to load. `NULL` (default) loads all.
#' @param flip_axis axis along which to reflect polygons and transcripts after
#'   loading: `"none"` (default), `"y"`, `"x"`, `"both"`. The reflection origin
#'   is the midpoint of the transcript coordinate range so both layers stay
#'   registered.
#' @param split_keyword keywords separating feature detections. Default
#'   `list("Blank")`, the MERSCOPE negative-control prefix.
#' @param spat_unit name for the cell spatial unit (default `"cell"`)
#' @param load_transcripts,load_expression,load_cellmeta logical
#' @param instructions Giotto instructions
#' @param verbose verbosity
#' @returns a `giotto` object backed by `backend`
#' @export
createGiottoMerscopeObjectDisk <- function(
    merscope_dir,
    backend,
    score_threshold = NULL,
    polygon_format = c("parquet", "hdf5"),
    data_to_use = c("subcellular", "aggregate"),
    poly_z_indices = NULL,
    FOVs = NULL,
    flip_axis = c("none", "y", "x", "both"),
    split_keyword = list("Blank"),
    spat_unit = "cell",
    load_transcripts = TRUE,
    load_expression = TRUE,
    load_cellmeta = FALSE,
    instructions = NULL,
    verbose = NULL
) {
    x <- importMerscopeDisk(
        merscope_dir = merscope_dir,
        backend = backend,
        score_threshold = score_threshold,
        polygon_format = polygon_format
    )
    x$create_gobject(
        data_to_use = data_to_use,
        split_keyword = split_keyword,
        poly_z_indices = poly_z_indices,
        FOVs = FOVs,
        flip_axis = flip_axis,
        spat_unit = spat_unit,
        load_transcripts = load_transcripts,
        load_expression = load_expression,
        load_cellmeta = load_cellmeta,
        instructions = instructions,
        verbose = verbose
    )
}



# PATHS ####

# MERSCOPE lays a region out as:
#   <experiment>/
#     experiment.json          <- one level ABOVE the region dir
#     region_R1/
#       detected_transcripts.parquet
#       cell_boundaries.parquet   (or cell_boundaries/ for hdf5 exports)
#       cell_by_gene.csv
#       cell_metadata.csv
# Xenium keeps everything in one flat directory, so unlike the Xenium reader
# this detection has to look upward as well.
.merscope_detect_paths <- function(merscope_dir) {
    f <- function(...) {
        p <- file.path(merscope_dir, ...)
        if (file.exists(p)) normalizePath(p) else NULL
    }
    parent <- dirname(merscope_dir)
    fp <- function(...) {
        p <- file.path(parent, ...)
        if (file.exists(p)) normalizePath(p) else NULL
    }

    # transcripts: prefer parquet over the (much larger) csv
    tx_path <- f("detected_transcripts.parquet") %||% f("detected_transcripts.csv")

    # boundaries: parquet file, else the legacy per-FOV hdf5 directory
    cell_bound_path <- f("cell_boundaries.parquet") %||% f("cell_boundaries")

    paths <- list(
        tx_path          = tx_path,
        cell_bound_path  = cell_bound_path,
        expr_path        = f("cell_by_gene.csv"),
        cell_meta_path   = f("cell_metadata.csv"),
        cell_cat_path    = f("cell_categories.csv"),
        cell_numcat_path = f("cell_numeric_categories.csv"),
        experiment_path  = fp("experiment.json"),
        image_dir        = f("images")
    )

    if (is.null(paths$tx_path) && is.null(paths$expr_path)) {
        stop("[merscope] no `detected_transcripts.*` or `cell_by_gene.csv` ",
             "found in:\n  ", merscope_dir,
             "\nIs this a MERSCOPE region directory?", call. = FALSE)
    }

    # keep NULLs -- downstream closures check for them -- but replace with
    # NA_character_ so list2env() produces bindings that fail loudly rather
    # than silently resolving to a parent-scope object of the same name.
    lapply(paths, function(p) p %||% NA_character_)
}



# MODULAR ####



## transcripts ####

# Columns dropped from MERSCOPE transcripts by default. Each was checked
# against the full 151M-row table rather than assumed:
#
#   barcode_id     1:1 bijection with `gene` (900/900/900). Recoverable.
#   transcript_id  A function of `gene`, but only 816 distinct values across
#                  151,017,583 rows -- NOT a per-detection id. Recoverable.
#   x, y           FOV-local pixels; global_x = 0.106 * x + <fov origin>, scale
#                  identical in every FOV tested, residuals ~5e-4 um (float32
#                  rounding). ONE-WAY: per-FOV origins survive in no retained
#                  column, so x/y cannot be reconstructed. Dropped anyway --
#                  FOV-local pixel space exists to register against per-FOV
#                  images and these exports carry none. `dropcols =
#                  character(0)` keeps them.
#
# Motivation is memory, not disk: only x/y meaningfully shrink the store
# (~1.1 GB; high-entropy float32 barely compresses) while the two id columns
# dictionary-encode to almost nothing. Four fewer columns is four fewer
# buffered per worker during the tile write, which is what exhausts memory on a
# full section.
.MERSCOPE_TX_DROPCOLS <- c("x", "y", "barcode_id", "transcript_id")

# Patch a MERSCOPE transcript schema.
#
# The vendor writes the pandas index as a leading column with a ZERO-LENGTH
# name. Any dplyr verb on a dataset carrying it fails with
#   "attempt to use zero-length variable name"
# because the arrow dplyr backend builds a data mask keyed by column name.
# Rename the empty field(s) so the mask is constructible; the caller then
# drops them.
.merscope_tx_schema <- function(path) {
    sc <- arrow::schema(arrow::open_dataset(path))
    nms <- names(sc)
    unnamed <- !nzchar(nms)
    if (any(unnamed)) {
        nms[unnamed] <- paste0(".merscope_unnamed_", seq_len(sum(unnamed)))
    }
    list(
        schema = do.call(arrow::schema,
            lapply(seq_along(nms), function(i) arrow::field(nms[i], sc[[i]]$type))),
        dropcols = nms[unnamed]
    )
}

# Lazy transcript query with the schema patched and index columns removed.
.merscope_tx_query <- function(path, score_threshold = NULL, FOVs = NULL,
                               dropcols = character(0L)) {
    sch <- .merscope_tx_schema(path)
    a <- arrow::open_dataset(sources = path, schema = sch$schema)
    drop <- unique(c(sch$dropcols, as.character(dropcols)))
    if (length(drop) > 0L) {
        a <- dplyr::select(a, -dplyr::any_of(drop))
    }
    if (!is.null(score_threshold)) {
        a <- dplyr::filter(a, transcript_score >= score_threshold)
    }
    if (!is.null(FOVs)) {
        fovs_ <- as.integer(FOVs)
        a <- dplyr::filter(a, fov %in% fovs_)
    }
    # MERSCOPE writes cell_id as int64, which parquetStore warns about and
    # which would not join cleanly against poly_ID (character, from EntityID).
    # Cast at read time so both sides agree.
    if ("cell_id" %in% names(a)) {
        a <- dplyr::mutate(a, cell_id = arrow::cast(cell_id, arrow::string()))
    }
    a
}

# Midpoint of the transcript coordinate range -- the reflection origin Giotto
# uses for flip_axis, computed lazily (min/max scan, no materialization).
.merscope_flip_anchor <- function(path, FOVs = NULL, verbose = NULL) {
    a <- .merscope_tx_query(path, FOVs = FOVs)
    e <- .dplyr_ext(a, sdimx = "global_x", sdimy = "global_y")
    ev <- .ext_to_num_vec(e)
    c(x0 = (ev[[1L]] + ev[[2L]]) / 2, y0 = (ev[[3L]] + ev[[4L]]) / 2)
}

# Apply Giotto's flip_axis semantics to a geom store via the lazy affine layer.
# Recorded in @post_ops and applied at materialization, so no data is rewritten.
.merscope_apply_flip <- function(store, flip_axis, flip_anchor) {
    if (identical(flip_axis, "none") || is.null(flip_anchor)) return(store)
    x0 <- flip_anchor[["x0"]]
    y0 <- flip_anchor[["y0"]]
    if (flip_axis %in% c("y", "both")) {
        store <- GiottoClass::flip(store, direction = "vertical", x0 = x0, y0 = y0)
    }
    if (flip_axis %in% c("x", "both")) {
        store <- GiottoClass::flip(store, direction = "horizontal", x0 = x0, y0 = y0)
    }
    store
}

# Disk-backed MERSCOPE transcript ingestion.
#
# Streams the source parquet through arrow with the index-column patch above,
# plus lazy score and FOV filters, then writes directly to a
# `parquetGeomTileStore`. As with Xenium transcripts, no `parquetStore`
# intermediate is needed -- points carry no row-ordering requirement.
.merscope_transcript_disk <- function(
    path,
    gsource,
    feat_type = c("rna", "Blank"),
    split_keyword = list("Blank"),
    FOVs = NULL,
    dropcols = .MERSCOPE_TX_DROPCOLS,
    score_threshold = NULL,
    flip_axis = "none",
    flip_anchor = NULL,
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

    # capture bare scalars only -- read_fun must be safe to ship to workers
    score_thr <- score_threshold
    fovs_     <- if (is.null(FOVs)) NULL else as.integer(FOVs)
    drop      <- as.character(dropcols)

    read_fun <- function(x, schema = NULL) {
        .merscope_tx_query(x, score_threshold = score_thr, FOVs = fovs_,
                           dropcols = drop)
    }

    fs <- fileStore(path = path, read_fun = read_fun)
    qs <- as(fs, "queryableStore")

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

    tx_store <- .merscope_apply_flip(tx_store, flip_axis, flip_anchor)

    if (output == "store") return(tx_store)

    gpoints <- GiottoClass::createGiottoPoints(
        x = tx_store,
        feat_type = feat_type,
        split_keyword = split_keyword
    )
    if (!inherits(gpoints, "list")) gpoints <- list(gpoints)
    gpoints
}



## polygons ####

# Determine the z-plane architecture of a MERSCOPE boundary table.
#
# Ports the v2 auto-detection from Giotto's createGiottoMerscopeObject():
# when every z-plane carries the same cell IDs *and* the same geometry, the
# export is a single 2D segmentation replicated across planes, and the
# redundant planes are dropped.
#
# Giotto compares the first 100 vertex coordinates between two planes. Here the
# geometries are stored as WKB, so the equivalent (and stricter) test is
# byte-equality of the serialized geometry for a sample of cells -- identical
# bytes imply identical vertices, and it avoids decoding anything.
#
# Returns a list(kind = "2d"|"replicated_2d"|"3d", z_indices, use_z).
#
# Why 3D is refused rather than reduced to one plane, and why the verdict rests
# on a sample of `n_sample` cells rather than all of them: adr/0010.
.merscope_zplane_architecture <- function(path, poly_z_indices = NULL,
                                          n_sample = 200L, verbose = NULL) {
    ds <- arrow::open_dataset(path)
    if (!"ZIndex" %in% names(ds)) {
        return(list(kind = "2d", z_indices = NA_integer_, use_z = NA_integer_))
    }

    z_all <- sort(ds |> dplyr::distinct(ZIndex) |> dplyr::collect() |>
                  dplyr::pull(ZIndex))
    if (!is.null(poly_z_indices)) {
        keep <- intersect(as.integer(poly_z_indices), as.integer(z_all))
        if (!length(keep)) {
            stop("[merscope] none of `poly_z_indices` (",
                 paste(poly_z_indices, collapse = ", "),
                 ") are present. Available: ", paste(z_all, collapse = ", "),
                 call. = FALSE)
        }
        z_all <- keep
    }

    if (length(z_all) == 1L) {
        GiottoUtils::vmsg(.v = verbose,
            "[merscope] single z-plane; 2D architecture")
        return(list(kind = "2d", z_indices = z_all, use_z = z_all[[1L]]))
    }

    zA <- z_all[[1L]]; zB <- z_all[[2L]]
    ids_A <- ds |> dplyr::filter(ZIndex == zA) |> dplyr::select(EntityID) |>
             dplyr::collect() |> dplyr::pull(EntityID)
    ids_B <- ds |> dplyr::filter(ZIndex == zB) |> dplyr::select(EntityID) |>
             dplyr::collect() |> dplyr::pull(EntityID)

    if (!setequal(ids_A, ids_B)) {
        return(.merscope_stop_3d(z_all, reason = "cell ID sets differ between z-planes"))
    }

    # sample cells and compare serialized geometry across all planes
    set.seed(1L)
    samp <- if (length(ids_A) > n_sample) sample(ids_A, n_sample) else ids_A
    g <- ds |> dplyr::filter(EntityID %in% samp) |>
         dplyr::select(EntityID, ZIndex, Geometry) |> dplyr::collect()
    sp <- split(g, g$EntityID)
    same <- vapply(sp, function(gi) {
        gi <- gi[order(gi$ZIndex), ]
        ref <- as.raw(gi$Geometry[[1L]])
        all(vapply(seq_len(nrow(gi)),
                   function(i) identical(as.raw(gi$Geometry[[i]]), ref),
                   logical(1L)))
    }, logical(1L))

    if (!all(same)) {
        return(.merscope_stop_3d(z_all, reason = sprintf(
            "geometry differs across z-planes for %d of %d sampled cells",
            sum(!same), length(same))))
    }

    GiottoUtils::vmsg(.v = verbose, sprintf(
        "[merscope] [Auto-Detection] Identical IDs AND vertices across %d z-planes (Replicated 2D Data).\n[merscope] [Auto-Correction] Using z-plane %s, dropping %d redundant plane(s).",
        length(z_all), zA, length(z_all) - 1L))
    list(kind = "replicated_2d", z_indices = z_all, use_z = zA)
}

# True 3D is out of scope for this reader. aggregateStacks() -- which combines
# per-plane spat_units into one aggregate unit -- has no GiottoDisk
# equivalent, so there is no correct way to reduce genuinely differing planes
# here. Fail loudly rather than ingest one arbitrary plane and hand back an
# object that silently discards most of the segmentation.
.merscope_stop_3d <- function(z_indices, reason) {
    stop(GiottoUtils::wrap_txt(sprintf(
        "[merscope] true 3D segmentation detected (%s).

This reader currently supports 2D and replicated-2D MERSCOPE exports only.
Genuine 3D requires aggregateStacks(), which has no GiottoDisk equivalent yet:
aggregating expression across per-plane spat_units is a store-level operation
that has not been implemented.

Available z-indices: %s

To proceed with a single plane anyway -- discarding the other planes, which is
NOT equivalent to Giotto's 3D handling -- pass an explicit single index, e.g.
  poly_z_indices = %s",
        reason, paste(z_indices, collapse = ", "), z_indices[[1L]]
    )), call. = FALSE)
}

# Disk-backed MERSCOPE polygon ingestion.
#
# MERSCOPE boundaries are already WKB MultiPolygon, the same representation
# GiottoDisk writes. Xenium's two-stage `parquetStore` -> `parquetGeomTile`
# route exists only to lock a `row_index` so ring winding survives arrow's
# threaded scan; here one row is one complete geometry, so ordering is
# meaningless.
#
# Route: WKB -> SpatVector -> parquetGeom -> parquetGeomTile. The decode/
# re-encode round trip is deliberate -- it reuses the existing SpatVector write
# path and costs ~75s for ~700k cells, below the aggregate step's noise floor.
# A byte-preserving streaming passthrough would need a new storeWrite method
# plus sedonadb derivation of x_index/y_index and max_poly_radius; not worth it
# at this scale.
.merscope_poly_disk <- function(
    path,
    gsource,
    name = "cell",
    poly_z_indices = NULL,
    polygon_format = c("parquet", "hdf5"),
    FOVs = NULL,
    metadata_path = NULL,
    threshold = 25000L,
    flip_axis = "none",
    flip_anchor = NULL,
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

    if (polygon_format == "hdf5") {
        stop(GiottoUtils::wrap_txt(
            "[merscope] `polygon_format = 'hdf5'` is not implemented yet.

Legacy MERSCOPE exports store boundaries as per-FOV HDF5
(cell_boundaries/feature_data_*.hdf5) rather than a single
cell_boundaries.parquet. Only the parquet path is supported for now."
        ), call. = FALSE)
    }
    checkmate::assert_file_exists(path)

    GiottoUtils::vmsg(
        sprintf("[merscope_poly_disk] loading boundary '%s'", name),
        .v = verbose
    )

    arch <- .merscope_zplane_architecture(path,
        poly_z_indices = poly_z_indices, verbose = verbose)

    ds <- arrow::open_dataset(path)
    q <- ds
    if (!is.na(arch$use_z) && "ZIndex" %in% names(ds)) {
        use_z <- arch$use_z
        q <- dplyr::filter(q, ZIndex == use_z)
    }

    # FOV subsetting. The boundary table carries no `fov` column -- only
    # cell_metadata.csv maps EntityID to FOV -- so the requested FOVs have to
    # be resolved to an id set first. Without this, `FOVs` would silently
    # subset transcripts while loading every polygon in the section.
    if (!is.null(FOVs)) {
        if (is.null(metadata_path) || is.na(metadata_path)) {
            stop("[merscope_poly_disk] `FOVs` requires cell_metadata.csv to ",
                 "map EntityID -> fov, but no metadata path was found.",
                 call. = FALSE)
        }
        cm <- data.table::fread(metadata_path, select = c("EntityID", "fov"))
        keep_ids <- cm[["EntityID"]][cm[["fov"]] %in% as.integer(FOVs)]
        if (!length(keep_ids)) {
            stop("[merscope_poly_disk] no cells found in FOV(s): ",
                 paste(FOVs, collapse = ", "), call. = FALSE)
        }
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "[merscope_poly_disk] FOV subset -> %s cells",
            format(length(keep_ids), big.mark = ",")))
        q <- dplyr::filter(q, EntityID %in% keep_ids)
    }

    x <- q |> dplyr::select(EntityID, Geometry) |> dplyr::collect()

    if (nrow(x) == 0L) {
        GiottoUtils::vmsg(.v = verbose, "[merscope_poly_disk] no polygons found")
        return(NULL)
    }
    GiottoUtils::vmsg(.v = verbose,
        sprintf("[merscope_poly_disk] decoding %s geometries",
                format(nrow(x), big.mark = ",")))

    # terra accepts a list of raw vectors as WKB natively; as.list() strips
    # arrow's binary class wrapper. Same idiom as methods-spatRelate.R.
    sv <- terra::vect(as.list(x$Geometry))
    terra::values(sv) <- data.frame(
        poly_ID = as.character(x$EntityID),
        stringsAsFactors = FALSE
    )

        #
    # No `verbose` on the first call: sourceWrite(gDirSource, SpatVector)
    # forwards `...` all the way into arrow::write_dataset(), which rejects
    # unknown arguments. The queryableStore path Xenium uses absorbs `verbose`
    # before that point, so the two are not interchangeable here.
    pg <- sourceWrite(gsource, sv, store_type = "parquetGeom")

    # Explicit `threshold`, because .auto_threshold() would guess in the wrong
    # unit here. Its polygon rule is 500k *vertex rows* per tile, which suits
    # the Xenium model where a geom store holds one row per vertex (162k cells
    # -> 4.05M vertex rows -> ~8 tiles). Our WKB rows are whole geometries, so
    # 703,879 cells is 703,879 rows and the same rule yields 2 tiles -- 4x
    # Xenium's cell count in a quarter of the partitions.
    #
    # That under-partitioning costs little at ingest but matters later:
    # calculateOverlap()/aggregateFeatures() fan out over tiles, so 2 tiles
    # caps a 700k-cell overlap at 2-way parallelism. 25k geometries/tile gives
    # ~28 tiles here, comparable to Xenium's cells-per-tile density.
    poly_store <- sourceWrite(gsource, pg,
        store_type = "parquetGeomTile", threshold = threshold,
        verbose = verbose, ...)

    poly_store <- .merscope_apply_flip(poly_store, flip_axis, flip_anchor)

    if (output == "store") return(poly_store)

    # Centroids are not pre-computed here; setGiotto(centroids_to_spatlocs =
    # TRUE) derives them lazily from the parquetGeomBase at attach time.
    GiottoClass::createGiottoPolygon(x = poly_store, name = name)
}



## expression ####

# MERSCOPE negative controls are identified by feature-name prefix. Unlike the
# 10x formats there is no features.tsv column 3 (or /features/feature_type in
# h5) carrying an authoritative class, so classification is by keyword.
.merscope_feat_classes <- function(feat_ids, split_keyword = "Blank") {
    kw <- as.character(unlist(split_keyword))
    out <- rep("rna", length(feat_ids))
    for (k in kw) {
        hit <- grepl(k, feat_ids, fixed = TRUE)
        out[hit] <- k
    }
    out
}

# Disk-backed MERSCOPE expression ingestion.
#
# `cell_by_gene.csv` is a dense wide cell x gene matrix, so it routes through
# `csvWideInput` rather than the 10x sparse `mtxInput`/`tenxH5Input` used by the
# Xenium reader. Post-processing otherwise matches: drop zero-detection
# features, split by feature class, name the classes to Giotto conventions.
.merscope_expression_disk <- function(
    path,
    gsource,
    split_by_type = TRUE,
    split_keyword = "Blank",
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
                                               split_keyword = split_keyword)

    pe <- sourceWrite(gsource, inp, store_type = "parquetExpr",
                      verbose = verbose, ...)

    # streaming Arrow distinct on col_id to find detected features, then
    # subset the store positionally (same approach as the Xenium reader)
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

    # Build exprObjs directly (bypass createExprObj / .evaluate_expr_matrix
    # since they don't recognize parquetExprStore).
    lapply(seq_along(store_list), function(i) {
        methods::new("exprObj",
            name       = "raw",
            exprMat    = store_list[[i]],
            spat_unit  = "cell",
            feat_type  = names(store_list)[[i]],
            provenance = "cell"
        )
    })
}



## cell metadata ####

# In-memory, matching the Xenium reader's treatment of cell metadata.
# center_x / center_y are the vendor centroids; they are retained as ordinary
# metadata columns since spatial locations come from the polygon store via
# setGiotto(centroids_to_spatlocs = TRUE).
.merscope_cellmeta <- function(path, FOVs = NULL, verbose = NULL, ...) {
    if (missing(path) || is.na(path)) {
        stop("[merscope_cellmeta] no cell metadata path provided", call. = FALSE)
    }
    checkmate::assert_file_exists(path)
    GiottoUtils::vmsg("[merscope_cellmeta] loading cell metadata", .v = verbose)

    dt <- data.table::fread(path)

    # Honour FOVs. Without this the object gets metadata for every cell in the
    # section while its polygons and points cover only the requested FOVs --
    # e.g. 703,879 metadata rows against 16,543 polygons -- and the two sides
    # disagree about which cells exist.
    if (!is.null(FOVs) && "fov" %in% names(dt)) {
        n_before <- nrow(dt)
        dt <- dt[dt[["fov"]] %in% as.integer(FOVs)]
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "[merscope_cellmeta] FOV subset -> %s of %s cells",
            format(nrow(dt), big.mark = ","), format(n_before, big.mark = ",")))
    }

    if ("EntityID" %in% names(dt)) {
        data.table::setnames(dt, "EntityID", "cell_ID")
        dt[, cell_ID := as.character(cell_ID)]
    }
    GiottoClass::createCellMetaObj(
        metadata = dt,
        spat_unit = "cell",
        feat_type = "rna",
        provenance = "cell"
    )
}
