# Xenium/Atera zarr -> parquet conversion: the exported standalone
# converter and the fingerprint-keyed cache wrapper the disk readers use.
#
# Cache layout (layer 1, reader path):
#   <getArtifactDumpDir()>/zarr_cache/<fingerprint>/
#     manifest.json                  completed products + provenance
#     cells.parquet                  10x cellmeta schema
#     cell_boundaries.parquet        cell_id, vertex_x, vertex_y
#     nucleus_boundaries.parquet
#     transcripts.parquet | transcripts/   (file, or shard dir if parallel)
# One fingerprint covers one source archive; the cells archive converts
# all three of its products in a single pass, so load_polys("cell"),
# load_polys("nucleus") and load_cellmeta share one conversion.

# bump when converter output changes; part of the cache fingerprint so
# stale conversions are redone rather than reused
.zarr_converter_version <- "0.2.0"

# Cheap content fingerprint of a zarr source (archive file or unzipped
# tree): path + size + mtime + layout/converter versions. Not a content
# hash -- matches the "size + mtime + partial identity" cache contract.
.zarr_fingerprint <- function(zarr_path) {
    if (dir.exists(zarr_path)) {
        files <- list.files(zarr_path, recursive = TRUE, full.names = TRUE)
        info <- file.info(files)
        size <- sum(info$size, na.rm = TRUE)
        mtime <- if (nrow(info)) max(as.numeric(info$mtime), na.rm = TRUE)
            else 0
    } else {
        info <- file.info(zarr_path)
        size <- as.numeric(info$size)
        mtime <- as.numeric(info$mtime)
    }
    .hash(list(
        path = normalizePath(zarr_path, winslash = "/", mustWork = FALSE),
        size = size,
        mtime = mtime,
        layout = .zarr_layout_version,
        converter = .zarr_converter_version
    ))
}

# Gene-identity lookup for transcripts: gene_panel.json beside the
# archive, else the CFM .zattrs feature catalog (ids align), else error.
.zarr_gene_lookup <- function(xenium_dir) {
    lookup <- .load_gene_panel_lookup(xenium_dir)
    if (length(lookup)) return(lookup)
    cfm_path <- .zarr_find_product(xenium_dir, "cell_feature_matrix")
    if (!is.null(cfm_path)) {
        src <- .zarr_open(cfm_path)
        on.exit(.zarr_close(src), add = TRUE)
        cat_info <- try(.load_cfm_feature_catalog(src), silent = TRUE)
        if (!inherits(cat_info, "try-error") &&
            length(cat_info$feature_ids)) {
            return(cat_info$feature_ids)
        }
    }
    stop("[zarr] no gene_panel.json in ", xenium_dir, " and no ",
        "cell_feature_matrix zarr to fall back to -- cannot resolve ",
        "transcript feature names", call. = FALSE)
}

# products of the cells archive (converted together in one pass)
.zarr_cells_products <- c(
    "cells", "cell_boundaries", "nucleus_boundaries"
)

# Convert-once wrapper used by the disk readers. Returns the converted
# parquet path (file, or shard dir for parallel transcripts) for `what`,
# converting only when the fingerprint-keyed cache misses.
.zarr_ensure_parquet <- function(zarr_path, what, workers = .par_workers(),
    compression = "zstd", verbose = NULL) {
    what <- match.arg(what, c("transcripts", .zarr_cells_products))
    cache_dir <- file.path(
        getArtifactDumpDir(), "zarr_cache", .zarr_fingerprint(zarr_path)
    )
    manifest_path <- file.path(cache_dir, "manifest.json")

    manifest <- if (file.exists(manifest_path)) {
        jsonlite::fromJSON(manifest_path, simplifyVector = TRUE)
    } else {
        NULL
    }
    hit <- !is.null(manifest) && what %in% names(manifest$products)
    if (hit) {
        target <- file.path(cache_dir, manifest$products[[what]])
        if (file.exists(target) || dir.exists(target)) {
            GiottoUtils::vmsg(
                "[zarr] using cached conversion:", target, .v = verbose
            )
            return(target)
        }
    }

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    products <- as.list(manifest$products)

    if (what %in% .zarr_cells_products) {
        src <- .zarr_open(zarr_path)
        on.exit(.zarr_close(src), add = TRUE)
        GiottoUtils::vmsg("[zarr] converting cells archive:", zarr_path,
            .v = verbose)
        cid <- .zarr_array(src, "cell_id")
        .zarr_cells_to_parquet(
            src, file.path(cache_dir, "cells.parquet"),
            compression = compression, verbose = verbose
        )
        .zarr_boundaries_to_parquet(
            src,
            out_cell = file.path(cache_dir, "cell_boundaries.parquet"),
            out_nuc = file.path(cache_dir, "nucleus_boundaries.parquet"),
            cell_id_lookup = list(prefix = cid[, 1L], suffix = cid[, 2L]),
            compression = compression, verbose = verbose
        )
        products[["cells"]] <- "cells.parquet"
        products[["cell_boundaries"]] <- "cell_boundaries.parquet"
        products[["nucleus_boundaries"]] <- "nucleus_boundaries.parquet"
    } else {
        src <- .zarr_open(zarr_path)
        on.exit(.zarr_close(src), add = TRUE)
        GiottoUtils::vmsg("[zarr] converting transcripts archive:",
            zarr_path, .v = verbose)
        layout <- .zarr_detect_product(zarr_path, "transcripts")
        gene_lookup <- if (identical(layout, "atera-zarr-v6")) NULL
            else .zarr_gene_lookup(dirname(zarr_path))
        out_rel <- if (workers > 1L) "transcripts" else "transcripts.parquet"
        .zarr_transcripts_to_parquet(
            src, file.path(cache_dir, out_rel),
            gene_lookup = gene_lookup,
            qv_threshold = NULL, # never filter in the cache; qv is lazy
            workers = workers, zarr_path = zarr_path,
            compression = compression, verbose = verbose, layout = layout
        )
        products[["transcripts"]] <- out_rel
    }

    .json_atomic_write(
        list(
            source = normalizePath(zarr_path, winslash = "/"),
            fingerprint = basename(cache_dir),
            layout_version = .zarr_layout_version,
            converter_version = .zarr_converter_version,
            created = .timestamp(),
            products = products
        ),
        file = manifest_path
    )
    file.path(cache_dir, products[[what]])
}

# Resolve an output path, refusing to replace an existing one unless
# `overwrite`.
.zarr_out_path <- function(p, overwrite, caller) {
    if (file.exists(p) || dir.exists(p)) {
        if (!overwrite) {
            stop("[", caller, "] output exists: ", p,
                "\n  pass overwrite = TRUE to replace.", call. = FALSE)
        }
        unlink(p, recursive = TRUE)
    }
    p
}

# cells.zarr -> cells / boundary parquet files in `output_dir`; one pass
# over the archive produces all requested files.
.zarr_convert_cells <- function(path, output_dir,
    what = c("cells", "cell_boundaries", "nucleus_boundaries"),
    overwrite = FALSE, compression = "zstd", verbose = NULL) {
    what <- match.arg(what, several.ok = TRUE)
    .zarr_detect_product(path, "cells")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    # resolve outputs (which may refuse) BEFORE opening the archive
    out <- lapply(stats::setNames(nm = what), function(w) {
        .zarr_out_path(file.path(output_dir, paste0(w, ".parquet")),
            overwrite, ".zarr_convert_cells")
    })
    src <- .zarr_open(path)
    on.exit(.zarr_close(src), add = TRUE)
    if (!is.null(out$cells)) {
        .zarr_cells_to_parquet(src, out$cells, compression = compression,
            verbose = verbose)
    }
    if (!is.null(out$cell_boundaries) || !is.null(out$nucleus_boundaries)) {
        cid <- .zarr_array(src, "cell_id")
        .zarr_boundaries_to_parquet(
            src, out_cell = out$cell_boundaries,
            out_nuc = out$nucleus_boundaries,
            cell_id_lookup = list(prefix = cid[, 1L], suffix = cid[, 2L]),
            compression = compression, verbose = verbose
        )
    }
    invisible(out)
}

# transcripts.zarr -> 10x-schema transcripts parquet (file, or shard
# directory with workers > 1).
.zarr_convert_transcripts <- function(path, output, qv_threshold = NULL,
    workers = .par_workers(), window_rows = 4194304L, gene_lookup = NULL,
    overwrite = FALSE, compression = "zstd", verbose = NULL) {
    layout <- .zarr_detect_product(path, "transcripts")
    if (is.null(gene_lookup) && !identical(layout, "atera-zarr-v6")) {
        gene_lookup <- .zarr_gene_lookup(dirname(path))
    }
    output <- .zarr_out_path(output, overwrite, ".zarr_convert_transcripts")
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    src <- .zarr_open(path)
    on.exit(.zarr_close(src), add = TRUE)
    .zarr_transcripts_to_parquet(
        src, output, gene_lookup = gene_lookup,
        qv_threshold = qv_threshold, workers = workers, zarr_path = path,
        compression = compression, verbose = verbose, layout = layout,
        window_rows = as.integer(window_rows)
    )
    invisible(output)
}

# cell_feature_matrix zarr -> long triplet parquet (row_id, col_id, value;
# sorted by row_id), drained from tenxZarrInput (either layout).
.zarr_cfm_to_parquet <- function(path, out_path, compression = "zstd",
    verbose = NULL) {
    inp <- tenxZarrInput(path)
    itr <- storeRead(inp)
    on.exit(itr$close(), add = TRUE)
    schema_cfm <- arrow::schema(
        row_id = arrow::int32(),
        col_id = arrow::int32(),
        value = arrow::float64()
    )
    sink <- arrow::FileOutputStream$create(out_path)
    props <- arrow::ParquetWriterProperties$create(
        column_names = schema_cfm$names,
        compression = .zarr_codec(compression)
    )
    writer <- arrow::ParquetFileWriter$create(
        schema = schema_cfm, sink = sink, properties = props
    )
    total <- 0
    repeat {
        dt <- itr$next_batch()
        if (is.null(dt)) break
        if (nrow(dt) == 0L) next
        writer$WriteTable(arrow::arrow_table(dt, schema = schema_cfm),
            chunk_size = nrow(dt))
        total <- total + nrow(dt)
    }
    writer$Close()
    sink$close()
    GiottoUtils::vmsg(sprintf(
        "  cell_feature_matrix.parquet: %d rows", total
    ), .v = verbose)
    out_path
}

#' @name xeniumZarrToParquet
#' @title Convert a Xenium/Atera zarr output directory to parquet
#' @description
#' Standalone converter from the `.zarr.zip` archives shipped by 10x
#' Xenium and Atera to parquet files in the 10x-shipped schemas:
#' `cells.parquet`, `cell_boundaries.parquet` / `nucleus_boundaries.parquet`
#' (from `cells.zarr.zip`), `transcripts.parquet` (from
#' `transcripts.zarr.zip`), and `cell_feature_matrix.parquet` (from
#' `cell_feature_matrix.zarr.zip`, written in the long triplet layout
#' consumed by [parquetExprStore()]: `row_id`, `col_id`, `value`, sorted by
#' `row_id`). An Atera `experiment.spatial` is copied as `experiment.xenium`.
#'
#' Output is unfiltered and unflipped. Transcripts with `workers > 1` are
#' written as a `transcripts/` directory of parquet shards, which
#' [arrow::open_dataset()] reads transparently. Atera v1 transcripts keep
#' their per-transcript `cell_id` and `overlaps_nucleus`; `transcript_id`
#' is synthesized and `fov_name` is `NA`.
#'
#' The disk readers ([importXeniumDisk()], [importAteraDisk()]) read zarr
#' input directly -- this entry point is for producing parquet outside of
#' a Giotto pipeline.
#' @param xenium_dir Xenium/Atera output directory containing the
#'   `.zarr.zip` archives (unzipped `.zarr` directories also work)
#' @param output_dir directory to write into (created if needed)
#' @param what which products to convert
#' @param qv_threshold optional Phred quality filter for transcripts
#'   (default `NULL`: keep everything, matching the 10x-shipped file)
#' @param overwrite overwrite existing outputs (default `FALSE`)
#' @param workers parallel workers for the transcripts step (default
#'   taken from `options("giottodisk.par_workers")`, or the future plan capped
#'   by the memory budget)
#' @param compression parquet compression codec (falls back to
#'   `"snappy"` when the arrow build lacks the requested codec)
#' @param verbose be verbose
#' @returns invisible named list of output paths for the converted
#'   products
#' @examples
#' \dontrun{
#' xeniumZarrToParquet("/path/to/xenium_output", "/path/to/parquet_out")
#' }
#' @export
xeniumZarrToParquet <- function(
    xenium_dir,
    output_dir,
    what = c("transcripts", "cells", "cell_boundaries",
        "nucleus_boundaries", "cell_feature_matrix"),
    qv_threshold = NULL,
    overwrite = FALSE,
    workers = .par_workers(),
    compression = "zstd",
    verbose = NULL) {
    checkmate::assert_directory_exists(xenium_dir)
    what <- match.arg(what, several.ok = TRUE)
    need <- c(
        cells = any(what %in% .zarr_cells_products),
        transcripts = "transcripts" %in% what,
        cell_feature_matrix = "cell_feature_matrix" %in% what
    )
    paths <- lapply(stats::setNames(nm = names(need)[need]), function(p) {
        zp <- .zarr_find_product(xenium_dir, p)
        if (is.null(zp)) {
            stop("[xeniumZarrToParquet] ", p, " zarr not found in ",
                xenium_dir, call. = FALSE)
        }
        zp
    })
    layouts <- lapply(stats::setNames(nm = names(paths)), function(p) {
        .zarr_detect_product(paths[[p]], p)
    })
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    out <- list()
    if (need[["cells"]]) {
        out <- c(out, .zarr_convert_cells(paths$cells, output_dir,
            what = intersect(what, .zarr_cells_products),
            overwrite = overwrite, compression = compression,
            verbose = verbose))
    }
    if (need[["transcripts"]]) {
        out$transcripts <- .zarr_convert_transcripts(
            paths$transcripts,
            file.path(output_dir,
                if (workers > 1L) "transcripts" else "transcripts.parquet"),
            qv_threshold = qv_threshold, workers = workers,
            overwrite = overwrite, compression = compression,
            verbose = verbose
        )
    }
    if (need[["cell_feature_matrix"]]) {
        p <- .zarr_out_path(
            file.path(output_dir, "cell_feature_matrix.parquet"),
            overwrite, "xeniumZarrToParquet"
        )
        out$cell_feature_matrix <- .zarr_cfm_to_parquet(
            paths$cell_feature_matrix, p, compression = compression,
            verbose = verbose
        )
    }

    exp_spatial <- file.path(xenium_dir, "experiment.spatial")
    exp_out <- file.path(output_dir, "experiment.xenium")
    if (file.exists(exp_spatial) && (overwrite || !file.exists(exp_out))) {
        file.copy(exp_spatial, exp_out, overwrite = TRUE)
    }

    ids <- unique(unlist(layouts))
    .json_atomic_write(
        list(
            source = normalizePath(xenium_dir, winslash = "/"),
            layout_version = if ("atera-zarr-v6" %in% ids) "atera-zarr-v6"
                else .zarr_layout_version,
            product_layouts = layouts,
            converter_version = .zarr_converter_version,
            created = .timestamp(),
            qv_threshold = qv_threshold,
            products = lapply(out, basename)
        ),
        file = file.path(output_dir, "conversion_manifest.json")
    )
    invisible(out)
}
