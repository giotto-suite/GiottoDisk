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
.zarr_converter_version <- "0.1.0"

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
        gene_lookup <- .zarr_gene_lookup(dirname(zarr_path))
        out_rel <- if (workers > 1L) "transcripts" else "transcripts.parquet"
        .zarr_transcripts_to_parquet(
            src, file.path(cache_dir, out_rel),
            gene_lookup = gene_lookup,
            qv_threshold = NULL, # never filter in the cache; qv is lazy
            workers = workers, zarr_path = zarr_path,
            compression = compression, verbose = verbose
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

#' @name xeniumZarrToParquet
#' @title Convert Xenium/Atera zarr output to parquet
#' @description
#' Standalone converter from the `.zarr.zip` archives shipped by 10x
#' Xenium (and Atera, which ships zarr only) to parquet files matching
#' the 10x-shipped schemas: `cells.parquet`,
#' `cell_boundaries.parquet` / `nucleus_boundaries.parquet` (from
#' `cells.zarr.zip`), `transcripts.parquet` (from `transcripts.zarr.zip`),
#' and `cell_feature_matrix.parquet` (from `cell_feature_matrix.zarr.zip`,
#' written in the long triplet layout consumed by [parquetExprStore()]:
#' `row_id`, `col_id`, `value`, sorted by `row_id`).
#'
#' Output is unfiltered and unflipped; the Giotto readers apply qv
#' filtering and the vertical flip lazily. Transcripts with
#' `workers > 1` are written as a `transcripts/` directory of parquet
#' shards, which [arrow::open_dataset()] (and the Giotto readers) consume
#' transparently.
#'
#' The disk readers ([importXeniumDisk()], [importAteraDisk()]) call the
#' converter automatically for zarr input -- this entry point is for
#' producing parquet outside of a Giotto pipeline.
#' @param xenium_dir Xenium/Atera output directory containing the
#'   `.zarr.zip` archives (unzipped `.zarr` directories also work)
#' @param output_dir directory to write parquet into (created if needed)
#' @param what which products to convert
#' @param qv_threshold optional Phred quality filter for transcripts
#'   (default `NULL`: keep everything, matching the 10x-shipped file)
#' @param overwrite overwrite existing outputs (default `FALSE`)
#' @param workers parallel workers for the transcripts step (default
#'   taken from `options("giottodisk.par_workers")` / the future plan)
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
    layout <- .zarr_layout_detect(xenium_dir, validate = TRUE)
    if (is.na(layout$layout_version)) {
        stop("[xeniumZarrToParquet] no zarr products found in ",
            xenium_dir, call. = FALSE)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    out <- list()
    out_file <- function(rel) {
        p <- file.path(output_dir, rel)
        if (!overwrite && (file.exists(p) || dir.exists(p))) {
            stop("[xeniumZarrToParquet] output exists: ", p,
                "\n  pass overwrite = TRUE to replace.", call. = FALSE)
        }
        if (dir.exists(p)) {
            unlink(p, recursive = TRUE)
        } else if (file.exists(p)) {
            unlink(p)
        }
        p
    }

    if (any(what %in% .zarr_cells_products)) {
        zp <- layout$products$cells
        if (is.null(zp)) {
            stop("[xeniumZarrToParquet] cells zarr not found in ",
                xenium_dir, call. = FALSE)
        }
        # resolve outputs (which may refuse without overwrite) BEFORE
        # opening the archive
        p <- if ("cells" %in% what) out_file("cells.parquet") else NULL
        pc <- if ("cell_boundaries" %in% what) {
            out_file("cell_boundaries.parquet")
        } else {
            NULL
        }
        pn <- if ("nucleus_boundaries" %in% what) {
            out_file("nucleus_boundaries.parquet")
        } else {
            NULL
        }
        src <- .zarr_open(zp)
        on.exit(.zarr_close(src), add = TRUE)
        cid <- .zarr_array(src, "cell_id")
        if (!is.null(p)) {
            .zarr_cells_to_parquet(src, p, compression = compression,
                verbose = verbose)
            out$cells <- p
        }
        if (!is.null(pc) || !is.null(pn)) {
            .zarr_boundaries_to_parquet(
                src, out_cell = pc, out_nuc = pn,
                cell_id_lookup = list(
                    prefix = cid[, 1L], suffix = cid[, 2L]
                ),
                compression = compression, verbose = verbose
            )
            out$cell_boundaries <- pc
            out$nucleus_boundaries <- pn
        }
        .zarr_close(src)
    }

    if ("transcripts" %in% what) {
        zp <- layout$products$transcripts
        if (is.null(zp)) {
            stop("[xeniumZarrToParquet] transcripts zarr not found in ",
                xenium_dir, call. = FALSE)
        }
        p <- out_file(
            if (workers > 1L) "transcripts" else "transcripts.parquet"
        )
        src <- .zarr_open(zp)
        on.exit(.zarr_close(src), add = TRUE)
        .zarr_transcripts_to_parquet(
            src, p,
            gene_lookup = .zarr_gene_lookup(xenium_dir),
            qv_threshold = qv_threshold,
            workers = workers, zarr_path = zp,
            compression = compression, verbose = verbose
        )
        .zarr_close(src)
        out$transcripts <- p
    }

    if ("cell_feature_matrix" %in% what) {
        zp <- layout$products$cell_feature_matrix
        if (is.null(zp)) {
            stop("[xeniumZarrToParquet] cell_feature_matrix zarr not ",
                "found in ", xenium_dir, call. = FALSE)
        }
        p <- out_file("cell_feature_matrix.parquet")
        inp <- tenxZarrInput(zp)
        itr <- storeRead(inp)
        schema_cfm <- arrow::schema(
            row_id = arrow::int32(),
            col_id = arrow::int32(),
            value = arrow::float64()
        )
        sink <- arrow::FileOutputStream$create(p)
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
        itr$close()
        GiottoUtils::vmsg(sprintf(
            "  cell_feature_matrix.parquet: %d rows", total
        ), .v = verbose)
        out$cell_feature_matrix <- p
    }

    .json_atomic_write(
        list(
            source = normalizePath(xenium_dir, winslash = "/"),
            layout_version = layout$layout_version,
            converter_version = .zarr_converter_version,
            created = .timestamp(),
            qv_threshold = qv_threshold,
            products = lapply(out, basename)
        ),
        file = file.path(output_dir, "conversion_manifest.json")
    )
    invisible(out)
}
