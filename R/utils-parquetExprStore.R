#' @include class-parquetExprStore.R
NULL

# mtx_to_parquetExprStore ####

#' @name mtx_to_parquetExprStore
#' @title Stream-convert a 10x / Xenium MatrixMarket triple into a parquetExprStore
#' @description
#' Reads a MatrixMarket sparse matrix (`.mtx` or `.mtx.gz`) line by line in
#' batches, swaps the (gene-row, cell-col) orientation to the streaming
#' convention (`row_id = cell`, `col_id = gene`), sorts each batch by
#' `row_id`, and writes one Parquet file per batch. **No dense or sparse
#' matrix is ever materialized in memory** — only `batch_lines` triplet
#' rows are held at a time.
#'
#' The 10x / Xenium MatrixMarket output is already column-major
#' (cell-sorted), so each batch is internally cell-sorted and Arrow's
#' file-level statistics enable predicate pushdown across the directory of
#' Parquet chunks.
#'
#' @param mtx_path character. Path to `matrix.mtx` or `matrix.mtx.gz`.
#'   The first two non-comment columns are interpreted as
#'   `(gene_index, cell_index)` per the 10x / Xenium convention.
#' @param barcodes_path character. Path to `barcodes.tsv` or
#'   `barcodes.tsv.gz` — one barcode per line; length must equal the
#'   number of columns in the mtx (n_cells).
#' @param features_path character. Path to `features.tsv` or
#'   `features.tsv.gz`. By default takes column 2 (gene symbol) as
#'   `feat_id`; pass `feature_id_col = 1L` to use the Ensembl id.
#' @param output_path character. Destination Parquet file (single file,
#'   used when `nnz <= batch_lines`) or directory (created if missing,
#'   one Parquet per batch).
#' @param feature_id_col integer. Which column of `features.tsv` to use as
#'   the gene identifier. Default `2L` (gene symbol).
#' @param batch_lines integer. Number of mtx triplet lines processed per
#'   batch. Default 5,000,000 — keeps peak RAM around ~120 MB per batch.
#' @param overwrite logical. If `TRUE`, replace `output_path` if it exists.
#' @return A [parquetExprStore-class] object pointing at `output_path` with
#'   `cell_ids` / `feat_ids` populated.
#' @seealso [xenium_to_parquetExprStore()] (convenience wrapper for the
#'   `cell_feature_matrix/` directory layout).
#' @export
mtx_to_parquetExprStore <- function(
    mtx_path,
    barcodes_path,
    features_path,
    output_path,
    feature_id_col = 2L,
    batch_lines    = 5000000L,
    overwrite      = FALSE
) {
    stopifnot(file.exists(mtx_path),
              file.exists(barcodes_path),
              file.exists(features_path))
    if (file.exists(output_path)) {
        if (!overwrite) {
            stop("[mtx_to_parquetExprStore] output already exists at ",
                 output_path, "\n  set overwrite = TRUE to replace.",
                 call. = FALSE)
        }
        unlink(output_path, recursive = TRUE)
    }

    # ---- Cell barcodes ----
    cell_ids <- readLines(barcodes_path)

    # ---- Gene features ----
    features <- data.table::fread(features_path, header = FALSE,
                                   sep = "\t", quote = "")
    if (feature_id_col > ncol(features)) {
        stop("[mtx_to_parquetExprStore] feature_id_col (", feature_id_col,
             ") exceeds number of columns in ", features_path,
             " (", ncol(features), ").", call. = FALSE)
    }
    feat_ids <- as.character(features[[feature_id_col]])

    # ---- Open mtx connection (handles both .gz and plain text) ----
    is_gz <- grepl("\\.gz$", mtx_path, ignore.case = TRUE)
    # Some 10x / Xenium files have .gz extension but are actually plain text;
    # detect by checking magic bytes 1f 8b.
    if (is_gz) {
        magic <- readBin(mtx_path, what = "raw", n = 2L)
        is_gz <- length(magic) == 2L && magic[1L] == as.raw(0x1f) &&
                  magic[2L] == as.raw(0x8b)
    }
    con <- if (is_gz) gzfile(mtx_path, "r") else file(mtx_path, "r")
    on.exit(close(con), add = TRUE)

    # ---- Skip header / comment lines, read dimensions ----
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (length(line) == 0L) {
            stop("[mtx_to_parquetExprStore] empty mtx file or no dimension line found.",
                 call. = FALSE)
        }
        if (!startsWith(line, "%")) break
    }
    dims <- as.integer(strsplit(trimws(line), "\\s+")[[1]])
    if (length(dims) != 3L || any(is.na(dims))) {
        stop("[mtx_to_parquetExprStore] could not parse dimension line: '",
             line, "'", call. = FALSE)
    }
    n_genes <- dims[1L]   # mtx rows  = features
    n_cells <- dims[2L]   # mtx cols  = cells
    nnz     <- dims[3L]

    if (length(cell_ids) != n_cells) {
        stop("[mtx_to_parquetExprStore] number of barcodes (",
             length(cell_ids), ") does not match mtx columns (",
             n_cells, ").", call. = FALSE)
    }
    if (length(feat_ids) != n_genes) {
        stop("[mtx_to_parquetExprStore] number of features (",
             length(feat_ids), ") does not match mtx rows (",
             n_genes, ").", call. = FALSE)
    }

    # ---- Decide single-file vs directory output ----
    use_dir <- nnz > batch_lines
    if (use_dir) {
        dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    }

    # ---- Stream-read in batches, swap orientation, sort, write ----
    lines_remaining <- nnz
    batch_idx       <- 0L
    while (lines_remaining > 0L) {
        n_to_read <- min(batch_lines, lines_remaining)
        raw_batch <- readLines(con, n = n_to_read, warn = FALSE)
        if (length(raw_batch) == 0L) break
        batch_idx <- batch_idx + 1L

        parsed <- data.table::fread(text = raw_batch, header = FALSE,
                                    sep = " ", quote = "",
                                    col.names = c("g", "c", "v"))

        # Swap orientation: row_id = cell, col_id = gene
        out <- data.table::data.table(
            row_id = as.integer(parsed$c),
            col_id = as.integer(parsed$g),
            value  = as.double(parsed$v)
        )
        # mtx is column-major (cell-sorted), so each batch is already grouped
        # by cell. Sort within batch to lock cell-order + within-cell gene order.
        data.table::setorder(out, row_id, col_id)

        if (use_dir) {
            file_name <- sprintf("chunk_%010d.parquet", batch_idx)
            arrow::write_parquet(out, file.path(output_path, file_name))
        } else {
            arrow::write_parquet(out, output_path)
        }

        lines_remaining <- lines_remaining - length(raw_batch)
    }

    parquetExprStore(
        path     = normalizePath(output_path),
        cell_ids = cell_ids,
        feat_ids = feat_ids,
        n_cells  = n_cells,
        n_genes  = n_genes
    )
}


# xenium_to_parquetExprStore ####

#' @name xenium_to_parquetExprStore
#' @title Convert a Xenium `cell_feature_matrix/` directory into a parquetExprStore
#' @description
#' Convenience wrapper around [mtx_to_parquetExprStore()] for the standard
#' Xenium output layout:
#' ```
#' cell_feature_matrix/
#'   ├── barcodes.tsv.gz
#'   ├── features.tsv.gz
#'   └── matrix.mtx.gz
#' ```
#' @param xenium_dir character. Path to a `cell_feature_matrix/` directory
#'   from a Xenium pipeline output.
#' @param output_path character. Destination Parquet file or directory.
#' @param ... passed to [mtx_to_parquetExprStore()] (e.g. `batch_lines`,
#'   `feature_id_col`, `overwrite`).
#' @return A [parquetExprStore-class] object.
#' @export
xenium_to_parquetExprStore <- function(xenium_dir, output_path, ...) {
    stopifnot(dir.exists(xenium_dir))
    mtx_to_parquetExprStore(
        mtx_path      = file.path(xenium_dir, "matrix.mtx.gz"),
        barcodes_path = file.path(xenium_dir, "barcodes.tsv.gz"),
        features_path = file.path(xenium_dir, "features.tsv.gz"),
        output_path   = output_path,
        ...
    )
}
