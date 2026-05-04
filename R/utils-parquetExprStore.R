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


# h5_to_parquetExprStore ####

#' @name h5_to_parquetExprStore
#' @title Stream-convert a 10x cell_feature_matrix.h5 into a parquetExprStore
#' @description
#' Reads a 10x HDF5 sparse-matrix file in cell-chunks (each chunk a
#' contiguous slice of CSC `data` / `indices` between two `indptr`
#' boundaries), swaps the (gene-row, cell-col) orientation to the
#' streaming convention (`row_id = cell`, `col_id = gene`), and writes
#' one Parquet file per chunk. **No dense or sparse matrix is ever
#' materialized in memory** — only `batch_cells` cells of nonzeros are
#' held at a time.
#'
#' The 10x H5 layout is well-documented:
#' ```
#' /<root>/data        : float, length nnz
#' /<root>/indices     : int,   length nnz   (0-based gene indices)
#' /<root>/indptr      : int,   length n_cells + 1   (CSC col offsets)
#' /<root>/shape       : int,   length 2     (n_genes, n_cells)
#' /<root>/barcodes    : char,  length n_cells
#' /<root>/features/id, /name, /feature_type
#' ```
#'
#' @param h5_path character. Path to `cell_feature_matrix.h5`.
#' @param output_path character. Destination Parquet file (single file,
#'   used when n_cells <= batch_cells) or directory (one Parquet per
#'   batch).
#' @param feature_id_col integer. `1L` for Ensembl id (`features/id`),
#'   `2L` for gene symbol (`features/name`). Default `2L`.
#' @param batch_cells integer. Number of cells processed per batch.
#'   Default 250,000 — keeps each batch's RAM around the per-cell mean
#'   nnz × batch_cells bytes (~tens of MB for typical Xenium panels).
#' @param overwrite logical. If `TRUE`, replace `output_path` if it exists.
#' @return A [parquetExprStore-class] object.
#' @seealso [xenium_h5_to_parquetExprStore()] (convenience wrapper).
#' @export
h5_to_parquetExprStore <- function(
    h5_path,
    output_path,
    feature_id_col = 2L,
    batch_cells    = 250000L,
    overwrite      = FALSE
) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("[h5_to_parquetExprStore] hdf5r is required to read 10x .h5 ",
             "files. Install with: install.packages(\"hdf5r\").",
             call. = FALSE)
    }
    stopifnot(file.exists(h5_path))
    if (file.exists(output_path)) {
        if (!overwrite) {
            stop("[h5_to_parquetExprStore] output already exists at ",
                 output_path, "\n  set overwrite = TRUE to replace.",
                 call. = FALSE)
        }
        unlink(output_path, recursive = TRUE)
    }

    h5 <- hdf5r::H5File$new(h5_path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)

    root <- names(h5)[1L]

    # ---- Read small metadata in full -----------------------------------------
    cell_ids   <- as.character(h5[[paste0(root, "/barcodes")]][])
    feat_id    <- as.character(h5[[paste0(root, "/features/id")]][])
    feat_name  <- as.character(h5[[paste0(root, "/features/name")]][])
    feat_ids   <- if (feature_id_col == 1L) feat_id else feat_name
    n_genes    <- as.integer(h5[[paste0(root, "/shape")]][])[1L]
    n_cells    <- as.integer(h5[[paste0(root, "/shape")]][])[2L]
    indptr     <- as.numeric(h5[[paste0(root, "/indptr")]][])

    if (length(cell_ids) != n_cells) {
        stop("[h5_to_parquetExprStore] barcode count (", length(cell_ids),
             ") != shape[2] (", n_cells, ").", call. = FALSE)
    }
    if (length(feat_ids) != n_genes) {
        stop("[h5_to_parquetExprStore] feature count (", length(feat_ids),
             ") != shape[1] (", n_genes, ").", call. = FALSE)
    }

    nnz_total <- as.numeric(indptr[length(indptr)])
    use_dir   <- nnz_total > 5e6 || n_cells > batch_cells
    if (use_dir) {
        dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    }

    data_ds    <- h5[[paste0(root, "/data")]]
    indices_ds <- h5[[paste0(root, "/indices")]]

    # ---- Stream cell-chunks -------------------------------------------------
    # CSC layout: cell c's nonzeros live at data[indptr[c]+1 : indptr[c+1]]
    # (R 1-based). hdf5r hyperslab reads accept `args = list(start:end)`.
    batch_idx <- 0L
    c_lo <- 1L
    while (c_lo <= n_cells) {
        c_hi <- min(c_lo + batch_cells - 1L, n_cells)
        # indptr is 0-based offsets; convert to R 1-based slice indices
        slice_lo <- as.numeric(indptr[c_lo]) + 1
        slice_hi <- as.numeric(indptr[c_hi + 1L])
        n_slice  <- slice_hi - slice_lo + 1
        if (n_slice <= 0L) {
            # all cells in this chunk are empty
            c_lo <- c_hi + 1L
            next
        }
        batch_idx <- batch_idx + 1L

        vals <- as.double(data_ds[seq.int(slice_lo, slice_hi)])
        gidx <- as.integer(indices_ds[seq.int(slice_lo, slice_hi)])

        # Reconstruct cell index per nnz from indptr deltas
        nnz_per_cell <- diff(indptr[c_lo:(c_hi + 1L)])
        cell_idx_local <- rep.int(seq.int(c_lo, c_hi), nnz_per_cell)

        out <- data.table::data.table(
            row_id = as.integer(cell_idx_local),
            col_id = gidx + 1L,             # h5 is 0-based; parquet is 1-based
            value  = vals
        )
        # CSC chunk is already cell-grouped; lock within-cell gene order too.
        data.table::setorder(out, row_id, col_id)

        if (use_dir) {
            file_name <- sprintf("chunk_%010d.parquet", batch_idx)
            arrow::write_parquet(out, file.path(output_path, file_name))
        } else {
            arrow::write_parquet(out, output_path)
        }

        c_lo <- c_hi + 1L
    }

    parquetExprStore(
        path     = normalizePath(output_path),
        cell_ids = cell_ids,
        feat_ids = feat_ids,
        n_cells  = n_cells,
        n_genes  = n_genes
    )
}


# xenium_h5_to_parquetExprStore ####

#' @name xenium_h5_to_parquetExprStore
#' @title Convert a Xenium `cell_feature_matrix.h5` into a parquetExprStore
#' @description
#' Convenience wrapper around [h5_to_parquetExprStore()] for the standard
#' Xenium output layout (`<xenium_dir>/cell_feature_matrix.h5`).
#' @param xenium_dir character. Path to a Xenium output directory.
#' @param output_path character. Destination Parquet file or directory.
#' @param ... passed to [h5_to_parquetExprStore()] (e.g. `batch_cells`,
#'   `feature_id_col`, `overwrite`).
#' @return A [parquetExprStore-class] object.
#' @export
xenium_h5_to_parquetExprStore <- function(xenium_dir, output_path, ...) {
    stopifnot(dir.exists(xenium_dir))
    h5_path <- file.path(xenium_dir, "cell_feature_matrix.h5")
    if (!file.exists(h5_path)) {
        stop("[xenium_h5_to_parquetExprStore] no cell_feature_matrix.h5 ",
             "found in: ", xenium_dir, call. = FALSE)
    }
    h5_to_parquetExprStore(
        h5_path     = h5_path,
        output_path = output_path,
        ...
    )
}
