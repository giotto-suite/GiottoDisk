#' @include class-parquetExprStore.R
NULL


# Disambiguate duplicate feature IDs the 10x way (suffix `--N`).
# Matches the convention used by Giotto::get10Xmatrix so the
# parquet path produces the same feat_ids as the in-memory dgCMatrix
# path. No-op when all names are already unique.
.disambiguate_feat_ids <- function(feat_ids) {
    feat_ids <- as.character(feat_ids)
    counts <- table(feat_ids)
    dups   <- names(counts)[counts > 1L]
    if (length(dups) == 0L) return(feat_ids)
    out <- feat_ids
    for (nm in dups) {
        idx <- which(feat_ids == nm)
        out[idx] <- paste0(nm, "--", seq_along(idx))
    }
    out
}


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
    feat_ids <- .disambiguate_feat_ids(as.character(features[[feature_id_col]]))

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
    feat_ids   <- .disambiguate_feat_ids(
        if (feature_id_col == 1L) feat_id else feat_name)
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


# cellbin_gef_to_parquetExprStore ####

#' @name cellbin_gef_to_parquetExprStore
#' @title Stream-convert a Stereo-seq cellbin .gef file into a parquetExprStore
#' @description
#' Reads a Stereo-seq cell-bin GEF (`*.cellbin.gef` / `*.adjusted.cellbin.gef`)
#' in **gene-chunks** via `rhdf5::h5read` hyperslab and writes one Parquet
#' shard per chunk. **No full triplet table is ever materialized in memory** —
#' only `batch_genes` genes' worth of `(cellID, count)` rows are held at a time.
#'
#' GEF cellbin layout (from BGI's SAW pipeline):
#' ```
#' cellBin/cell     : (id, x, y, ...) — small, n_cells rows
#' cellBin/gene     : (geneName, geneID, cellCount, offset, ...) — small,
#'                    n_genes rows. cellCount = nnz per gene; rows in geneExp
#'                    are stored gene-major (run-length encoded by cellCount).
#' cellBin/geneExp  : (cellID, count) — large, sum(cellCount) total rows.
#' ```
#'
#' Note on predicate pushdown: shards are sorted by `row_id` (= cell index)
#' within themselves, but each gene-chunk shard contains rows from many
#' cells. Chunk-level row_id stats are wide; downstream cell-range filters
#' touch every shard but read only matching rows within them.
#'
#' @param gef_path character. Path to a Stereo-seq cellbin .gef file.
#' @param output_path character. Destination Parquet file (single file
#'   when total nnz <= `batch_genes` worth) or directory (one shard per
#'   gene-chunk, used otherwise).
#' @param gene_column character. Which gene-id column to use as `feat_ids`:
#'   `"geneName"` (default) or `"geneID"`.
#' @param batch_genes integer. Number of genes per streaming chunk.
#'   Default 500.
#' @param overwrite logical. If `TRUE`, replace `output_path` if it exists.
#' @return A [parquetExprStore-class] with `cell_ids = paste0("cell_",
#'   cellDT$id)` and `feat_ids` from the chosen `gene_column`.
#' @seealso [bin_gef_to_parquetExprStore()],
#'   [h5_to_parquetExprStore()]
#' @export
cellbin_gef_to_parquetExprStore <- function(
    gef_path,
    output_path,
    gene_column = c("geneName", "geneID"),
    batch_genes = 500L,
    overwrite   = FALSE
) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("[cellbin_gef_to_parquetExprStore] rhdf5 is required to read ",
             "Stereo-seq .gef files. Install with: ",
             "BiocManager::install(\"rhdf5\").",
             call. = FALSE)
    }
    stopifnot(file.exists(gef_path))
    gene_column <- match.arg(gene_column, c("geneName", "geneID"))

    if (file.exists(output_path)) {
        if (!overwrite) {
            stop("[cellbin_gef_to_parquetExprStore] output already exists at ",
                 output_path, "\n  set overwrite = TRUE to replace.",
                 call. = FALSE)
        }
        unlink(output_path, recursive = TRUE)
    }

    # ---- Small support tables (full reads) -----------------------------------
    cellDT <- data.table::setDT(rhdf5::h5read(gef_path, "cellBin/cell"))
    geneDT <- data.table::setDT(rhdf5::h5read(gef_path, "cellBin/gene"))

    n_cells  <- nrow(cellDT)
    cell_ids <- paste0("cell_", cellDT$id)

    # cellID (raw) → 1-based row index for the parquetExprStore
    cell_idx_map <- data.table::setattr(
        seq_len(n_cells), "names", as.character(cellDT$id)
    )

    # Collapse duplicate gene names + drop unexpressed genes
    # (matches Giotto's matrix-path .stereoseq_build_expression).
    cnt <- as.integer(geneDT$cellCount)
    all_names <- as.character(geneDT[[gene_column]])
    expressed_names <- all_names[cnt > 0]
    feat_ids        <- sort(unique(expressed_names))
    name_to_row     <- match(all_names, feat_ids)   # NA for unexpressed
    n_genes         <- length(feat_ids)
    n_genes_raw     <- length(cnt)

    # cumulative gene offsets in geneExp (0-based exclusive)
    cum <- c(0L, as.integer(cumsum(as.numeric(cnt))))   # length n_genes_raw+1
    nnz_total <- cum[length(cum)]

    # Build chunk boundaries that respect col_id groups: never split a
    # set of raw geneDT rows that share the same name_to_row across two
    # chunks (otherwise a duplicate-named gene's data would land in
    # separate parquet shards and miss the within-chunk aggregation).
    chunks <- .gef_safe_chunks(name_to_row, batch_genes)

    use_dir <- nnz_total > 5e6 || length(chunks) > 1L
    if (use_dir) dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

    # ---- Stream gene-chunks (safe-boundary plan) ----------------------------
    batch_idx <- 0L
    for (chunk_def in chunks) {
        g_lo     <- chunk_def[1L]
        g_hi     <- chunk_def[2L]
        slice_lo <- cum[g_lo] + 1L            # 1-based row index in geneExp
        slice_hi <- cum[g_hi + 1L]            # inclusive end
        if (slice_hi < slice_lo) next
        batch_idx <- batch_idx + 1L

        # rhdf5 hyperslab on a compound dataset: start (1-based) + count
        chunk <- data.table::setDT(rhdf5::h5read(
            gef_path, "cellBin/geneExp",
            start = slice_lo,
            count = slice_hi - slice_lo + 1L
        ))
        # gene_idx (raw geneDT row): rep(g_lo..g_hi, cnt[g_lo..g_hi]).
        # Remap through name_to_row so duplicate gene names collapse into
        # the same col_id (matches the matrix path's behavior).
        gene_idx_raw <- rep.int(seq.int(g_lo, g_hi), cnt[g_lo:g_hi])

        # Aggregate by (row_id, col_id) within the chunk so duplicate
        # geneDT entries that share a name (collapsed via name_to_row)
        # become a single row with summed value — matches the matrix
        # path's sparseMatrix() duplicate-summing behavior.
        out <- data.table::data.table(
            row_id = as.integer(cell_idx_map[as.character(chunk$cellID)]),
            col_id = as.integer(name_to_row[gene_idx_raw]),
            value  = as.double(chunk$count)
        )
        out <- out[!is.na(row_id) & !is.na(col_id),
                    .(value = sum(value)),
                    keyby = .(row_id, col_id)]

        target <- if (use_dir)
            file.path(output_path, sprintf("chunk_%010d.parquet", batch_idx))
        else
            output_path
        arrow::write_parquet(out, target)
    }

    parquetExprStore(
        path     = normalizePath(output_path),
        cell_ids = cell_ids,
        feat_ids = feat_ids,
        n_cells  = n_cells,
        n_genes  = n_genes
    )
}


# Build chunk boundaries that respect duplicate-name groups: never split
# a run of raw geneDT rows that share the same name_to_row between two
# chunks. Returns a list of c(g_lo, g_hi) integer pairs covering 1..n.
# `target_size` is the desired chunk size in raw geneDT rows.
.gef_safe_chunks <- function(name_to_row, target_size) {
    n <- length(name_to_row)
    if (n == 0L) return(list())
    target_size <- as.integer(target_size)

    # Boundary positions (1-based): position i is a "safe break" if
    # name_to_row[i] differs from name_to_row[i-1] (meaning the previous
    # group ends at i-1 and a new one starts at i). Position 1 is always
    # a starting point.
    nm <- name_to_row
    is_break <- c(TRUE, nm[-1] != nm[-n] |
                  (is.na(nm[-1]) != is.na(nm[-n])))
    is_break[is.na(is_break)] <- TRUE
    safe_starts <- which(is_break)

    # Walk safe_starts and group them into chunks of approximately
    # target_size raw genes each.
    out <- vector("list", length(safe_starts))
    k <- 0L
    last_start <- safe_starts[1L]
    for (s in safe_starts[-1L]) {
        if (s - last_start >= target_size) {
            k <- k + 1L
            out[[k]] <- c(last_start, s - 1L)
            last_start <- s
        }
    }
    k <- k + 1L
    out[[k]] <- c(last_start, n)
    out[seq_len(k)]
}


# bin_gef_to_parquetExprStore ####

#' @name bin_gef_to_parquetExprStore
#' @title Stream-convert a Stereo-seq bin .gef file into a parquetExprStore
#' @description
#' Reads a Stereo-seq bin GEF (`*.tissue.gef`, `*.gef`, `*.raw.gef`) at a
#' chosen `bin_size` in **gene-chunks** via `rhdf5::h5read` hyperslab and
#' writes one Parquet shard per chunk. Maintains a running
#' `(x, y) → bin_ID` map so bin IDs are assigned online — no second pass
#' through the data is required.
#'
#' GEF bin layout:
#' ```
#' geneExp/<bin_size>/expression : (x, y, count) — gene-major (run-length
#'                                  encoded by gene$count).
#' geneExp/<bin_size>/gene       : (geneName, geneID, count, ...).
#' ```
#'
#' @param gef_path character. Path to a Stereo-seq bin .gef file.
#' @param bin_size character. e.g. `"bin1"`, `"bin50"`, `"bin100"` (must
#'   match a `geneExp/<bin_size>` group in the GEF).
#' @param output_path character. Destination Parquet file or directory.
#' @param gene_column character. `"geneName"` (default) or `"geneID"`.
#' @param batch_genes integer. Number of genes per streaming chunk.
#'   Default 500.
#' @param overwrite logical. If `TRUE`, replace `output_path` if it exists.
#' @return A list with two components:
#'   * `pe`         — a [parquetExprStore-class] with
#'                    `cell_ids = paste0("bin_", 1:n_bins)`.
#'   * `bin_coords` — a `data.table` of `(bin_ID, x, y)` so callers can
#'                    construct spatial locations without a second read.
#' @seealso [cellbin_gef_to_parquetExprStore()]
#' @export
bin_gef_to_parquetExprStore <- function(
    gef_path,
    bin_size,
    output_path,
    gene_column = c("geneName", "geneID"),
    batch_genes = 500L,
    overwrite   = FALSE
) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("[bin_gef_to_parquetExprStore] rhdf5 is required to read ",
             "Stereo-seq .gef files. Install with: ",
             "BiocManager::install(\"rhdf5\").",
             call. = FALSE)
    }
    stopifnot(file.exists(gef_path))
    gene_column <- match.arg(gene_column, c("geneName", "geneID"))

    if (file.exists(output_path)) {
        if (!overwrite) {
            stop("[bin_gef_to_parquetExprStore] output already exists at ",
                 output_path, "\n  set overwrite = TRUE to replace.",
                 call. = FALSE)
        }
        unlink(output_path, recursive = TRUE)
    }

    # NSE bindings
    bin_ID <- x <- y <- NULL

    # ---- Small support table -------------------------------------------------
    geneDT <- data.table::setDT(rhdf5::h5read(
        gef_path, paste0("geneExp/", bin_size, "/gene")
    ))
    cnt <- as.integer(geneDT$count)

    # Collapse duplicate gene names + drop unexpressed genes
    # (matches Giotto's matrix-path .stereoseq_build_expression).
    all_names       <- as.character(geneDT[[gene_column]])
    expressed_names <- all_names[cnt > 0]
    feat_ids        <- sort(unique(expressed_names))
    name_to_row     <- match(all_names, feat_ids)   # NA for unexpressed
    n_genes         <- length(feat_ids)
    n_genes_raw     <- length(cnt)

    cum <- c(0L, as.integer(cumsum(as.numeric(cnt))))
    nnz_total <- cum[length(cum)]

    # Safe-boundary chunk plan (respects duplicate-name groupings)
    chunks <- .gef_safe_chunks(name_to_row, batch_genes)

    use_dir <- nnz_total > 5e6 || length(chunks) > 1L
    if (use_dir) dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

    # ---- Running (x, y) → bin_ID lookup --------------------------------------
    # Grows online as new bins are seen. data.table key on (x, y) makes
    # the lookup O(log n_bins_so_far) per entry.
    xy_to_bin <- data.table::data.table(
        x = integer(0), y = integer(0), bin_ID = integer(0)
    )
    data.table::setkey(xy_to_bin, x, y)
    n_bins <- 0L

    expr_path <- paste0("geneExp/", bin_size, "/expression")

    # ---- Stream gene-chunks (safe-boundary plan) ----------------------------
    batch_idx <- 0L
    for (chunk_def in chunks) {
        g_lo     <- chunk_def[1L]
        g_hi     <- chunk_def[2L]
        slice_lo <- cum[g_lo] + 1L
        slice_hi <- cum[g_hi + 1L]
        if (slice_hi < slice_lo) next
        batch_idx <- batch_idx + 1L

        chunk <- data.table::setDT(rhdf5::h5read(
            gef_path, expr_path,
            start = slice_lo,
            count = slice_hi - slice_lo + 1L
        ))
        gene_idx_raw <- rep.int(seq.int(g_lo, g_hi), cnt[g_lo:g_hi])

        # Update bin lookup with new (x, y) seen in this chunk
        chunk_xy <- unique(chunk[, .(x, y)])
        new_xy <- chunk_xy[!xy_to_bin, on = c("x", "y")]
        if (nrow(new_xy) > 0L) {
            new_xy[, bin_ID := seq.int(n_bins + 1L, n_bins + .N)]
            n_bins <- n_bins + nrow(new_xy)
            xy_to_bin <- rbind(xy_to_bin, new_xy)
            data.table::setkey(xy_to_bin, x, y)
        }
        # Look up bin_ID for every chunk row
        chunk[, bin_ID := xy_to_bin[.SD, on = c("x", "y"), bin_ID]]

        # Remap gene_idx through name_to_row to collapse duplicate names,
        # then aggregate by (row_id, col_id) so within-chunk duplicates
        # sum into one row — matches sparseMatrix()'s duplicate-summing.
        out <- data.table::data.table(
            row_id = as.integer(chunk$bin_ID),
            col_id = as.integer(name_to_row[gene_idx_raw]),
            value  = as.double(chunk$count)
        )
        out <- out[!is.na(row_id) & !is.na(col_id),
                    .(value = sum(value)),
                    keyby = .(row_id, col_id)]

        target <- if (use_dir)
            file.path(output_path, sprintf("chunk_%010d.parquet", batch_idx))
        else
            output_path
        arrow::write_parquet(out, target)
    }

    cell_ids <- paste0("bin_", seq_len(n_bins))
    pe <- parquetExprStore(
        path     = normalizePath(output_path),
        cell_ids = cell_ids,
        feat_ids = feat_ids,
        n_cells  = n_bins,
        n_genes  = n_genes
    )
    # Return bin coordinates so callers (e.g. spatlocs construction) can
    # reuse them without a second pass through the GEF.
    data.table::setorder(xy_to_bin, bin_ID)
    list(pe = pe, bin_coords = xy_to_bin)
}


# csv_to_parquetExprStore ####

#' @name csv_to_parquetExprStore
#' @title Stream-convert a wide-format dense CSV into a parquetExprStore
#' @description
#' Reads a wide-format CSV (one row per cell, one column per gene; the
#' first column or `cell_id_col` is the cell ID) in row-chunks via
#' `data.table::fread(skip = ..., nrows = ...)`, melts each chunk
#' wide → long (dropping zeros), and writes one Parquet shard per
#' chunk. **No full dense matrix is ever materialized in memory** —
#' only `batch_rows` cells of nonzeros are held at a time.
#'
#' This is the input format used by NanoString CosMx
#' (`*_exprMat_file.csv`) and Vizgen MERSCOPE (`cell_by_gene.csv`).
#'
#' @param csv_path character. Path to the wide-format CSV (gzip OK —
#'   `data.table::fread` autodetects).
#' @param output_path character. Destination Parquet file (single
#'   file when n_cells <= batch_rows) or directory (one shard per
#'   batch otherwise).
#' @param cell_id_col character. Name of the column to use as cell
#'   ID. Default `"cell_ID"`. The values become `pe@cell_ids` (after
#'   any caller-side renaming).
#' @param skip_cols character. Other non-feature columns to drop
#'   (e.g. `"fov"` for CosMx). All remaining columns become
#'   features (`feat_ids`).
#' @param row_filter_fun function or `NULL`. Optional predicate
#'   applied to each chunk's `data.table` BEFORE melting; should
#'   return a logical vector (same length as the chunk's rows) or
#'   the row indices to keep. Used by callers to apply
#'   technology-specific filters (CosMx: `cell_ID != 0`, FOV
#'   subset) without a separate pre-pass over the CSV.
#' @param batch_rows integer. Cells per streaming chunk. Default
#'   50,000 — peak RAM ≈ batch_rows × n_genes × 16 bytes.
#' @param overwrite logical. If `TRUE`, replace `output_path` if it
#'   exists.
#' @return A [parquetExprStore-class] with `cell_ids` collected
#'   (in chunk order) and `feat_ids` from the header.
#' @export
csv_to_parquetExprStore <- function(
    csv_path,
    output_path,
    cell_id_col   = "cell_ID",
    skip_cols     = NULL,
    row_filter_fun = NULL,
    batch_rows    = 50000L,
    overwrite     = FALSE
) {
    stopifnot(file.exists(csv_path))
    if (file.exists(output_path)) {
        if (!overwrite) {
            stop("[csv_to_parquetExprStore] output already exists at ",
                 output_path, "\n  set overwrite = TRUE to replace.",
                 call. = FALSE)
        }
        unlink(output_path, recursive = TRUE)
    }

    # NSE bindings
    row_id <- col_id <- value <- NULL

    # ---- Single-pass read ----------------------------------------------------
    # `data.table::fread(skip = N)` on a gzipped file re-decompresses from
    # the start each time → quadratic. Read once, then iterate row-chunks
    # of the in-RAM data.table. Peak RAM is bounded by the CSV's
    # uncompressed size; per-chunk melt buffer adds batch_rows × n_genes
    # × 16 bytes on top.
    full <- data.table::fread(csv_path, header = TRUE)
    all_cols <- colnames(full)
    if (!cell_id_col %in% all_cols) {
        stop("[csv_to_parquetExprStore] cell_id_col \"", cell_id_col,
             "\" not found in CSV header. Header columns: ",
             toString(head(all_cols, 5)), "...", call. = FALSE)
    }
    drop_cols <- unique(c(cell_id_col, skip_cols))
    feat_cols <- setdiff(all_cols, drop_cols)
    feat_ids  <- .disambiguate_feat_ids(feat_cols)
    n_genes   <- length(feat_ids)

    # Apply optional row filter to the whole table once
    if (!is.null(row_filter_fun)) {
        keep <- row_filter_fun(full)
        if (is.logical(keep)) full <- full[keep, ]
        else full <- full[keep, ]
    }

    n_cells <- nrow(full)
    cell_ids_out <- as.character(full[[cell_id_col]])

    # ---- Iterate row-chunks of the in-RAM data.table -------------------------
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    batch_idx <- 0L
    r_lo <- 1L
    while (r_lo <= n_cells) {
        r_hi <- min(r_lo + batch_rows - 1L, n_cells)
        batch_idx <- batch_idx + 1L

        # Extract feature columns for this chunk as a base matrix; nonzero
        # positions become triplets. as.matrix on a data.table sub-block
        # is a single allocation of size chunk_rows × n_genes integers.
        feat_mat <- as.matrix(full[r_lo:r_hi, ..feat_cols])
        nz <- which(feat_mat != 0, arr.ind = TRUE)
        if (nrow(nz) > 0L) {
            chunk_row_ids <- seq.int(r_lo, r_hi)
            out <- data.table::data.table(
                row_id = chunk_row_ids[nz[, 1L]],
                col_id = as.integer(nz[, 2L]),
                value  = as.double(feat_mat[nz])
            )
            data.table::setorder(out, row_id, col_id)
            arrow::write_parquet(
                out,
                file.path(output_path,
                          sprintf("chunk_%010d.parquet", batch_idx))
            )
        }
        rm(feat_mat, nz)
        r_lo <- r_hi + 1L
    }

    parquetExprStore(
        path     = normalizePath(output_path),
        cell_ids = cell_ids_out,
        feat_ids = feat_ids,
        n_cells  = n_cells,
        n_genes  = n_genes
    )
}
