#' @include class-fileInputs.R class-parquetExprStore.R
NULL

# dim / dimnames ####
# Matrix-like introspection for any exprInput subclass. Mirrors the
# parquetExprStore convention: rows = genes (features), cols = cells.
# For accumulating-metadata formats (binGefInput, csvWideInput), n_cells
# and cell_ids reflect zero / empty until storeRead() has been driven to
# completion.

#' @export
setMethod("nrow", "exprInput", function(x) x@n_genes)

#' @export
setMethod("ncol", "exprInput", function(x) x@n_cells)

#' @export
setMethod("dim", "exprInput", function(x) c(x@n_genes, x@n_cells))

#' @export
setMethod("dimnames", "exprInput",
    function(x) list(x@feat_ids, x@cell_ids)
)


# storeRead — batch iterators ####
# Each exprInput subclass returns a `list(next_batch, close)` from
# storeRead(). Calling next_batch() returns a data.table with
# (row_id, col_id, value) for the next batch, or NULL on EOF.
# close() releases any open connection or file handle. Iterator state
# (connection position, cursor) lives in the closure environment.

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("mtxInput"), function(store, ...) {
    mtx_path <- store@path
    n_cells  <- store@n_cells
    n_genes  <- store@n_genes
    batch    <- store@batch_lines

    con <- if (.is_real_gz(mtx_path)) gzfile(mtx_path, "r") else file(mtx_path, "r")

    # Skip header / comment lines and read the dimension line.
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (length(line) == 0L) {
            close(con)
            stop("[storeRead] empty mtx or no dimension line at ", mtx_path,
                 call. = FALSE)
        }
        if (!startsWith(line, "%")) break
    }
    dims <- as.integer(strsplit(trimws(line), "\\s+")[[1L]])
    if (length(dims) != 3L || any(is.na(dims))) {
        close(con)
        stop("[storeRead] could not parse mtx dimension line: '", line, "'",
             call. = FALSE)
    }
    if (dims[1L] != n_genes || dims[2L] != n_cells) {
        close(con)
        stop("[storeRead] mtx header (genes=", dims[1L], ", cells=", dims[2L],
             ") disagrees with sidecar metadata (genes=", n_genes,
             ", cells=", n_cells, ")", call. = FALSE)
    }
    lines_remaining <- dims[3L]
    closed <- FALSE

    close_fn <- function() {
        if (!closed) {
            try(close(con), silent = TRUE)
            closed <<- TRUE
        }
        invisible(NULL)
    }

    next_batch <- function() {
        if (closed || lines_remaining <= 0L) {
            close_fn()
            return(NULL)
        }
        n_to_read <- min(batch, lines_remaining)
        raw_batch <- readLines(con, n = n_to_read, warn = FALSE)
        if (length(raw_batch) == 0L) {
            close_fn()
            return(NULL)
        }
        parsed <- data.table::fread(text = raw_batch, header = FALSE,
                                    sep = " ", quote = "",
                                    col.names = c("g", "c", "v"))
        # mtx is gene-row / cell-col; parquetExprStore is cell-row / gene-col.
        out <- data.table::data.table(
            row_id = as.integer(parsed$c),
            col_id = as.integer(parsed$g),
            value  = as.double(parsed$v)
        )
        # mtx is column-major (cell-sorted); sort within batch to lock
        # cell + within-cell gene order.
        data.table::setorder(out, row_id, col_id)
        lines_remaining <<- lines_remaining - length(raw_batch)
        out
    }

    list(
        next_batch = next_batch,
        close      = close_fn,
        # Eager-known metadata — accessors just return store slots.
        cell_ids   = function() store@cell_ids,
        feat_ids   = function() store@feat_ids,
        n_cells    = function() store@n_cells,
        n_genes    = function() store@n_genes
    )
})


# storeRead — tenxH5Input ####

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("tenxH5Input"), function(store, ...) {
    h5      <- hdf5r::H5File$new(store@path, mode = "r")
    root    <- store@params$root
    indptr  <- as.numeric(h5[[paste0(root, "/indptr")]][])
    data_ds <- h5[[paste0(root, "/data")]]
    indices_ds <- h5[[paste0(root, "/indices")]]
    n_cells <- store@n_cells
    batch   <- store@batch_cells

    closed <- FALSE
    c_lo   <- 1L

    close_fn <- function() {
        if (!closed) {
            try(h5$close_all(), silent = TRUE)
            closed <<- TRUE
        }
        invisible(NULL)
    }

    next_batch <- function() {
        if (closed || c_lo > n_cells) { close_fn(); return(NULL) }
        # Advance past any all-empty cell runs at the start of this batch.
        repeat {
            if (c_lo > n_cells) { close_fn(); return(NULL) }
            c_hi     <- min(c_lo + batch - 1L, n_cells)
            slice_lo <- as.numeric(indptr[c_lo]) + 1
            slice_hi <- as.numeric(indptr[c_hi + 1L])
            if (slice_hi >= slice_lo) break
            c_lo <<- c_hi + 1L
        }
        vals <- as.double(data_ds[seq.int(slice_lo, slice_hi)])
        gidx <- as.integer(indices_ds[seq.int(slice_lo, slice_hi)])
        nnz_per_cell   <- diff(indptr[c_lo:(c_hi + 1L)])
        cell_idx_local <- rep.int(seq.int(c_lo, c_hi), nnz_per_cell)
        out <- data.table::data.table(
            row_id = as.integer(cell_idx_local),
            col_id = gidx + 1L,       # h5 indices are 0-based
            value  = vals
        )
        data.table::setorder(out, row_id, col_id)
        c_lo <<- c_hi + 1L
        out
    }

    list(
        next_batch = next_batch,
        close      = close_fn,
        cell_ids   = function() store@cell_ids,
        feat_ids   = function() store@feat_ids,
        n_cells    = function() store@n_cells,
        n_genes    = function() store@n_genes
    )
})


# storeRead — cellbinGefInput ####

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("cellbinGefInput"), function(store, ...) {
    gef_path    <- store@path
    cell_id_map <- store@params$cell_id_map
    cnt         <- store@params$gene_cnt
    name_to_row <- store@params$name_to_row
    cum         <- store@params$cum_offsets

    chunks <- .gef_safe_chunks(name_to_row, store@batch_genes)
    chunk_i <- 0L
    closed  <- FALSE

    close_fn <- function() { closed <<- TRUE; invisible(NULL) }

    next_batch <- function() {
        repeat {
            if (closed || chunk_i >= length(chunks)) {
                close_fn(); return(NULL)
            }
            chunk_i <<- chunk_i + 1L
            chunk_def <- chunks[[chunk_i]]
            g_lo <- chunk_def[1L]; g_hi <- chunk_def[2L]
            slice_lo <- cum[g_lo] + 1L
            slice_hi <- cum[g_hi + 1L]
            if (slice_hi >= slice_lo) break
        }
        chunk <- data.table::setDT(rhdf5::h5read(
            gef_path, "cellBin/geneExp",
            start = slice_lo,
            count = slice_hi - slice_lo + 1L
        ))
        gene_idx_raw <- rep.int(seq.int(g_lo, g_hi), cnt[g_lo:g_hi])
        out <- data.table::data.table(
            row_id = as.integer(cell_id_map[as.character(chunk$cellID)]),
            col_id = as.integer(name_to_row[gene_idx_raw]),
            value  = as.double(chunk$count)
        )
        # Aggregate duplicate-name collapses within the chunk.
        out <- out[!is.na(row_id) & !is.na(col_id),
                   .(value = sum(value)), keyby = .(row_id, col_id)]
        out
    }

    list(
        next_batch = next_batch,
        close      = close_fn,
        cell_ids   = function() store@cell_ids,
        feat_ids   = function() store@feat_ids,
        n_cells    = function() store@n_cells,
        n_genes    = function() store@n_genes
    )
})


# storeRead — binGefInput ####

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("binGefInput"), function(store, ...) {
    gef_path    <- store@path
    cnt         <- store@params$gene_cnt
    name_to_row <- store@params$name_to_row
    cum         <- store@params$cum_offsets
    expr_path   <- paste0("geneExp/", store@bin_size, "/expression")

    chunks <- .gef_safe_chunks(name_to_row, store@batch_genes)
    chunk_i <- 0L
    closed  <- FALSE

    # Running (x, y) -> bin_ID lookup. Persists across batches; published
    # to the iterator's metadata accessors when iteration completes.
    xy_to_bin <- data.table::data.table(
        x = integer(0), y = integer(0), bin_ID = integer(0)
    )
    data.table::setkey(xy_to_bin, x, y)
    n_bins <- 0L

    close_fn <- function() { closed <<- TRUE; invisible(NULL) }

    next_batch <- function() {
        # NSE bindings
        x <- y <- bin_ID <- NULL
        repeat {
            if (closed || chunk_i >= length(chunks)) {
                close_fn(); return(NULL)
            }
            chunk_i <<- chunk_i + 1L
            chunk_def <- chunks[[chunk_i]]
            g_lo <- chunk_def[1L]; g_hi <- chunk_def[2L]
            slice_lo <- cum[g_lo] + 1L
            slice_hi <- cum[g_hi + 1L]
            if (slice_hi >= slice_lo) break
        }
        chunk <- data.table::setDT(rhdf5::h5read(
            gef_path, expr_path,
            start = slice_lo,
            count = slice_hi - slice_lo + 1L
        ))
        gene_idx_raw <- rep.int(seq.int(g_lo, g_hi), cnt[g_lo:g_hi])

        # Assign bin_IDs for any (x, y) coords not seen before.
        chunk_xy <- unique(chunk[, .(x, y)])
        new_xy <- chunk_xy[!xy_to_bin, on = c("x", "y")]
        if (nrow(new_xy) > 0L) {
            new_xy[, bin_ID := seq.int(n_bins + 1L, n_bins + .N)]
            n_bins <<- n_bins + nrow(new_xy)
            xy_to_bin <<- rbind(xy_to_bin, new_xy)
            data.table::setkey(xy_to_bin, x, y)
        }
        chunk[, bin_ID := xy_to_bin[.SD, on = c("x", "y"), bin_ID]]

        out <- data.table::data.table(
            row_id = as.integer(chunk$bin_ID),
            col_id = as.integer(name_to_row[gene_idx_raw]),
            value  = as.double(chunk$count)
        )
        out <- out[!is.na(row_id) & !is.na(col_id),
                   .(value = sum(value)), keyby = .(row_id, col_id)]
        out
    }

    # Accessors. cell_ids / n_cells reflect accumulated state, so they
    # change as iteration advances. Callers should query them after the
    # iterator is exhausted.
    list(
        next_batch = next_batch,
        close      = close_fn,
        cell_ids   = function() {
            if (n_bins == 0L) return(character(0L))
            data.table::setorder(xy_to_bin, bin_ID)
            paste0("bin_", seq_len(n_bins))
        },
        feat_ids   = function() store@feat_ids,
        n_cells    = function() n_bins,
        n_genes    = function() store@n_genes,
        # Bonus: expose the (x, y) lookup for downstream spatial use.
        bin_coords = function() {
            data.table::setorder(xy_to_bin, bin_ID)
            data.table::copy(xy_to_bin)
        }
    )
})


# storeRead — csvWideInput ####

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("csvWideInput"), function(store, ...) {
    csv_path  <- store@path
    all_cols  <- store@params$all_cols
    feat_cols <- store@params$feat_cols
    cell_col  <- store@cell_id_col
    batch_rows<- store@batch_rows
    row_filt  <- store@row_filter_fun

    con <- if (.is_real_gz(csv_path)) gzfile(csv_path, "r") else file(csv_path, "r")
    # Skip the header line (we already parsed it at construction).
    readLines(con, n = 1L, warn = FALSE)

    closed         <- FALSE
    cell_ids_acc   <- character(0L)
    n_cells_so_far <- 0L

    close_fn <- function() {
        if (!closed) {
            try(close(con), silent = TRUE)
            closed <<- TRUE
        }
        invisible(NULL)
    }

    next_batch <- function() {
        repeat {
            if (closed) return(NULL)
            raw_lines <- readLines(con, n = batch_rows, warn = FALSE)
            if (length(raw_lines) == 0L) { close_fn(); return(NULL) }
            chunk <- data.table::fread(
                text = raw_lines, header = FALSE, col.names = all_cols
            )
            if (!is.null(row_filt)) {
                keep <- row_filt(chunk)
                chunk <- chunk[keep, ]
            }
            if (length(raw_lines) < batch_rows) {
                # Close eagerly; one more batch may still return below.
                close_fn()
            }
            if (nrow(chunk) > 0L) break
            if (closed) return(NULL)
        }
        chunk_cell_ids <- as.character(chunk[[cell_col]])
        chunk_row_ids  <- seq.int(n_cells_so_far + 1L,
                                  n_cells_so_far + nrow(chunk))
        # Convert wide → triplets via dense-matrix coercion of feat cols.
        # Bounded per-chunk by batch_rows × n_feats × 8 bytes.
        feat_mat <- as.matrix(chunk[, ..feat_cols])
        nz <- which(feat_mat != 0, arr.ind = TRUE)
        cell_ids_acc   <<- c(cell_ids_acc, chunk_cell_ids)
        n_cells_so_far <<- n_cells_so_far + nrow(chunk)
        if (nrow(nz) == 0L) {
            # No nonzeros in this chunk — return an empty data.table with
            # the right schema so the writer doesn't trip on it.
            return(data.table::data.table(
                row_id = integer(0), col_id = integer(0), value = double(0)
            ))
        }
        out <- data.table::data.table(
            row_id = chunk_row_ids[nz[, 1L]],
            col_id = as.integer(nz[, 2L]),
            value  = as.double(feat_mat[nz])
        )
        data.table::setorder(out, row_id, col_id)
        out
    }

    list(
        next_batch = next_batch,
        close      = close_fn,
        cell_ids   = function() cell_ids_acc,
        feat_ids   = function() store@feat_ids,
        n_cells    = function() n_cells_so_far,
        n_genes    = function() store@n_genes
    )
})


# storeWrite — universal exprInput consumer ####

#' @rdname storeWrite
#' @export
setMethod(
    "storeWrite",
    signature("parquetExprStore", "exprInput"),
    function(store, data, ...) {
        # Pre-allocate destination as a directory; one chunk_*.parquet
        # per batch. Single-file is just a 1-batch special case under
        # this scheme — open_dataset() handles both transparently.
        out_path <- store@path
        if (dir.exists(out_path)) {
            # accept existing empty dir; refuse to silently merge into a
            # populated one.
            existing <- list.files(out_path, pattern = "\\.parquet$",
                                   full.names = FALSE)
            if (length(existing) > 0L) {
                stop("[storeWrite] output directory is not empty: ", out_path,
                     "\n  remove it or pass a fresh path before writing.",
                     call. = FALSE)
            }
        } else if (file.exists(out_path)) {
            stop("[storeWrite] output path exists as a file: ", out_path,
                 "\n  pre-allocated store path must be a directory or absent.",
                 call. = FALSE)
        } else {
            dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
        }

        itr <- storeRead(data)
        on.exit(itr$close(), add = TRUE)

        batch_idx <- 0L
        repeat {
            dt <- itr$next_batch()
            if (is.null(dt)) break
            if (nrow(dt) == 0L) next   # skip empty batches (no on-disk chunk)
            batch_idx <- batch_idx + 1L
            arrow::write_parquet(
                dt,
                file.path(out_path,
                          sprintf("chunk_%010d.parquet", batch_idx))
            )
        }

        # Stamp metadata from the iterator's accessors (single source of
        # truth — works for both eagerly-known formats and ones that
        # accumulate cell_ids during iteration, like csvWideInput).
        store@cell_ids <- itr$cell_ids()
        store@feat_ids <- itr$feat_ids()
        store@n_cells  <- as.numeric(itr$n_cells())
        store@n_genes  <- as.numeric(itr$n_genes())
        store
    }
)
