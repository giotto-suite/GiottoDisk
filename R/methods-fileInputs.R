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
        # `lo:hi`, never `seq.int(lo, hi)` -- adr/0007.
        vals <- as.double(data_ds[slice_lo:slice_hi])
        gidx <- as.integer(indices_ds[slice_lo:slice_hi])
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

    # Columns whose raw gene rows may straddle chunks; their records are
    # held back and flushed as one aggregated batch. See .gef_dup_cols().
    dup_cols <- .gef_dup_cols(name_to_row)
    deferred <- list()
    flushed  <- FALSE

    close_fn <- function() { closed <<- TRUE; invisible(NULL) }

    next_batch <- function() {
        row_id <- col_id <- value <- NULL  # NSE bindings
        repeat {
            if (closed || chunk_i >= length(chunks)) {
                if (!closed && !flushed) {
                    flushed <<- TRUE
                    fl <- .gef_flush_deferred(deferred)
                    deferred <<- list()
                    if (!is.null(fl)) return(fl)
                }
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
        if (length(dup_cols)) {
            hold <- out[col_id %in% dup_cols]
            if (nrow(hold)) {
                deferred[[length(deferred) + 1L]] <<- hold
                out <- out[!col_id %in% dup_cols]
            }
        }
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
    coord_env   <- store@params$coord_env
    expr_path   <- paste0("geneExp/", store@bin_size, "/expression")

    chunks <- .gef_safe_chunks(name_to_row, store@batch_genes)
    chunk_i <- 0L
    closed  <- FALSE
    exhausted <- FALSE

    # Columns whose raw gene rows may straddle chunks; their records are
    # held back and flushed as one aggregated batch. See .gef_dup_cols().
    dup_cols <- .gef_dup_cols(name_to_row)
    deferred <- list()
    flushed  <- FALSE

    # Running (x, y) -> bin_ID lookup. Persists across batches; published
    # to the iterator's metadata accessors when iteration completes.
    xy_to_bin <- data.table::data.table(
        x = integer(0), y = integer(0), bin_ID = integer(0)
    )
    data.table::setkey(xy_to_bin, x, y)
    n_bins <- 0L

    # Hand the finished coordinate map back to the input object, which
    # outlives this iterator. Only on a full pass -- a partial map would
    # silently produce spatial locations for a subset of the bins. See
    # binGefInput()'s `coord_env` note.
    .publish_coords <- function() {
        if (!exhausted || is.null(coord_env)) return(invisible(NULL))
        data.table::setorder(xy_to_bin, bin_ID)
        coord_env$bin_coords <- data.table::copy(xy_to_bin)
        invisible(NULL)
    }

    close_fn <- function() {
        closed <<- TRUE
        .publish_coords()
        invisible(NULL)
    }

    next_batch <- function() {
        # NSE bindings
        x <- y <- bin_ID <- row_id <- col_id <- value <- NULL
        repeat {
            if (closed || chunk_i >= length(chunks)) {
                if (chunk_i >= length(chunks)) exhausted <<- TRUE
                # Flush held-back duplicate-name records before closing. The
                # coordinate map is already complete -- bin_IDs are assigned
                # on the full chunk, ahead of the deferral split.
                if (!closed && !flushed) {
                    flushed <<- TRUE
                    fl <- .gef_flush_deferred(deferred)
                    deferred <<- list()
                    if (!is.null(fl)) return(fl)
                }
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
        if (length(dup_cols)) {
            hold <- out[col_id %in% dup_cols]
            if (nrow(hold)) {
                deferred[[length(deferred) + 1L]] <<- hold
                out <- out[!col_id %in% dup_cols]
            }
        }
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
        # parquetExprStore lays out chunks under a `source_id=<uid>/`
        # hive partition (.idpath helper). Mirrors the parquetStore
        # source_id convention and lets a future unionParquetExprStore
        # aggregate substores by hardlinking their partition dirs.
        out_path <- store@path
        if (file.exists(out_path) && !dir.exists(out_path)) {
            stop("[storeWrite] output path exists as a file: ", out_path,
                 "\n  pre-allocated store path must be a directory or absent.",
                 call. = FALSE)
        }
        partition_dir <- .idpath(out_path, store@uid)
        if (dir.exists(partition_dir)) {
            existing <- list.files(partition_dir, pattern = "\\.parquet$",
                                   full.names = FALSE)
            if (length(existing) > 0L) {
                stop("[storeWrite] partition is not empty: ", partition_dir,
                     "\n  remove it or pass a fresh store before writing.",
                     call. = FALSE)
            }
        } else {
            dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)
        }

        # tenxH5Input has an embarrassingly-parallel batch structure
        # (disjoint cell ranges, each worker opens its own h5 handle and
        # writes its own shard). Route through .storewrite_h5_parallel
        # when a parallel future plan is set. All other inputs, and
        # tenxH5Input under a sequential plan, use the serial iterator —
        # avoids lapply_flex ceremony (and the sequential-plan warning)
        # when there's no throughput to gain.
        if (inherits(data, "tenxH5Input") && .have_parallel_plan()) {
            .storewrite_h5_parallel(data, partition_dir)
            store@cell_ids <- data@cell_ids
            store@feat_ids <- data@feat_ids
            store@n_cells  <- as.numeric(data@n_cells)
            store@n_genes  <- as.numeric(data@n_genes)
            .pestore_finalize_stats(store)
            return(invisible(store))
        }

        itr <- storeRead(data)
        on.exit(itr$close(), add = TRUE)

        batch_idx <- 0L
        repeat {
            dt <- itr$next_batch()
            if (is.null(dt)) break
            if (nrow(dt) == 0L) next   # skip empty batches (no on-disk chunk)
            batch_idx <- batch_idx + 1L
            .write_parquet_file(
                dt,
                file.path(partition_dir,
                          sprintf("part-%d.parquet", batch_idx - 1L))
            )
        }

        # Stamp metadata from the iterator's accessors (single source of
        # truth — works for both eagerly-known formats and ones that
        # accumulate cell_ids during iteration, like csvWideInput).
        store@cell_ids <- itr$cell_ids()
        store@feat_ids <- itr$feat_ids()
        store@n_cells  <- as.numeric(itr$n_cells())
        store@n_genes  <- as.numeric(itr$n_genes())
        .pestore_finalize_stats(store)
    }
)


# ---- tenxH5Input parallel storeWrite ---------------------------------------

# Number of workers this package should use for band-parallel steps.
# Priority order:
#   1. Explicit `options("giottodisk.par_workers")` (int).
#   2. `future::nbrOfWorkers()` if future is installed AND a non-uniprocess
#      plan is set (backwards-compatible with users who set only a plan).
#   3. Otherwise 1 (serial).
# On Windows, `parallel::mclapply(mc.cores = n)` silently degrades to
# sequential, so callers must guard their own fork-path branches.
.par_workers <- function() {
    n <- getOption("giottodisk.par_workers", NULL)
    if (!is.null(n)) return(max(1L, as.integer(n)))
    if (requireNamespace("future", quietly = TRUE) &&
        !inherits(future::plan(), "uniprocess")) {
        return(max(1L, as.integer(future::nbrOfWorkers())))
    }
    1L
}

# Legacy name; still used by `.storewrite_h5_parallel` on the lapply_flex
# path. Now defined in terms of `.par_workers()` so both stay in sync.
.have_parallel_plan <- function() .par_workers() > 1L

# Read one cell-range slice from a 10x .h5 (CSC-by-cell layout) and return
# a sorted (row_id, col_id, value) data.table. `indptr` must be supplied
# (a numeric vector of length n_cells+1). Reading it in the caller keeps
# workers from paying 8x for the same full-array read.
.tenxh5_read_range <- function(h5_path, root, c_lo, c_hi, indptr) {
    slice_lo <- as.numeric(indptr[c_lo]) + 1
    slice_hi <- as.numeric(indptr[c_hi + 1L])
    if (slice_hi < slice_lo) return(NULL)   # empty batch
    h5 <- hdf5r::H5File$new(h5_path, mode = "r")
    on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
    # Index as `lo:hi`, NOT `seq.int(lo, hi)`. `:` is an ALTREP compact
    # sequence that hdf5r serves as one hyperslab; a materialized vector falls
    # back to point selection, at 4x the time and memory for identical output.
    # Load-bearing, and invisible at the call site -- see adr/0007.
    vals <- as.double(h5[[paste0(root, "/data")]][slice_lo:slice_hi])
    gidx <- as.integer(h5[[paste0(root, "/indices")]][slice_lo:slice_hi])
    nnz_per_cell   <- diff(indptr[c_lo:(c_hi + 1L)])
    cell_idx_local <- rep.int(seq.int(c_lo, c_hi), nnz_per_cell)
    dt <- data.table::data.table(
        row_id = as.integer(cell_idx_local),
        col_id = gidx + 1L,        # h5 indices are 0-based
        value  = vals
    )
    data.table::setorder(dt, row_id, col_id)
    dt
}

# Parallel h5 → parquet writer.  Fork-based via `parallel::mclapply` on
# Unix — workers COW-inherit the parent's loaded namespace and the
# pre-read indptr;
# each opens its own hdf5r handle (required — hdf5r C handles aren't
# fork-safe across parent+child sharing).  On Windows, mclapply degrades
# to sequential (no fork available).
.storewrite_h5_parallel <- function(data, partition_dir) {
    n_cells <- as.integer(data@n_cells)
    bc      <- as.integer(data@batch_cells)
    starts  <- seq.int(1L, n_cells, by = bc)
    h5_path <- data@path
    root    <- data@params$root
    # Read indptr once in the parent (small; n_cells+1 doubles) and pass it
    # to each batch as an argument, so no worker repeats the full-array read.
    h5     <- hdf5r::H5File$new(h5_path, mode = "r")
    indptr <- as.numeric(h5[[paste0(root, "/indptr")]][])
    h5$close_all()

    n_workers <- .par_workers()
    args <- list(h5_path = h5_path, root = root, starts = starts, bc = bc,
                 n_cells = n_cells, indptr = indptr,
                 partition_dir = partition_dir)

    # Fork (`mclapply`) where available, socket workers only as the fallback.
    # Measured on Atera (170k cells, 8 workers): fork 9.0 s / +16.1 GB peak
    # versus lapply_flex -> mirai 17.8 s / +8.7 GB with daemons already warm.
    # Fork wins on wall-clock by ~2x because each task ships `indptr` and its
    # arguments over a socket and pays per-daemon namespace loading, while a
    # fork inherits everything for free.  The socket path does roughly halve
    # peak RSS, so it is the better choice when memory, not time, is the
    # binding constraint -- and it is the only option on Windows, where
    # `mclapply(mc.cores > 1)` silently degrades to sequential.
    #
    # `.storewrite_h5_batch` is deliberately a TOP-LEVEL internal rather than
    # a local closure so it works under both: its environment is the
    # GiottoDisk namespace, which serializes by name reference (~3 KB)
    # instead of by value.  A closure defined here would carry this whole
    # call frame -- including `indptr` -- into every task.  That also makes
    # `future.packages` unnecessary: deserializing a namespace-parented
    # function loads the namespace in the worker, so sibling internals
    # resolve lexically and every cross-package call here is `::`-qualified.
    if (n_workers > 1L && .Platform$OS.type == "unix") {
        res <- do.call(parallel::mclapply, c(
            list(X = seq_along(starts), FUN = .storewrite_h5_batch,
                 mc.cores = n_workers, mc.preschedule = TRUE),
            args))
        errs <- vapply(res, inherits, logical(1L), "try-error")
        if (any(errs)) {
            msg <- attr(res[[which(errs)[1L]]], "condition")$message
            stop("[.storewrite_h5_parallel] worker failed: ", msg,
                 call. = FALSE)
        }
    } else if (n_workers > 1L) {
        do.call(GiottoUtils::lapply_flex, c(
            list(X = seq_along(starts), FUN = .storewrite_h5_batch,
                 cores = n_workers, future.seed = NULL),
            args))
    } else {
        do.call(lapply, c(
            list(X = seq_along(starts), FUN = .storewrite_h5_batch),
            args))
    }
    invisible(NULL)
}

# One h5 cell-batch -> one parquet shard. Top-level (not a closure) so it
# serializes as a namespace reference; every input arrives as an argument.

#' @keywords internal
#' @noRd
.storewrite_h5_batch <- function(b, h5_path, root, starts, bc, n_cells,
                                  indptr, partition_dir) {
    c_lo <- starts[b]
    c_hi <- min(c_lo + bc - 1L, n_cells)
    dt   <- .tenxh5_read_range(h5_path, root, c_lo, c_hi, indptr)
    if (is.null(dt) || nrow(dt) == 0L) return(invisible(NULL))
    .write_parquet_file(
        dt,
        file.path(partition_dir, sprintf("part-%d.parquet", b - 1L)))
    invisible(NULL)
}


# storeRead -- tenxZarrInput ####

# The zarr cell_feature_matrix is CSC by FEATURE (indptr over features,
# `indices` = cell ids), but the iterator contract is ascending-cell
# batches. Two strategies, picked at read time:
#
# full      two-pass cell-major placement ("B2"): pass 1 streams `indices`
#           chunks into a per-cell nnz histogram; pass 2 streams
#           `indices` + `data` in lockstep, placing each triplet at its
#           cell's cursor in preallocated int buffers. Features are
#           walked in ascending order, so within-cell col_id is ascending
#           by construction -- no sort. Peak RAM ~ 12 B x nnz.
# cellblock bounded windows over the cell axis: each window rescans the
#           `indices`/`data` chunks, masks entries to the window, and
#           sorts the bounded result. RAM is bounded by the window
#           regardless of nnz (rescans over spill, cf. adr/0011); forced
#           when nnz exceeds int32 (full-mode buffers are impossible).
#
# Both modes emit strictly ascending cell ranges; an out-of-order batch
# is a bug and stops rather than corrupting the sorted-parquet invariant.

#' @rdname storeRead
#' @export
setMethod("storeRead", signature("tenxZarrInput"), function(store, ...) {
    src <- .zarr_open(store@path)
    indptr <- store@params$indptr # double, 0-based offsets, length n_feat+1
    keep <- store@params$keep
    col_remap <- store@params$col_remap
    n_cells <- as.integer(store@n_cells)
    n_feat <- length(keep)
    nnz <- store@nnz
    batch_rows <- 5000000L

    ind_r <- .zarr_chunk_reader(src, "cell_features/indices",
        u4_as_integer = TRUE)
    dat_r <- .zarr_chunk_reader(src, "cell_features/data",
        u4_as_integer = TRUE)
    if (ind_r$total != dat_r$total || ind_r$n_chunks != dat_r$n_chunks) {
        .zarr_close(src)
        stop("[storeRead] cell_features indices/data disagree in shape ",
            "or chunking", call. = FALSE)
    }

    # call body(j, s, e) for every feature j whose nnz range overlaps the
    # chunk at 0-based offset chunk_off; s:e are 1-based local indices
    for_features_in_chunk <- function(chunk_off, chunk_len, body) {
        chunk_end <- chunk_off + chunk_len
        first_j <- max(1L, findInterval(chunk_off, indptr))
        last_j <- min(n_feat, findInterval(chunk_end - 1, indptr))
        if (last_j < first_j) return(invisible(NULL))
        for (j in first_j:last_j) {
            s <- max(0, indptr[j] - chunk_off) + 1
            e <- min(chunk_len, indptr[j + 1L] - chunk_off)
            if (e < s) next
            body(j, as.integer(s), as.integer(e))
        }
        invisible(NULL)
    }

    # resolve mode
    mode <- store@mode
    block <- as.integer(store@cells_per_block)
    if (block <= 0L) {
        density <- if (n_cells > 0L && store@n_genes > 0L) {
            nnz / (as.double(n_cells) * as.double(store@n_genes))
        } else {
            1
        }
        block <- .recommend_chunk_size(
            n_cells = n_cells, n_genes = store@n_genes,
            density = density, k = 0L, bytes_per_nz = 16
        )
    }
    if (mode == "auto") {
        mode <- if (nnz > .Machine$integer.max || block < n_cells) {
            "cellblock"
        } else {
            "full"
        }
    }
    if (mode == "full" && nnz > .Machine$integer.max) {
        .zarr_close(src)
        stop("[storeRead] nnz (", format(nnz, big.mark = ","),
            ") exceeds int32; use mode = \"cellblock\"", call. = FALSE)
    }

    closed <- FALSE
    close_fn <- function() {
        if (!closed) {
            .zarr_close(src)
            closed <<- TRUE
        }
        invisible(NULL)
    }
    last_row_emitted <- 0L
    check_ascending <- function(dt) {
        if (nrow(dt) && dt$row_id[1L] <= last_row_emitted) {
            close_fn()
            stop("[storeRead] internal: batch starts at cell ",
                dt$row_id[1L], " but cells through ", last_row_emitted,
                " were already emitted -- ascending-cell invariant broken",
                call. = FALSE)
        }
        if (nrow(dt)) last_row_emitted <<- dt$row_id[nrow(dt)]
        dt
    }

    if (mode == "full") {
        built <- FALSE
        row_ids <- col_ids <- values <- NULL
        off <- NULL
        lo_cell <- 1L

        build <- function() {
            # pass 1: per-cell nnz over kept features
            per_cell_nnz <- integer(n_cells)
            for (k in seq_len(ind_r$n_chunks)) {
                ic <- ind_r$read_chunk(k - 1L)
                chunk_off <- (k - 1) * ind_r$chunk_size
                for_features_in_chunk(chunk_off, length(ic),
                    function(j, s, e) {
                        if (!keep[j]) return(invisible(NULL))
                        cells_1b <- ic[s:e] + 1L
                        per_cell_nnz[cells_1b] <<-
                            per_cell_nnz[cells_1b] + 1L
                    })
            }
            off <<- c(0L, cumsum(per_cell_nnz))
            total <- off[n_cells + 1L]
            row_ids <<- integer(total)
            col_ids <<- integer(total)
            values <<- integer(total)
            cursor <- integer(n_cells)
            # pass 2: lockstep placement
            for (k in seq_len(ind_r$n_chunks)) {
                ic <- ind_r$read_chunk(k - 1L)
                dc <- dat_r$read_chunk(k - 1L)
                chunk_off <- (k - 1) * ind_r$chunk_size
                for_features_in_chunk(chunk_off, length(ic),
                    function(j, s, e) {
                        if (!keep[j]) return(invisible(NULL))
                        cells_1b <- ic[s:e] + 1L
                        slots <- off[cells_1b] + cursor[cells_1b] + 1L
                        row_ids[slots] <<- cells_1b
                        col_ids[slots] <<- col_remap[j]
                        values[slots] <<- dc[s:e]
                        cursor[cells_1b] <<- cursor[cells_1b] + 1L
                    })
            }
            if (!identical(cursor, per_cell_nnz)) {
                close_fn()
                stop("[storeRead] internal: cell-major placement drifted ",
                    "(cursor != per-cell nnz)", call. = FALSE)
            }
            built <<- TRUE
            invisible(NULL)
        }

        next_batch <- function() {
            if (closed) return(NULL)
            if (!built) build()
            repeat {
                if (lo_cell > n_cells) {
                    close_fn()
                    return(NULL)
                }
                # widest hi with off[hi + 1] - off[lo] <= batch_rows,
                # but at least one cell
                hi_cell <- findInterval(off[lo_cell] + batch_rows, off) - 1L
                hi_cell <- max(lo_cell, min(hi_cell, n_cells))
                s <- off[lo_cell] + 1L
                e <- off[hi_cell + 1L]
                lo_cell <<- hi_cell + 1L
                if (e >= s) {
                    return(check_ascending(data.table::data.table(
                        row_id = row_ids[s:e],
                        col_id = col_ids[s:e],
                        value = as.double(values[s:e])
                    )))
                }
            }
        }
    } else { # cellblock
        lo_cell <- 1L
        next_batch <- function() {
            if (closed) return(NULL)
            repeat {
                if (lo_cell > n_cells) {
                    close_fn()
                    return(NULL)
                }
                hi_cell <- min(lo_cell + block - 1L, n_cells)
                parts <- list()
                for (k in seq_len(ind_r$n_chunks)) {
                    ic <- ind_r$read_chunk(k - 1L)
                    chunk_len <- length(ic)
                    # dense col id per entry (0 = dropped feature)
                    cols <- integer(chunk_len)
                    chunk_off <- (k - 1) * ind_r$chunk_size
                    for_features_in_chunk(chunk_off, chunk_len,
                        function(j, s, e) {
                            cols[s:e] <<- col_remap[j]
                        })
                    cells_1b <- ic + 1L
                    m <- cols > 0L & cells_1b >= lo_cell &
                        cells_1b <= hi_cell
                    if (!any(m)) next
                    dc <- dat_r$read_chunk(k - 1L)
                    parts[[length(parts) + 1L]] <- data.table::data.table(
                        row_id = cells_1b[m],
                        col_id = cols[m],
                        value = as.double(dc[m])
                    )
                }
                lo_cell <<- hi_cell + 1L
                if (length(parts)) {
                    out <- data.table::rbindlist(parts)
                    data.table::setorder(out, row_id, col_id)
                    return(check_ascending(out))
                }
            }
        }
    }

    list(
        next_batch = next_batch,
        close = close_fn,
        cell_ids = function() store@cell_ids,
        feat_ids = function() store@feat_ids,
        n_cells = function() store@n_cells,
        n_genes = function() store@n_genes
    )
})
