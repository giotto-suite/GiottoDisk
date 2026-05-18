#' @include class-fileInputs.R class-parquetExprStore.R
NULL

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

    list(next_batch = next_batch, close = close_fn)
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
            batch_idx <- batch_idx + 1L
            arrow::write_parquet(
                dt,
                file.path(out_path,
                          sprintf("chunk_%010d.parquet", batch_idx))
            )
        }

        # Stamp metadata from the input. cell_ids / feat_ids etc. are
        # authoritative on the input (read eagerly at construction or
        # accumulated during iteration) — adopt them into the store.
        store@cell_ids <- data@cell_ids
        store@feat_ids <- data@feat_ids
        store@n_cells  <- as.numeric(data@n_cells)
        store@n_genes  <- as.numeric(data@n_genes)
        store
    }
)
