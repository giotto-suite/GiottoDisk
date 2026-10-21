# Zarr array readers over a `zarr_source` (see zarr-source.R).
#
# `.zarr_array()` reads a whole 1D/2D array or a contiguous row range of
# one. The range read decodes ONLY the chunks overlapping the requested
# rows -- this is what keeps repeated block reads over a large array (the
# boundary vertices at 700k+ cells) linear instead of quadratic.
#
# `.zarr_chunk_reader()` is a streaming iterator over a 1D (or 2D Nx1)
# array, used where an array should never be fully materialized (the CFM
# `indices`/`data` walks). Peak memory is one decoded chunk.

# Read a 1D or 2D zarr array.
#
# `range = c(first, last)` (1-based, inclusive) restricts the FIRST axis;
# trailing axes are always read in full. NULL reads everything. Returns a
# vector (1D) or matrix (2D). Missing chunks fill with the metadata
# fill_value. `u4_as_integer` is forwarded to the value decoder (see
# `.zarr_decode_values()`).
.zarr_array <- function(src, prefix, range = NULL, u4_as_integer = FALSE) {
    prefix <- sub("/+$", "", prefix)
    meta <- .zarr_meta(src, prefix)
    dt <- .zarr_meta_dtype(meta, prefix)
    shape <- unlist(meta$shape)
    chunks <- unlist(meta$chunks)
    order <- meta$order %||% "C"

    if (!length(shape) %in% 1:2) {
        stop("[zarr] only 1D/2D arrays supported (got ", length(shape),
            "D at ", prefix, ")", call. = FALSE)
    }

    # doubles: row ranges of large arrays (e.g. nnz) exceed int32
    r0 <- 1
    r1 <- shape[1L]
    if (!is.null(range)) {
        r0 <- as.numeric(range[1L])
        r1 <- as.numeric(range[2L])
        if (is.na(r0) || is.na(r1) || r0 < 1L || r1 > shape[1L] || r1 < r0) {
            stop("[zarr] bad range [", range[1L], ", ", range[2L],
                "] for axis of length ", shape[1L], " at ", prefix,
                call. = FALSE)
        }
    }
    n_rows <- r1 - r0 + 1

    decode <- function(raw_bytes) {
        .zarr_decode_values(
            .zarr_blosc_decompress(raw_bytes, meta, dt),
            dt, u4_as_integer = u4_as_integer
        )
    }
    fill_val <- .zarr_fill_value(meta$fill_value, dt, u4_as_integer)

    # chunk index range overlapping [r0, r1] on the first axis (0-based)
    ci0 <- (r0 - 1) %/% chunks[1L]
    ci1 <- (r1 - 1) %/% chunks[1L]

    if (length(shape) == 1L) {
        out <- vector(typeof(fill_val), n_rows)
        for (ci in ci0:ci1) {
            cs <- ci * chunks[1L] + 1L # global first row of chunk
            ce <- min(cs + chunks[1L] - 1L, shape[1L])
            s <- max(cs, r0)
            e <- min(ce, r1)
            raw_bytes <- .zarr_read_raw(src, paste0(prefix, "/", ci))
            if (is.null(raw_bytes)) {
                out[(s - r0 + 1L):(e - r0 + 1L)] <- fill_val
                next
            }
            vals <- decode(raw_bytes)
            out[(s - r0 + 1L):(e - r0 + 1L)] <- vals[(s - cs + 1L):(e - cs + 1L)]
        }
        return(out)
    }

    # 2D: rows ranged, columns full (possibly chunked on the trailing axis)
    n_col_chunks <- as.integer(ceiling(shape[2L] / chunks[2L]))
    out <- matrix(vector(typeof(fill_val), 1L), nrow = n_rows,
        ncol = shape[2L])
    for (ci in ci0:ci1) {
        cs <- ci * chunks[1L] + 1L
        ce <- min(cs + chunks[1L] - 1L, shape[1L])
        s <- max(cs, r0)
        e <- min(ce, r1)
        chunk_nrow <- chunks[1L] # stored chunk is always full-sized
        for (cj in seq_len(n_col_chunks) - 1L) {
            gs <- cj * chunks[2L] + 1L
            ge <- min(gs + chunks[2L] - 1L, shape[2L])
            raw_bytes <- .zarr_read_raw(src,
                paste0(prefix, "/", ci, ".", cj))
            if (is.null(raw_bytes)) {
                out[(s - r0 + 1L):(e - r0 + 1L), gs:ge] <- fill_val
                next
            }
            vals <- decode(raw_bytes)
            m <- matrix(vals, nrow = chunk_nrow, ncol = chunks[2L],
                byrow = (order == "C"))
            out[(s - r0 + 1L):(e - r0 + 1L), gs:ge] <-
                m[(s - cs + 1L):(e - cs + 1L), seq_len(ge - gs + 1L),
                    drop = FALSE]
        }
    }
    out
}

# Read a row range of a 2D array as a list of column vectors. When the
# trailing axis is a single chunk, each column is copied straight out of
# the decoded chunk (a contiguous slice in F order), with no intermediate
# matrix; otherwise falls back to `.zarr_array()`.
.zarr_columns <- function(src, prefix, range = NULL, u4_as_integer = FALSE) {
    prefix <- sub("/+$", "", prefix)
    meta <- .zarr_meta(src, prefix)
    shape <- unlist(meta$shape)
    chunks <- unlist(meta$chunks)
    if (length(shape) != 2L || chunks[2L] != shape[2L]) {
        m <- .zarr_array(src, prefix, range = range,
            u4_as_integer = u4_as_integer)
        if (!is.matrix(m)) return(list(m))
        return(lapply(seq_len(ncol(m)), function(j) m[, j]))
    }
    dt <- .zarr_meta_dtype(meta, prefix)
    order <- meta$order %||% "C"
    nc <- shape[2L]
    cn <- chunks[1L]
    r0 <- if (is.null(range)) 1 else as.numeric(range[1L])
    r1 <- if (is.null(range)) shape[1L] else as.numeric(range[2L])
    if (r0 < 1 || r1 > shape[1L] || r1 < r0) {
        stop("[zarr] bad range [", r0, ", ", r1, "] for axis of length ",
            shape[1L], " at ", prefix, call. = FALSE)
    }
    fill_val <- .zarr_fill_value(meta$fill_value, dt, u4_as_integer)
    n <- r1 - r0 + 1
    out <- lapply(seq_len(nc), function(j) vector(typeof(fill_val), n))
    for (ci in seq.int((r0 - 1) %/% cn, (r1 - 1) %/% cn)) {
        cs <- ci * cn + 1
        s <- max(cs, r0)
        e <- min(cs + cn - 1, r1)
        dst <- seq.int(s - r0 + 1, e - r0 + 1)
        raw_bytes <- .zarr_read_raw(src, paste0(prefix, "/", ci, ".0"))
        if (is.null(raw_bytes)) {
            for (j in seq_len(nc)) out[[j]][dst] <- fill_val
            next
        }
        vals <- .zarr_decode_values(
            .zarr_blosc_decompress(raw_bytes, meta, dt),
            dt, u4_as_integer = u4_as_integer
        )
        rows <- seq.int(s - cs + 1, e - cs + 1)
        for (j in seq_len(nc)) {
            out[[j]][dst] <- if (order == "F") {
                vals[(j - 1) * cn + rows]
            } else {
                vals[(rows - 1) * nc + j]
            }
        }
    }
    out
}

# Streaming chunk reader for 1D or 2D-Nx1 arrays. Returns
# list(n_chunks, total, chunk_size, read_chunk) where read_chunk(i0b)
# yields the decoded values of chunk i0b (missing chunk -> fill_value,
# truncated final chunk trimmed to the array length).
.zarr_chunk_reader <- function(src, prefix, u4_as_integer = FALSE) {
    prefix <- sub("/+$", "", prefix)
    meta <- .zarr_meta(src, prefix)
    dt <- .zarr_meta_dtype(meta, prefix)
    shape <- unlist(meta$shape)
    chunks <- unlist(meta$chunks)

    if (!(length(shape) == 1L || (length(shape) == 2L && shape[2L] == 1L))) {
        stop("[zarr] chunk reader supports 1D or 2D-Nx1 arrays only; got ",
            paste(shape, collapse = "x"), " at ", prefix, call. = FALSE)
    }
    total <- shape[1L]
    chunk_axis <- chunks[1L]
    n_chunks <- as.integer(ceiling(total / chunk_axis))
    fill_val <- .zarr_fill_value(meta$fill_value, dt, u4_as_integer)
    is_2d <- length(shape) == 2L

    read_chunk <- function(chunk_idx_0b) {
        fpath <- if (is_2d) paste0(prefix, "/", chunk_idx_0b, ".0")
            else paste0(prefix, "/", chunk_idx_0b)
        raw_bytes <- .zarr_read_raw(src, fpath)
        cs <- min(chunk_axis, total - chunk_idx_0b * chunk_axis)
        if (is.null(raw_bytes)) {
            return(rep(fill_val, cs))
        }
        vals <- .zarr_decode_values(
            .zarr_blosc_decompress(raw_bytes, meta, dt),
            dt, u4_as_integer = u4_as_integer
        )
        if (length(vals) > cs) vals <- vals[seq_len(cs)]
        vals
    }

    list(
        n_chunks = n_chunks,
        total = total,
        chunk_size = chunk_axis,
        read_chunk = read_chunk
    )
}
