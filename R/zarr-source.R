# Zarr v2 source layer: read-only access to Xenium/Atera `.zarr.zip`
# archives (in place, no extraction) and unzipped `.zarr` directory trees.
#
# A `zarr_source` wraps either a filesystem directory or a `.zarr.zip`
# archive. The zip variant caches the central directory once via
# `.zip_central_directory()` and serves per-entry reads with a seek + 30-byte
# local-header parse + raw read on a single long-lived connection. This
# relies on Xenium archives storing entries UNCOMPRESSED in the zip (each
# chunk is already blosc-compressed) -- `.zarr_open()` refuses archives
# that violate this.
#
# Blosc/zlib chunk decompression goes through `.zarr_blosc_decompress()`,
# the single call site of the internal `Rarr:::.decompress_chunk`.
#
# Rarr and zip live in Suggests; `.zarr_open()` is the gate.
#
# All functions here take zarr paths as forward-slash strings without a
# leading slash, e.g. "polygon_sets/0/vertices" or "cell_features/indptr".

# open / close ####

# Verify the internal Rarr API this layer depends on. Tested against
# Rarr 1.10.1; the signature is (compressed_chunk, metadata).
.zarr_check_rarr <- function() {
    GiottoUtils::package_check("Rarr", repository = "Bioc")
    if (!is.function(get0(".decompress_chunk", envir = asNamespace("Rarr")))) {
        stop("[zarr] installed Rarr (", utils::packageVersion("Rarr"),
            ") does not provide the internal `.decompress_chunk` this ",
            "reader was built against (tested with Rarr 1.10.1). ",
            "Please report this to the GiottoDisk maintainers.",
            call. = FALSE)
    }
    invisible(NULL)
}

# Central directory entries: filename, local-header offset and compressed size.
# Replaces zip::zip_list(), which cannot parse a zip64 central directory.
# Xenium writes one for transcripts.zarr.zip whatever the archive size.
.zip_central_directory <- function(path) {
    size <- file.info(path)$size
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)

    rd <- function(pos, n) {
        seek(con, pos)
        readBin(con, "raw", n)
    }
    u16 <- function(r, i) sum(as.numeric(r[i + 1:2]) * 256^(0:1))
    u32 <- function(r, i) sum(as.numeric(r[i + 1:4]) * 256^(0:3))
    u64 <- function(r, i) sum(as.numeric(r[i + 1:8]) * 256^(0:7))
    zip64_na <- 4294967295

    tail_n <- min(size, 65557)
    tl <- rd(size - tail_n, tail_n)
    n <- length(tl)
    hits <- which(
        tl[1:(n - 3)] == as.raw(0x50) & tl[2:(n - 2)] == as.raw(0x4b) &
            tl[3:(n - 1)] == as.raw(0x05) & tl[4:n] == as.raw(0x06)
    )
    if (!length(hits)) {
        stop("[zip] no end-of-central-directory record in ", path,
            call. = FALSE)
    }
    e <- hits[length(hits)] - 1L
    n_entries <- u16(tl, e + 10)
    cd_size <- u32(tl, e + 12)
    cd_offset <- u32(tl, e + 16)

    if (e >= 20) {
        loc <- tl[(e - 20) + 1:20]
        if (identical(loc[1:4], as.raw(c(0x50, 0x4b, 0x06, 0x07)))) {
            z <- rd(u64(loc, 8), 56)
            if (identical(z[1:4], as.raw(c(0x50, 0x4b, 0x06, 0x06)))) {
                n_entries <- u64(z, 32)
                cd_size <- u64(z, 40)
                cd_offset <- u64(z, 48)
            }
        }
    }

    cd <- rd(cd_offset, cd_size)
    fn <- character(n_entries)
    off <- numeric(n_entries)
    csz <- numeric(n_entries)
    p <- 0
    for (i in seq_len(n_entries)) {
        if (!identical(cd[p + 1:4], as.raw(c(0x50, 0x4b, 0x01, 0x02)))) {
            stop("[zip] malformed central directory entry ", i, " in ", path,
                call. = FALSE)
        }
        c_size <- u32(cd, p + 20)
        u_size <- u32(cd, p + 24)
        n_len <- u16(cd, p + 28)
        e_len <- u16(cd, p + 30)
        k_len <- u16(cd, p + 32)
        lho <- u32(cd, p + 42)
        fn[i] <- rawToChar(cd[p + 46 + seq_len(n_len)])
        if (e_len > 0 &&
            (c_size == zip64_na || u_size == zip64_na || lho == zip64_na)) {
            ex <- cd[p + 46 + n_len + seq_len(e_len)]
            q <- 0
            while (q + 4 <= e_len) {
                if (u16(ex, q) == 1) {
                    r <- q + 4
                    if (u_size == zip64_na) {
                        u_size <- u64(ex, r)
                        r <- r + 8
                    }
                    if (c_size == zip64_na) {
                        c_size <- u64(ex, r)
                        r <- r + 8
                    }
                    if (lho == zip64_na) lho <- u64(ex, r)
                    break
                }
                q <- q + 4 + u16(ex, q + 2)
            }
        }
        off[i] <- lho
        csz[i] <- c_size
        p <- p + 46 + n_len + e_len + k_len
    }
    list(filename = fn, offset = off, compressed_size = csz)
}

.zarr_open <- function(path) {
    .zarr_check_rarr()
    meta_cache <- new.env(parent = emptyenv())
    if (grepl("\\.zip$", path, ignore.case = TRUE)) {
        if (!file.exists(path)) {
            stop("[zarr] archive does not exist: ", path, call. = FALSE)
        }
        entries <- .zip_central_directory(path)
        con <- file(path, open = "rb")
        structure(
            list(
                kind = "zip", path = path, con = con,
                filenames = as.character(entries$filename),
                offset = as.numeric(entries$offset),
                compressed_size = as.numeric(entries$compressed_size),
                meta_cache = meta_cache
            ),
            class = "zarr_source"
        )
    } else {
        if (!dir.exists(path)) {
            stop("[zarr] source does not exist: ", path, call. = FALSE)
        }
        structure(
            list(kind = "dir", path = path, meta_cache = meta_cache),
            class = "zarr_source"
        )
    }
}

.zarr_close <- function(src) {
    if (!is.null(src) && identical(src$kind, "zip") && !is.null(src$con)) {
        # isOpen() errors (rather than returning FALSE) once a connection
        # has been destroyed, so closing must stay idempotent through it
        is_open <- tryCatch(isOpen(src$con), error = function(e) FALSE)
        if (is_open) try(close(src$con), silent = TRUE)
    }
    invisible(NULL)
}

# entry-level access ####

.zarr_exists <- function(src, entry) {
    if (src$kind == "dir") {
        file.exists(file.path(src$path, entry)) ||
            dir.exists(file.path(src$path, entry))
    } else {
        entry %in% src$filenames ||
            any(startsWith(src$filenames, paste0(entry, "/")))
    }
}

# Direct children under `prefix` ("" = root). `dirs_only` keeps only
# entries that are themselves groups/arrays (zip: have a slash below).
.zarr_list <- function(src, prefix, dirs_only = FALSE) {
    prefix <- sub("/+$", "", prefix)
    if (src$kind == "dir") {
        p <- if (nzchar(prefix)) file.path(src$path, prefix) else src$path
        if (!dir.exists(p)) return(character(0L))
        if (dirs_only) {
            basename(list.dirs(p, recursive = FALSE))
        } else {
            list.files(p, recursive = FALSE)
        }
    } else {
        if (nzchar(prefix)) {
            pfx <- paste0(prefix, "/")
            hit <- startsWith(src$filenames, pfx)
            inner <- substring(src$filenames[hit], nchar(pfx) + 1L)
        } else {
            inner <- src$filenames
        }
        is_dir <- grepl("/", inner, fixed = TRUE)
        children <- ifelse(is_dir, sub("/.*$", "", inner), inner)
        if (dirs_only) {
            unique(children[is_dir])
        } else {
            unique(children[nzchar(children)])
        }
    }
}

.zarr_read_raw <- function(src, entry) {
    if (src$kind == "dir") {
        fpath <- file.path(src$path, entry)
        if (!file.exists(fpath)) return(NULL)
        readBin(fpath, what = "raw", n = file.info(fpath)$size)
    } else {
        hit <- which(src$filenames == entry)
        if (length(hit) == 0L) return(NULL)
        local_off <- src$offset[hit[1L]]
        csize <- src$compressed_size[hit[1L]]
        con <- src$con
        # Local file header layout (PKZIP APPNOTE 4.3.7):
        #   sig(4) ver(2) flag(2) method(2) time(2) date(2) crc(4)
        #   csize(4) usize(4) fname_len(2) extra_len(2) [name] [extra] [data]
        # One read from +8 covers method (word 1) and the name/extra
        # lengths (words 10, 11). `compressed_size` bytes go straight to
        # the blosc decoder, which is only valid for STORED (method 0)
        # entries -- Xenium/Atera archives always store; refuse anything
        # else rather than decode garbage.
        seek(con, where = local_off + 8, origin = "start", rw = "read")
        hdr <- readBin(con, what = "integer", size = 2L, n = 11L,
            signed = FALSE, endian = "little")
        if (hdr[1L] != 0L) {
            stop("[zarr] ", src$path, " entry '", entry, "' is ",
                "compressed inside the zip (method ", hdr[1L], "); only ",
                "stored .zarr.zip archives are supported. Unzip the ",
                "archive and pass the directory instead.", call. = FALSE)
        }
        seek(con, where = local_off + 30 + hdr[10L] + hdr[11L],
            origin = "start", rw = "read")
        readBin(con, what = "raw", n = csize)
    }
}

.zarr_read_text <- function(src, entry) {
    raw_vec <- .zarr_read_raw(src, entry)
    if (is.null(raw_vec)) return(NULL)
    rawToChar(raw_vec)
}

# metadata (memoized per source handle) ####

.zarr_meta <- function(src, prefix) {
    prefix <- sub("/+$", "", prefix)
    key <- paste0("zarray::", prefix)
    hit <- get0(key, envir = src$meta_cache, inherits = FALSE)
    if (!is.null(hit)) return(hit)
    json <- .zarr_read_text(src, paste0(prefix, "/.zarray"))
    if (is.null(json)) {
        stop("[zarr] .zarray not found at ", prefix, call. = FALSE)
    }
    meta <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    assign(key, meta, envir = src$meta_cache)
    meta
}

.zarr_attrs <- function(src, prefix) {
    prefix <- sub("/+$", "", prefix)
    key <- paste0("zattrs::", prefix)
    hit <- get0(key, envir = src$meta_cache, inherits = FALSE)
    if (!is.null(hit)) return(hit)
    entry <- if (nzchar(prefix)) paste0(prefix, "/.zattrs") else ".zattrs"
    json <- .zarr_read_text(src, entry)
    if (is.null(json)) return(NULL)
    attrs <- jsonlite::fromJSON(json, simplifyVector = TRUE)
    assign(key, attrs, envir = src$meta_cache)
    attrs
}

# dtype / decompression ####

# Parse a zarr v2 dtype string ("<u4", "<f4", "|u1", "u1", ...) into readBin
# parameters. `kind` is one of u/i/f/b; `size` in bytes. The byte-order prefix
# is optional for single-byte types.
.parse_zarr_dtype <- function(dtype) {
    bo <- substr(dtype, 1L, 1L)
    has_bo <- bo %in% c("<", ">", "|", "=")
    rest <- if (has_bo) substring(dtype, 2L) else dtype
    list(
        kind = substr(rest, 1L, 1L),
        size = as.integer(substring(rest, 2L)),
        # "|" (n/a) and "=" (native) behave as little on Xenium hardware
        endian = if (has_bo && bo == ">") "big" else "little"
    )
}

# Decompress one chunk's raw bytes via Rarr's internal blosc/zlib decoder.
# THE single call site of `Rarr:::.decompress_chunk` (internal API;
# existence verified in `.zarr_open()`). Returns raw bytes ready for
# readBin() reinterpretation.
.zarr_blosc_decompress <- function(raw_bytes, meta, dt) {
    base_type <- switch(dt$kind,
        "f" = "float", "i" = "int", "u" = "uint", "b" = "uint",
        "unknown")
    rarr_meta <- list(
        datatype = list(
            nbytes = dt$size, base_type = base_type,
            is_signed = (dt$kind == "i")
        ),
        compressor = meta$compressor,
        chunks = meta$chunks,
        order = meta$order %||% "C",
        fill_value = meta$fill_value %||% 0
    )
    Rarr:::.decompress_chunk(raw_bytes, rarr_meta)
}

# Reinterpret decompressed chunk bytes as an R vector per the dtype.
#   f4/f8      -> double
#   u4         -> double with the uint32 sign shift (values > 2^31 exact),
#                 or R integer without the shift when `u4_as_integer`
#                 (caller asserts values < 2^31; halves per-chunk memory)
#   u8         -> double via lo/hi uint32 split (exact to 2^53)
#   u1/u2/i1/i2/i4 -> integer
#   b1         -> logical
.zarr_decode_values <- function(decomp, dt, u4_as_integer = FALSE) {
    n_vals <- length(decomp) / dt$size
    if (dt$kind == "f") {
        readBin(decomp, what = "double", size = dt$size, n = n_vals,
            endian = dt$endian)
    } else if (dt$kind == "u" && dt$size == 4L) {
        ints <- readBin(decomp, what = "integer", size = 4L, n = n_vals,
            signed = TRUE, endian = dt$endian)
        if (isTRUE(u4_as_integer)) return(ints)
        out <- as.double(ints)
        neg <- !is.na(out) & out < 0
        out[neg] <- out[neg] + 4294967296
        out
    } else if (dt$kind == "u" && dt$size == 8L) {
        # uint64 as double: exact up to 2^53, ample for nnz offsets.
        words <- readBin(decomp, what = "integer", size = 4L, n = n_vals * 2,
            signed = TRUE, endian = dt$endian)
        u32 <- as.double(words)
        neg <- !is.na(u32) & u32 < 0
        u32[neg] <- u32[neg] + 4294967296
        if (dt$endian == "little") {
            lo <- u32[seq.int(1L, length(u32), by = 2L)]
            hi <- u32[seq.int(2L, length(u32), by = 2L)]
        } else {
            hi <- u32[seq.int(1L, length(u32), by = 2L)]
            lo <- u32[seq.int(2L, length(u32), by = 2L)]
        }
        lo + hi * 4294967296
    } else if (dt$kind == "u" || dt$kind == "i") {
        if (dt$size == 8L) {
            stop("[zarr] int64 arrays are not supported", call. = FALSE)
        }
        readBin(decomp, what = "integer", size = dt$size, n = n_vals,
            signed = dt$kind == "i", endian = dt$endian)
    } else if (dt$kind == "b") {
        readBin(decomp, what = "logical", size = 1L, n = n_vals)
    } else {
        stop("[zarr] unsupported dtype kind: ", dt$kind, call. = FALSE)
    }
}
