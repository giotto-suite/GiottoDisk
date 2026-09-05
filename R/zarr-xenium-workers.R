# Xenium/Atera zarr -> parquet conversion workers.
#
# These produce UNFILTERED, UNFLIPPED parquet in the 10x-shipped schemas
# (cells.parquet, {cell,nucleus}_boundaries.parquet, transcripts.parquet).
# Downstream qv filtering and the y flip stay lazy in the disk readers'
# read_funs (convenience-xenium.R), so converted output is independent of
# reader parameters and directly comparable to 10x-shipped files.
#
# Boundary output is written sequentially in cell order as a SINGLE
# parquet file — the polygon ingest locks vertex order via a row_index
# intermediate, and a multi-file dataset does not guarantee scan order.
#
# Known divergences from 10x-shipped parquet (not derivable from zarr):
# transcripts `cell_id` ("UNASSIGNED"; reassigned downstream by polygon
# overlap), `overlaps_nucleus` (0), `nucleus_distance` (NA), and
# `fov_name` (synthetic "FOV%03d" — the alphanumeric codes 10x ships come
# from per-instrument config absent from the zarr).

# shared helpers ####

# Encode uint32 prefixes to the 8-char a-p alphabet (10x cell_id
# encoding), suffixed with the dataset index. Vectorized; full uint32
# range via hi/lo split.
.encode_xenium_id <- function(prefix, suffix = 1L) {
    prefix <- as.double(prefix)
    hi <- as.integer(prefix %/% 65536)
    lo <- as.integer(prefix %% 65536)
    hex <- paste0(sprintf("%04x", hi), sprintf("%04x", lo))
    paste0(
        chartr("0123456789abcdef", "abcdefghijklmnop", hex),
        "-", as.integer(suffix)
    )
}

# gene_identity (0-based) -> feature name lookup from gene_panel.json.
# Returns character(0) when the panel file is absent (caller falls back
# to the CFM feature catalog, then errors).
.load_gene_panel_lookup <- function(xenium_dir) {
    gp_path <- file.path(xenium_dir, "gene_panel.json")
    if (!file.exists(gp_path)) return(character(0L))
    gp <- jsonlite::fromJSON(gp_path, simplifyVector = FALSE)
    vapply(gp$payload$targets, function(t) {
        nm <- t$type$data$name
        id <- t$type$data$id
        if (!is.null(nm) && nzchar(nm)) nm
        else if (!is.null(id) && nzchar(id)) id
        else paste0(t$type$descriptor, "_unknown")
    }, character(1L))
}

# cell_feature_matrix .zattrs catalog.
.load_cfm_feature_catalog <- function(src) {
    attrs <- .zarr_attrs(src, "cell_features")
    if (is.null(attrs)) {
        stop("[zarr] cell_features/.zattrs not found in source",
            call. = FALSE)
    }
    list(
        feature_ids = as.character(attrs$feature_ids),
        feature_keys = as.character(attrs$feature_keys),
        feature_types = as.character(attrs$feature_types),
        n_cells = as.integer(attrs$number_cells)
    )
}

# Compression codec with graceful fallback for arrow builds without zstd.
# zstd is the default: measured ~27% smaller than snappy at equal write
# speed on the reference Xenium dataset, with negligible scan cost.
.zarr_codec <- function(compression) {
    if (arrow::codec_is_available(compression)) return(compression)
    "snappy"
}

# cells ####

# cells.zarr -> cells.parquet (10x cellmeta schema + segmentation_method).
.zarr_cells_to_parquet <- function(src, out_path, compression = "zstd",
    verbose = NULL) {
    t0 <- Sys.time()
    cid <- .zarr_array(src, "cell_id")
    cs <- .zarr_array(src, "cell_summary")

    # Per-cell segmentation method: polygon_sets/1 is the CELL set (1:1
    # with cells; the nucleus set may hold multiple polygons per cell).
    # `method` indexes the `segmentation_methods` strings on the root
    # .zattrs.
    seg_method <- rep(NA_character_, nrow(cs))
    if (.zarr_exists(src, "polygon_sets/1/method")) {
        zattrs <- .zarr_attrs(src, "")
        methods_vec <- zattrs$segmentation_methods
        if (length(methods_vec) > 0L) {
            m_idx <- as.integer(.zarr_array(src, "polygon_sets/1/method"))
            ok <- !is.na(m_idx) & m_idx >= 0L & m_idx < length(methods_vec)
            seg_method[ok] <- methods_vec[m_idx[ok] + 1L]
        }
    }

    dt <- data.table::data.table(
        cell_id = .encode_xenium_id(cid[, 1L], cid[, 2L]),
        x_centroid = cs[, 1L],
        y_centroid = cs[, 2L],
        cell_area = cs[, 3L],
        nucleus_centroid_x = cs[, 4L],
        nucleus_centroid_y = cs[, 5L],
        nucleus_area = cs[, 6L],
        z_level = cs[, 7L],
        nucleus_count = as.integer(cs[, 8L]),
        segmentation_method = seg_method
    )
    arrow::write_parquet(dt, out_path,
        compression = .zarr_codec(compression))
    GiottoUtils::vmsg(sprintf(
        "  cells.parquet: %d rows in %.2fs", nrow(dt),
        as.numeric(Sys.time() - t0, units = "secs")
    ), .v = verbose)
    list(rows = nrow(dt), path = out_path)
}

# boundaries ####

# cells.zarr polygon_sets -> {cell,nucleus}_boundaries.parquet
# (cell_id, vertex_x, vertex_y; single file, cell-ordered).
#
# Per block of `block_size` polygons: ONE range read of the vertices
# array (only overlapping zarr chunks are decoded) and a fully vectorized
# expansion of the (n x 2*max_v) interleaved x,y rows into vertex runs.
.zarr_boundaries_to_parquet <- function(src, out_cell, out_nuc,
    cell_id_lookup, block_size = 100000L, compression = "zstd",
    verbose = NULL) {
    results <- list()
    for (kind in c("nucleus", "cell")) {
        pset_idx <- if (kind == "nucleus") "0" else "1"
        out_path <- if (kind == "nucleus") out_nuc else out_cell
        if (is.null(out_path)) next
        base <- paste0("polygon_sets/", pset_idx)
        if (!.zarr_exists(src, base)) {
            GiottoUtils::vmsg(sprintf(
                "  %s_boundaries: polygon_set %s missing, skipped",
                kind, pset_idx
            ), .v = verbose)
            next
        }
        t0 <- Sys.time()
        nv_full <- as.integer(.zarr_array(src, paste0(base, "/num_vertices")))
        ci_full <- as.integer(.zarr_array(src, paste0(base, "/cell_index")))
        n_poly <- length(nv_full)

        schema_b <- arrow::schema(
            cell_id = arrow::string(),
            vertex_x = arrow::float32(),
            vertex_y = arrow::float32()
        )
        sink <- arrow::FileOutputStream$create(out_path)
        props <- arrow::ParquetWriterProperties$create(
            column_names = schema_b$names,
            compression = .zarr_codec(compression)
        )
        writer <- arrow::ParquetFileWriter$create(
            schema = schema_b, sink = sink, properties = props
        )
        total_rows <- 0
        for (start in seq.int(1L, n_poly, by = block_size)) {
            end <- min(start + block_size - 1L, n_poly)
            nv <- nv_full[start:end]
            total_v <- sum(nv)
            if (total_v == 0L) next
            v <- .zarr_array(src, paste0(base, "/vertices"),
                range = c(start, end))
            ci <- ci_full[start:end]
            cids <- .encode_xenium_id(
                cell_id_lookup$prefix[ci + 1L],
                cell_id_lookup$suffix[ci + 1L]
            )
            # vectorized expansion: polygon i contributes its first
            # 2*nv[i] floats, interleaved x,y. rep.int/sequence handle
            # nv == 0 natively.
            poly <- rep.int(seq_along(nv), nv)
            k <- sequence(nv)
            tbl <- arrow::arrow_table(
                cell_id = rep.int(cids, nv),
                vertex_x = arrow::Array$create(
                    v[cbind(poly, 2L * k - 1L)], type = arrow::float32()),
                vertex_y = arrow::Array$create(
                    v[cbind(poly, 2L * k)], type = arrow::float32())
            )
            writer$WriteTable(tbl, chunk_size = nrow(tbl))
            total_rows <- total_rows + total_v
        }
        writer$Close()
        sink$close()
        GiottoUtils::vmsg(sprintf(
            "  %s_boundaries.parquet: %d rows in %.2fs", kind, total_rows,
            as.numeric(Sys.time() - t0, units = "secs")
        ), .v = verbose)
        results[[kind]] <- list(rows = total_rows, path = out_path)
    }
    results
}

# transcripts ####

# Canonical transcripts.parquet schema (matches the 10x-shipped file).
.xenium_tx_schema <- function() {
    arrow::schema(
        transcript_id = arrow::uint64(),
        cell_id = arrow::string(),
        overlaps_nucleus = arrow::uint8(),
        feature_name = arrow::string(),
        x_location = arrow::float32(),
        y_location = arrow::float32(),
        z_location = arrow::float32(),
        qv = arrow::float32(),
        fov_name = arrow::string(),
        nucleus_distance = arrow::float32(),
        codeword_index = arrow::int32()
    )
}

# Lookups for the transcript arrays, resolved per source. Format v4
# archives carry both on the root .zattrs: `gene_names` (gene_identity ->
# feature name, INCLUDING the UnassignedCodeword_* entries absent from
# gene_panel.json) and `fov_names` (the instrument's alphanumeric FOV
# codes). Fallbacks: the supplied panel lookup, and synthetic "FOV%03d".
.zarr_tx_lookups <- function(src, gene_lookup = character(0L)) {
    attrs <- .zarr_attrs(src, "")
    gl <- attrs$gene_names
    if (!length(gl)) gl <- gene_lookup
    list(
        gene_names = as.character(gl),
        fov_names = as.character(attrs$fov_names %||% character(0L))
    )
}

# Build the arrow table for one grid tile (no I/O beyond the reads).
# Pure per-tile kernel shared by the sequential and parallel paths.
# NULL when the tile is empty (or empty after the optional qv filter).
.xenium_tx_tile_table <- function(src, tile, gene_lookup,
    qv_threshold = NULL, fov_names = character(0L)) {
    pref <- paste0("grids/0/", tile)
    loc <- .zarr_array(src, paste0(pref, "/location"))
    if (!is.matrix(loc)) return(NULL)
    n <- nrow(loc)
    if (n == 0L) return(NULL)
    .col1 <- function(x) if (is.matrix(x)) x[, 1L] else x
    ge <- as.integer(.col1(.zarr_array(src, paste0(pref, "/gene_identity"))))
    cw <- as.integer(.col1(.zarr_array(src,
        paste0(pref, "/codeword_identity"))))
    qv_ <- as.numeric(.col1(.zarr_array(src, paste0(pref, "/quality_score"))))
    va <- as.integer(.col1(.zarr_array(src, paste0(pref, "/valid"))))
    id_ <- .zarr_array(src, paste0(pref, "/id"))

    keep <- (va == 1L)
    if (!is.null(qv_threshold)) keep <- keep & (qv_ > qv_threshold)
    if (!any(keep)) return(NULL)
    if (!all(keep)) {
        loc <- loc[keep, , drop = FALSE]
        ge <- ge[keep]
        cw <- cw[keep]
        qv_ <- qv_[keep]
        id_ <- id_[keep, , drop = FALSE]
        n <- sum(keep)
    }

    # 10x packs transcript_id as (2^16 + fov_index) in the high 32 bits
    # and the within-FOV decode counter in the low 32 bits — verified
    # against shipped transcripts.parquet (id column 1 = counter,
    # column 2 = fov index)
    fov_idx <- as.integer(id_[, 2L])
    tid <- bit64::as.integer64(65536 + fov_idx) * 4294967296 +
        bit64::as.integer64(id_[, 1L])
    fov_name <- if (length(fov_names) && all(fov_idx < length(fov_names))) {
        fov_names[fov_idx + 1L]
    } else {
        # no fov_names on the source: format each distinct index once
        ufov <- unique(fov_idx)
        sprintf("FOV%03d", ufov)[match(fov_idx, ufov)]
    }

    fname <- if (length(gene_lookup)) {
        idx1 <- ge + 1L
        out <- rep(NA_character_, length(idx1))
        ok <- !is.na(idx1) & idx1 >= 1L & idx1 <= length(gene_lookup)
        out[ok] <- gene_lookup[idx1[ok]]
        out
    } else {
        paste0("gene_", ge)
    }

    f32 <- function(x) arrow::Array$create(x, type = arrow::float32())
    arrow::arrow_table(
        transcript_id = arrow::Array$create(tid, type = arrow::uint64()),
        cell_id = rep("UNASSIGNED", n),
        overlaps_nucleus = arrow::Array$create(
            rep(0L, n), type = arrow::uint8()),
        feature_name = fname,
        x_location = f32(loc[, 1L]),
        y_location = f32(loc[, 2L]),
        z_location = f32(loc[, 3L]),
        qv = f32(qv_),
        fov_name = fov_name,
        nucleus_distance = f32(rep(NA_real_, n)),
        codeword_index = arrow::Array$create(cw, type = arrow::int32())
    )
}

# One worker: open an own zarr source (zip connections don't survive
# fork/serialization), process a tile subset, write one parquet shard.
# Tile tables are buffered and flushed at >= flush_rows so row groups are
# ~1-4M rows instead of one tiny group per tile. Top-level internal so it
# ships to workers as a namespace reference, never a closure.
#' @keywords internal
#' @noRd
.xenium_tx_worker <- function(tiles, zarr_path, gene_lookup, qv_threshold,
    part_path, flush_rows = 2000000L, compression = "zstd") {
    src <- .zarr_open(zarr_path)
    on.exit(.zarr_close(src), add = TRUE)
    lk <- .zarr_tx_lookups(src, gene_lookup)
    schema_tx <- .xenium_tx_schema()
    sink <- arrow::FileOutputStream$create(part_path)
    props <- arrow::ParquetWriterProperties$create(
        column_names = schema_tx$names,
        compression = .zarr_codec(compression)
    )
    writer <- arrow::ParquetFileWriter$create(
        schema = schema_tx, sink = sink, properties = props
    )
    on.exit({
        try(writer$Close(), silent = TRUE)
        try(sink$close(), silent = TRUE)
    }, add = TRUE)

    buf <- list()
    buf_rows <- 0L
    flush <- function() {
        if (!length(buf)) return(invisible(NULL))
        tbl <- if (length(buf) == 1L) buf[[1L]] else
            do.call(arrow::concat_tables, buf)
        writer$WriteTable(tbl, chunk_size = nrow(tbl))
        buf <<- list()
        buf_rows <<- 0L
        invisible(NULL)
    }
    n_total <- 0
    for (tile in tiles) {
        tbl <- .xenium_tx_tile_table(src, tile, lk$gene_names,
            qv_threshold, fov_names = lk$fov_names)
        if (is.null(tbl)) next
        buf[[length(buf) + 1L]] <- tbl
        buf_rows <- buf_rows + nrow(tbl)
        n_total <- n_total + nrow(tbl)
        if (buf_rows >= flush_rows) flush()
    }
    flush()
    n_total
}

# transcripts.zarr -> transcripts parquet. `workers = 1` writes a single
# out_path file; `workers > 1` writes a directory of part-*.parquet
# shards (read transparently by arrow::open_dataset). Scheduling follows
# .storewrite_h5_parallel: fork on unix, lapply_flex elsewhere, plain
# lapply when serial — one worker function for all three.
.zarr_transcripts_to_parquet <- function(src, out_path, gene_lookup,
    qv_threshold = NULL, workers = 1L, zarr_path = NULL,
    flush_rows = 2000000L, compression = "zstd", verbose = NULL) {
    t0 <- Sys.time()
    tiles <- .zarr_list(src, "grids/0", dirs_only = TRUE)
    if (!length(tiles)) {
        stop("[zarr] no transcript tiles under grids/0 — corrupt or ",
            "unsupported archive layout", call. = FALSE)
    }
    workers <- max(1L, as.integer(workers))
    workers <- min(workers, length(tiles))

    if (workers > 1L) {
        if (is.null(zarr_path)) {
            stop("[zarr] workers > 1 requires `zarr_path` so each worker ",
                "can open its own source", call. = FALSE)
        }
        if (file.exists(out_path) && !dir.exists(out_path)) unlink(out_path)
        dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
        # round-robin striping balances uneven tile sizes across workers
        assign <- ((seq_along(tiles) - 1L) %% workers) + 1L
        tile_groups <- split(tiles, assign)
        part_paths <- file.path(
            out_path, sprintf("part-%03d.parquet", seq_len(workers))
        )
        run_one <- function(i) {
            .xenium_tx_worker(
                tiles = tile_groups[[i]], zarr_path = zarr_path,
                gene_lookup = gene_lookup, qv_threshold = qv_threshold,
                part_path = part_paths[i], flush_rows = flush_rows,
                compression = compression
            )
        }
        counts <- if (.Platform$OS.type == "unix") {
            parallel::mclapply(seq_len(workers), run_one,
                mc.cores = workers, mc.preschedule = FALSE)
        } else {
            GiottoUtils::lapply_flex(seq_len(workers), run_one,
                cores = workers, future.seed = NULL)
        }
        errs <- vapply(counts, inherits, logical(1L), "try-error")
        if (any(errs)) {
            stop("[zarr] transcript worker failed: ",
                attr(counts[[which(errs)[1L]]], "condition")$message,
                call. = FALSE)
        }
        total_rows <- sum(unlist(counts))
        GiottoUtils::vmsg(sprintf(
            "  transcripts/ (dataset): %d rows in %.2fs (%d shards)",
            total_rows, as.numeric(Sys.time() - t0, units = "secs"), workers
        ), .v = verbose)
        return(list(rows = total_rows, path = out_path, n_parts = workers))
    }

    # serial: single file through the same worker kernel
    total_rows <- .xenium_tx_worker(
        tiles = tiles,
        zarr_path = if (!is.null(zarr_path)) zarr_path else src$path,
        gene_lookup = gene_lookup, qv_threshold = qv_threshold,
        part_path = out_path, flush_rows = flush_rows,
        compression = compression
    )
    GiottoUtils::vmsg(sprintf(
        "  transcripts.parquet: %d rows in %.2fs", total_rows,
        as.numeric(Sys.time() - t0, units = "secs")
    ), .v = verbose)
    list(rows = total_rows, path = out_path)
}
