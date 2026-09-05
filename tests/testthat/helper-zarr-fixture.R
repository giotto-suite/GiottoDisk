# Builds a synthetic Xenium-layout zarr fixture at test time (no binary
# fixtures in the repo — same convention as the GEF tests). Ground truth
# is the tables the fixture was built from, never a second run of the
# same read path.
#
# Layout matches Xenium Onboard Analysis zarr exports: cells.zarr.zip
# (cell_id, cell_summary, polygon_sets/{0,1}), transcripts.zarr.zip
# (grids/0/<tile>/ arrays), cell_feature_matrix.zarr.zip (CSC-by-feature
# cell_features), gene_panel.json. Edge cases baked in: uint32 values
# > 2^31, uint64 indptr, multi-chunk 1D arrays, trailing-dim-chunked 2D
# vertices, a zero-vertex polygon, an all-invalid transcript tile, an
# aggregate_gene feature, duplicate feature symbols.
#
# Most arrays are hand-encoded with `"compressor": null` (raw chunks) so
# every dtype path is exercised precisely; cell_summary is written
# through Rarr with blosc so the real decompression path runs too.

skip_if_no_zarr_deps <- function() {
    testthat::skip_if_not_installed("Rarr")
    testthat::skip_if_not_installed("zip")
}

# encode a value vector as raw little-endian bytes per zarr dtype
.zf_encode <- function(x, dtype) {
    switch(dtype,
        "<f4" = writeBin(as.double(x), raw(), size = 4L, endian = "little"),
        "<f8" = writeBin(as.double(x), raw(), size = 8L, endian = "little"),
        "<u1" = ,
        "|u1" = writeBin(as.integer(x), raw(), size = 1L, endian = "little"),
        "<u2" = writeBin(as.integer(x), raw(), size = 2L, endian = "little"),
        "<i4" = writeBin(as.integer(x), raw(), size = 4L, endian = "little"),
        "<u4" = {
            v <- as.double(x)
            v[v >= 2^31] <- v[v >= 2^31] - 2^32
            writeBin(as.integer(v), raw(), size = 4L, endian = "little")
        },
        "<u8" = {
            lo <- as.double(x) %% 2^32
            hi <- as.double(x) %/% 2^32
            lo[lo >= 2^31] <- lo[lo >= 2^31] - 2^32
            inter <- as.integer(rbind(lo, hi))
            writeBin(inter, raw(), size = 4L, endian = "little")
        },
        stop("unsupported fixture dtype: ", dtype)
    )
}

.zf_write_json <- function(root, rel, obj) {
    path <- file.path(root, rel)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(
        jsonlite::toJSON(obj, auto_unbox = TRUE, null = "null", digits = NA),
        path
    )
}

.zf_write_group <- function(root, prefix = "") {
    rel <- if (nzchar(prefix)) file.path(prefix, ".zgroup") else ".zgroup"
    .zf_write_json(root, rel, list(zarr_format = 2L))
}

# Write a 1D vector or 2D matrix as a raw-chunk (compressor null) zarr
# array in C order, padding edge chunks with `fill` per the spec.
.zf_write_array <- function(root, prefix, x, dtype, chunks, fill = 0) {
    adir <- file.path(root, prefix)
    dir.create(adir, recursive = TRUE, showWarnings = FALSE)
    shape <- if (is.matrix(x)) dim(x) else length(x)
    chunks <- as.integer(chunks)
    .zf_write_json(root, file.path(prefix, ".zarray"), list(
        chunks = as.list(chunks),
        compressor = NULL,
        dtype = dtype,
        fill_value = fill,
        filters = NULL,
        order = "C",
        shape = as.list(as.integer(shape)),
        zarr_format = 2L
    ))
    if (length(shape) == 1L) {
        n_chunks <- ceiling(shape / chunks)
        for (ci in seq_len(n_chunks) - 1L) {
            s <- ci * chunks[1L] + 1L
            e <- min(s + chunks[1L] - 1L, shape[1L])
            vals <- rep(fill, chunks[1L])
            vals[seq_len(e - s + 1L)] <- x[s:e]
            writeBin(.zf_encode(vals, dtype), file.path(adir, ci))
        }
    } else {
        for (ci in seq_len(ceiling(shape[1L] / chunks[1L])) - 1L) {
            for (cj in seq_len(ceiling(shape[2L] / chunks[2L])) - 1L) {
                block <- matrix(fill, nrow = chunks[1L], ncol = chunks[2L])
                rs <- ci * chunks[1L] + 1L
                re <- min(rs + chunks[1L] - 1L, shape[1L])
                cs <- cj * chunks[2L] + 1L
                ce <- min(cs + chunks[2L] - 1L, shape[2L])
                block[seq_len(re - rs + 1L), seq_len(ce - cs + 1L)] <-
                    x[rs:re, cs:ce, drop = FALSE]
                # C order = row-major
                writeBin(.zf_encode(as.vector(t(block)), dtype),
                    file.path(adir, paste0(ci, ".", cj)))
            }
        }
    }
    invisible(NULL)
}

# zip a fixture tree with STORED entries (the seek-based in-place reader
# requires method 0, matching real Xenium archives)
.zf_zip <- function(root, zipfile) {
    files <- list.files(root, recursive = TRUE, all.files = TRUE,
        no.. = TRUE)
    if (file.exists(zipfile)) unlink(zipfile)
    zip::zip(zipfile, files = files, root = root,
        compression_level = 0, include_directories = FALSE,
        mode = "mirror")
    invisible(zipfile)
}

# 10x a-p barcode encoding per the documented spec (uint32 -> 8 hex
# chars -> a-p alphabet), written out here so fixture truth does not
# come from the code under test
.zf_barcode <- function(prefix, suffix = 1L) {
    hi <- as.integer(prefix %/% 65536)
    lo <- as.integer(prefix %% 65536)
    hex <- sprintf("%04x%04x", hi, lo)
    paste0(chartr("0123456789abcdef", "abcdefghijklmnop", hex), "-", suffix)
}

# Build the full mini fixture. Returns list(dir, paths, truth).
make_zarr_fixture <- function(dir = file.path(tempdir(),
    paste0("mini_xenium_", as.integer(stats::runif(1L, 1, 1e8))))) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    n_cells <- 12L
    cid_prefix <- c(101, 102, 103, 104, 105, 3000000000, 107, 108, 109,
        110, 111, 112) # one > 2^31: the uint32 edge
    cell_ids <- .zf_barcode(cid_prefix)

    # ---- cells.zarr ----
    croot <- file.path(dir, "_build_cells")
    dir.create(croot)
    .zf_write_group(croot)
    seg_methods <- c("dapi", "boundary", "interior")
    .zf_write_json(croot, ".zattrs",
        list(segmentation_methods = seg_methods))
    cid <- cbind(cid_prefix, rep(1, n_cells))
    .zf_write_array(croot, "cell_id", cid, "<u4", chunks = c(7L, 1L))
    cs <- cbind(
        x_centroid = seq_len(n_cells) * 10,
        y_centroid = seq_len(n_cells) * 20,
        cell_area = seq_len(n_cells) + 0.5,
        nucleus_centroid_x = seq_len(n_cells) * 10 + 1,
        nucleus_centroid_y = seq_len(n_cells) * 20 + 1,
        nucleus_area = seq_len(n_cells) / 2,
        z_level = rep(1, n_cells),
        nucleus_count = rep(c(1, 2), 6)
    )
    # blosc-compressed via Rarr so the real decompression path runs.
    # F order: Rarr's writer round-trips F correctly (its C-order write
    # is transposed on disk); the reader under test handles both orders.
    p_cs <- file.path(croot, "cell_summary")
    Rarr::create_empty_zarr_array(p_cs, dim = dim(cs),
        chunk_dim = c(7L, 8L), data_type = "double", order = "F",
        compressor = Rarr::use_blosc(), fill_value = 0)
    Rarr::update_zarr_array(p_cs, cs,
        index = list(seq_len(nrow(cs)), seq_len(ncol(cs))))

    # polygon sets: 0 = nucleus (has a zero-vertex polygon), 1 = cell.
    # vertices (12 x 10) f4, chunked (5, 4): multi-chunk on BOTH axes so
    # range reads cross chunk boundaries on the trailing dim too.
    mk_polyset <- function(pset, nv, with_method = FALSE) {
        base <- file.path("polygon_sets", pset)
        .zf_write_array(croot, file.path(base, "num_vertices"),
            nv, "<u4", chunks = 5L)
        .zf_write_array(croot, file.path(base, "cell_index"),
            seq_len(n_cells) - 1L, "<u4", chunks = 5L)
        v <- matrix(0, nrow = n_cells, ncol = 10L)
        for (i in seq_len(n_cells)) {
            k <- nv[i]
            if (k > 0L) {
                v[i, seq_len(2L * k)] <- as.vector(rbind(
                    i * 10 + seq_len(k) + 0.5, # x, exactly f4-representable
                    i * 100 + seq_len(k) + 0.25 # y
                ))
            }
        }
        .zf_write_array(croot, file.path(base, "vertices"),
            v, "<f4", chunks = c(5L, 4L))
        if (with_method) {
            .zf_write_array(croot, file.path(base, "method"),
                rep(0:2, 4), "<u4", chunks = 5L)
        }
        # truth table: one row per vertex, polygons in order
        poly <- rep.int(seq_len(n_cells), nv)
        k <- sequence(nv)
        data.table::data.table(
            cell_id = rep.int(cell_ids, nv),
            vertex_x = poly * 10 + k + 0.5,
            vertex_y = poly * 100 + k + 0.25
        )
    }
    nv_nuc <- c(3L, 3L, 0L, 4L, 3L, 4L, 5L, 3L, 3L, 4L, 4L, 3L)
    nv_cell <- c(3L, 4L, 5L, 3L, 4L, 5L, 3L, 4L, 5L, 3L, 4L, 4L)
    truth_nuc <- mk_polyset("0", nv_nuc)
    truth_cell <- mk_polyset("1", nv_cell, with_method = TRUE)
    cells_zip <- file.path(dir, "cells.zarr.zip")
    .zf_zip(croot, cells_zip)
    unlink(croot, recursive = TRUE)

    truth_cells <- data.table::data.table(
        cell_id = cell_ids,
        x_centroid = cs[, 1L], y_centroid = cs[, 2L],
        cell_area = cs[, 3L],
        nucleus_centroid_x = cs[, 4L], nucleus_centroid_y = cs[, 5L],
        nucleus_area = cs[, 6L], z_level = cs[, 7L],
        nucleus_count = as.integer(cs[, 8L]),
        segmentation_method = seg_methods[rep(0:2, 4) + 1L]
    )

    # ---- gene panel ----
    panel_names <- c("GA", "GB", "GC", "GD",
        "NegControlProbe_1", "UnassignedCodeword_0001")
    .zf_write_json(dir, "gene_panel.json", list(
        payload = list(targets = lapply(panel_names, function(nm) {
            list(type = list(
                descriptor = "gene",
                data = list(id = paste0("ID_", nm), name = nm)
            ))
        }))
    ))

    # ---- transcripts.zarr ----
    troot <- file.path(dir, "_build_tx")
    dir.create(troot)
    .zf_write_group(troot)
    .zf_write_group(troot, "grids")
    .zf_write_group(troot, "grids/0")
    tx_truth <- list()
    tid_counter <- 0L
    mk_tile <- function(tile, n, fov, valid, qv) {
        pref <- file.path("grids/0", tile)
        loc <- cbind(
            x = seq_len(n) + tid_counter + 0.5,
            y = seq_len(n) * 2 + tid_counter + 0.5,
            z = rep(1.5, n)
        )
        ge <- (seq_len(n) + tid_counter) %% length(panel_names) # 0-based
        cw <- seq_len(n) + 100L
        # counters above 2^31 exercise the uint32 decode edge in the low
        # word; the packed id (2^16 + fov) * 2^32 + counter stays far
        # inside int64
        prefix32 <- 2999999990 + seq_len(n) + tid_counter
        .zf_write_array(troot, file.path(pref, "location"),
            loc, "<f4", chunks = c(max(4L, n), 3L))
        .zf_write_array(troot, file.path(pref, "gene_identity"),
            matrix(ge, ncol = 1L), "<u2", chunks = c(7L, 1L))
        .zf_write_array(troot, file.path(pref, "codeword_identity"),
            matrix(cw, ncol = 1L), "<u4", chunks = c(7L, 1L))
        .zf_write_array(troot, file.path(pref, "quality_score"),
            matrix(qv, ncol = 1L), "<f4", chunks = c(7L, 1L))
        .zf_write_array(troot, file.path(pref, "valid"),
            matrix(as.integer(valid), ncol = 1L), "|u1",
            chunks = c(7L, 1L))
        .zf_write_array(troot, file.path(pref, "id"),
            cbind(prefix32, rep(fov, n)), "<u4", chunks = c(7L, 2L))
        tid_counter <<- tid_counter + n
        keep <- valid == 1L
        data.table::data.table(
            # 10x id packing: (2^16 + fov) in the high word, counter low
            transcript_id = bit64::as.integer64(65536 + fov) *
                4294967296 + bit64::as.integer64(prefix32[keep]),
            feature_name = panel_names[ge[keep] + 1L],
            x_location = loc[keep, 1L],
            y_location = loc[keep, 2L],
            z_location = loc[keep, 3L],
            qv = qv[keep],
            fov_name = sprintf("FOV%03d", fov),
            codeword_index = as.integer(cw[keep])
        )
    }
    tx_truth[["0,0"]] <- mk_tile("0,0", 20L, fov = 1L,
        valid = rep(c(1L, 1L, 1L, 0L), 5), qv = rep(c(30, 10), 10))
    tx_truth[["0,1"]] <- mk_tile("0,1", 15L, fov = 2L,
        valid = rep(1L, 15), qv = seq(5, 40, length.out = 15))
    # all-invalid tile: must contribute zero rows
    tx_truth[["1,1"]] <- mk_tile("1,1", 4L, fov = 3L,
        valid = rep(0L, 4), qv = rep(30, 4))
    tx_zip <- file.path(dir, "transcripts.zarr.zip")
    .zf_zip(troot, tx_zip)
    unlink(troot, recursive = TRUE)
    truth_tx <- data.table::rbindlist(tx_truth)

    # ---- cell_feature_matrix.zarr ----
    # 7 features: 4 genes, 1 neg control probe,
    # 1 unassigned codeword, 1 aggregate_gene (colsums; dropped by input).
    # Keys match the panel names, as in real exports.
    feat_ids_ens <- c("ENSG01", "ENSG02", "ENSG03", "ENSG04",
        "NegControlProbe_1", "UnassignedCodeword_0001", "Total transcripts")
    feat_keys <- c("GA", "GB", "GC", "GD",
        "NegControlProbe_1", "UnassignedCodeword_0001", "Total transcripts")
    feat_types <- c("gene", "gene", "gene", "gene",
        "negative_control_probe", "unassigned_codeword", "aggregate_gene")
    set.seed(42)
    m <- matrix(0, nrow = 6L, ncol = n_cells) # kept features x cells
    m[cbind(
        sample(1:6, 34, replace = TRUE),
        sample(seq_len(n_cells), 34, replace = TRUE)
    )] <- sample(1:20, 34, replace = TRUE)
    m_all <- rbind(m, colSums(m)) # + aggregate row
    # CSC by feature: cells ascending within each feature
    indices <- integer(0)
    values <- integer(0)
    indptr <- 0
    for (j in seq_len(nrow(m_all))) {
        nz <- which(m_all[j, ] > 0)
        indices <- c(indices, nz - 1L)
        values <- c(values, m_all[j, nz])
        indptr <- c(indptr, length(indices))
    }
    froot <- file.path(dir, "_build_cfm")
    dir.create(froot)
    .zf_write_group(froot)
    .zf_write_group(froot, "cell_features")
    .zf_write_json(froot, "cell_features/.zattrs", list(
        feature_ids = feat_ids_ens,
        feature_keys = feat_keys,
        feature_types = feat_types,
        major_version = 1L, minor_version = 0L,
        number_cells = n_cells,
        number_features = length(feat_ids_ens)
    ))
    # uint64 indptr (the Atera-scale dtype) + multi-chunk indices/data
    .zf_write_array(froot, "cell_features/indptr",
        indptr, "<u8", chunks = length(indptr))
    .zf_write_array(froot, "cell_features/indices",
        indices, "<u4", chunks = 7L)
    .zf_write_array(froot, "cell_features/data",
        as.integer(values), "<u4", chunks = 7L)
    .zf_write_array(froot, "cell_features/cell_id",
        cid, "<u4", chunks = c(7L, 1L))
    cfm_zip <- file.path(dir, "cell_feature_matrix.zarr.zip")
    .zf_zip(froot, cfm_zip)
    unlink(froot, recursive = TRUE)

    list(
        dir = dir,
        paths = list(
            cells = cells_zip,
            transcripts = tx_zip,
            cell_feature_matrix = cfm_zip
        ),
        truth = list(
            cell_ids = cell_ids,
            cells = truth_cells,
            cell_boundaries = truth_cell,
            nucleus_boundaries = truth_nuc,
            nv_cell = nv_cell,
            nv_nucleus = nv_nuc,
            transcripts = truth_tx,
            cfm = m, # kept features x cells (aggregate dropped)
            cfm_feat_keys = feat_keys[1:6],
            cfm_feat_ids = feat_ids_ens[1:6],
            cfm_feat_types = feat_types[1:6],
            panel = panel_names
        )
    )
}
