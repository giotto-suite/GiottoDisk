# Atera v1 (atera-zarr-v6) archives: decoders, the direct-to-tile
# transcripts writer, the CSR expression input, the standalone converter,
# and importAteraDisk on the raw output folder. The fixture is a real
# 13-cell Atera run (annotation arrays stripped from the matrix) plus
# stubs of the archives detection must skip.

skip_if_no_zarr_deps()

av6 <- test_path("fixtures", "atera_v6")
av6_zip <- function(name) file.path(av6, paste0(name, ".zarr.zip"))

# transcripts written as a 10x-schema table, for reference
av6_tx_table <- function() {
    src <- GiottoDisk:::.zarr_open(av6_zip("transcripts"))
    on.exit(GiottoDisk:::.zarr_close(src))
    lk <- GiottoDisk:::.atera_tx_lookups(src)
    units <- GiottoDisk:::.atera_tx_units(src)
    tbl <- do.call(arrow::concat_tables, lapply(units, function(u) {
        GiottoDisk:::.atera_tx_unit_table(src, u, lk)
    }))
    data.table::as.data.table(as.data.frame(tbl))
}

test_that("half, int64 and vlen-utf8 decoders", {
    expect_identical(
        GiottoDisk:::.half_to_double(
            c(0L, 15360L, 49152L, 1L, 31743L, 31744L, 64512L)
        ),
        c(0, 1, -2, 2^-24, 65504, Inf, -Inf)
    )
    expect_true(is.nan(GiottoDisk:::.half_to_double(32256L)))

    i8 <- writeBin(as.integer(c(5L, 0L, -1L, -1L, 0L, 1L)), raw(),
        size = 4L, endian = "little")
    expect_identical(
        GiottoDisk:::.zarr_decode_values(i8,
            GiottoDisk:::.parse_zarr_dtype("<i8")),
        c(5, -1, 4294967296)
    )

    strs <- c("ATOX1", "", "CD19")
    vl <- c(writeBin(3L, raw(), size = 4L, endian = "little"),
        unlist(lapply(strs, function(s) c(
            writeBin(nchar(s, "bytes"), raw(), size = 4L, endian = "little"),
            charToRaw(s)))))
    expect_identical(GiottoDisk:::.zarr_decode_vlen_utf8(vl), strs)
})

test_that("each archive is detected on its own", {
    det <- function(name) {
        GiottoDisk:::.zarr_detect_product(av6_zip(name), name)
    }
    expect_identical(det("cells"), "xenium-zarr-v1")
    expect_identical(det("transcripts"), "atera-zarr-v6")
    expect_identical(det("cell_feature_matrix"), "atera-zarr-v6")
    expect_identical(detectZarrLayout(av6)$layout_version, "atera-zarr-v6")
})

test_that("v6 transcripts keep cell ids and gene names", {
    tx <- av6_tx_table()
    expect_identical(nrow(tx), 1023L)
    expect_identical(sum(tx$cell_id != "UNASSIGNED"), 674L)
    expect_false(anyDuplicated(tx$transcript_id) > 0L)
    # gene-sorted within the grid cell: gene_offset gives the same names
    src <- GiottoDisk:::.zarr_open(av6_zip("transcripts"))
    on.exit(GiottoDisk:::.zarr_close(src))
    go <- GiottoDisk:::.zarr_array(src, "grid/0,0/gene_offset")
    gn <- GiottoDisk:::.zarr_attrs(src, "")$gene_names
    expect_identical(tx$feature_name, rep(gn, go[, 2L] - go[, 1L]))
})

test_that("transcripts stream straight into balanced tile partitions", {
    src <- gDirSource(path = withr::local_tempdir())
    st <- sourceWrite(src, GiottoDisk:::tenxZarrTxInput(av6_zip("transcripts")),
        store_type = "parquetGeomTile", threshold = 100, qv_threshold = 30,
        workers = 1L, window_rows = 300L, buffer_rows = 500L)
    expect_s4_class(st, "parquetGeomTileStore")

    df <- as.data.frame(dplyr::collect(storeRead(st, omit_internals = FALSE)))
    tx <- av6_tx_table()
    expect_identical(nrow(df), sum(tx$qv >= 30))
    expect_setequal(df$transcript_id, as.character(tx$transcript_id[tx$qv >= 30]))
    # y flipped like the reader does
    expect_true(all(df$y_index <= 0))

    nr <- st@tiles@metadata$n_records
    expect_identical(length(nr), 16L) # s = ceiling(sqrt(1023 / 100))
    expect_lt(max(nr) / min(nr), 1.5)
    # row_index contiguous per tile, every point inside its tile bounds
    by_tile <- split(df$row_index, df$tile_index)
    expect_true(all(vapply(by_tile, function(r) {
        identical(sort(r), seq_along(r))
    }, logical(1L))))
    b <- st@tiles@bounds[df$tile_index, , drop = FALSE]
    expect_true(all(df$x_index >= b[, 1L] & df$x_index <= b[, 2L] &
        df$y_index >= b[, 3L] & df$y_index <= b[, 4L]))

    # crop prunes tiles and returns exactly the points inside
    e <- terra::ext(st)
    q <- terra::ext(e$xmin, e$xmin + (e$xmax - e$xmin) / 3,
        e$ymin, e$ymin + (e$ymax - e$ymin) / 3)
    cr <- crop(st, q)
    expect_gt(length(cr@tile_filter), 0L)
    expect_lt(length(cr@tile_filter), 16L)
    expect_identical(nrow(dplyr::collect(storeRead(cr))), sum(
        df$x_index >= q$xmin & df$x_index <= q$xmax &
            df$y_index >= q$ymin & df$y_index <= q$ymax))

    # the same store from parallel workers
    st2 <- sourceWrite(src, GiottoDisk:::tenxZarrTxInput(av6_zip("transcripts")),
        store_type = "parquetGeomTile", threshold = 100, qv_threshold = 30,
        workers = 2L, window_rows = 300L, buffer_rows = 500L)
    expect_identical(st2@tiles@metadata$n_records, nr)
})

test_that("the CSR matrix input streams cell-major and matches transcripts", {
    inp <- tenxZarrInput(av6_zip("cell_feature_matrix"))
    expect_s4_class(inp, "tenxZarrCsrInput")
    expect_identical(c(inp@n_cells, inp@n_genes), c(13L, 132L))
    expect_setequal(unique(inp@feat_types), c("Gene Expression",
        "Negative Control Codeword", "Negative Control Probe",
        "Unassigned Codeword"))
    withr::local_options(giottodisk.h5_batch_rows = 50)
    itr <- storeRead(inp)
    batches <- list()
    repeat {
        x <- itr$next_batch()
        if (is.null(x)) break
        batches[[length(batches) + 1L]] <- x
    }
    itr$close()
    expect_gt(length(batches), 1L)
    tri <- data.table::rbindlist(batches)
    expect_identical(nrow(tri), 315L)
    expect_false(is.unsorted(tri$row_id))

    src <- GiottoDisk:::.zarr_open(av6_zip("cell_feature_matrix"))
    on.exit(GiottoDisk:::.zarr_close(src))
    totals <- GiottoDisk:::.zarr_array(src, "obs/total_counts")
    expect_equal(tri[, sum(value), by = row_id][order(row_id)]$V1,
        as.numeric(totals))

    # every (cell, gene) count equals the transcripts assigned to it
    tx <- av6_tx_table()
    agg <- tx[cell_id != "UNASSIGNED", .N, by = .(cell_id, feature_name)]
    tri[, `:=`(cell_id = inp@cell_ids[row_id],
        feature_name = inp@feat_ids[col_id])]
    m <- merge(agg, tri, by = c("cell_id", "feature_name"), all = TRUE)
    expect_false(anyNA(m$N) || anyNA(m$value))
    expect_identical(as.numeric(m$N), m$value)
})

test_that("xeniumZarrToParquet converts a v6 folder", {
    out <- withr::local_tempdir()
    res <- xeniumZarrToParquet(av6, out, workers = 1L, verbose = FALSE)
    expect_identical(nrow(arrow::read_parquet(res$cells)), 13L)
    expect_identical(nrow(arrow::read_parquet(res$transcripts)), 1023L)
    expect_identical(nrow(arrow::read_parquet(res$cell_feature_matrix)), 315L)
    expect_true(file.exists(file.path(out, "experiment.xenium")))
    man <- jsonlite::fromJSON(file.path(out, "conversion_manifest.json"))
    expect_identical(man$layout_version, "atera-zarr-v6")
})

test_that("importAteraDisk builds a gobject from the raw Atera folder", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("tilework")
    withr::local_options(giottodisk.artifact_dump = withr::local_tempdir())
    withr::local_options(tilework.warn_sequential = FALSE)
    withr::local_options(giotto.warn_sequential = FALSE)
    withr::local_options(giottodisk.par_workers = 1L)

    rdr <- importAteraDisk(av6, backend = withr::local_tempdir(),
        qv_threshold = 20)
    expect_s4_class(rdr, "AteraDiskReader")
    expect_false(methods::is(rdr, "XeniumReader"))
    # the exact archives, not binned_transcripts / csc_cell_feature_matrix
    expect_identical(basename(rdr@paths$tx_path), "transcripts.zarr.zip")
    expect_identical(basename(rdr@paths$expr_path),
        "cell_feature_matrix.zarr.zip")
    expect_equal(rdr$micron,
        jsonlite::fromJSON(file.path(av6, "experiment.spatial"))$pixel_size)

    g <- suppressWarnings(
        rdr$create_gobject(load_images = NULL, verbose = FALSE)
    )
    expect_s4_class(g, "giotto")
    expect_setequal(
        GiottoClass::list_spatial_info_names(g), c("cell", "nucleus")
    )
    pts <- GiottoClass::getFeatureInfo(g, feat_type = "rna",
        return_giottoPoints = TRUE)[]
    expect_s4_class(pts, "parquetGeomTileStore")
    expect_identical(nrow(dplyr::collect(storeRead(pts))), 1011L)
    pe <- GiottoClass::getExpression(g, spat_unit = "cell",
        feat_type = "rna", values = "raw")[]
    expect_length(pe@cell_ids, 13L)
    expect_identical(nrow(dplyr::collect(storeRead(pe))), 308L)
})

test_that("create_gobject loads the morphology_2d images by default", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("tilework")
    withr::local_options(giottodisk.artifact_dump = withr::local_tempdir())
    withr::local_options(tilework.warn_sequential = FALSE)
    withr::local_options(giotto.warn_sequential = FALSE)
    withr::local_options(giottodisk.par_workers = 1L)
    # image conversion writes next to the images: work on a copy
    dir <- withr::local_tempdir()
    file.copy(list.files(av6, full.names = TRUE), dir, recursive = TRUE)

    rdr <- importAteraDisk(dir, backend = withr::local_tempdir())
    expect_identical(names(rdr@paths$img_paths),
        c("dapi", "atp1a1_cd45_e-cadherin", "18s", "alphasma_vimentin"))
    g <- suppressWarnings(rdr$create_gobject(load_transcripts = FALSE,
        verbose = FALSE))
    expect_setequal(GiottoClass::list_images_names(g),
        names(rdr@paths$img_paths))
})

test_that("Giotto's createGiottoAteraObject builds the disk object", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("tilework")
    skip_if_not(exists("createGiottoAteraObject", asNamespace("Giotto")))
    withr::local_options(giottodisk.artifact_dump = withr::local_tempdir())
    withr::local_options(tilework.warn_sequential = FALSE)
    withr::local_options(giotto.warn_sequential = FALSE)
    withr::local_options(giottodisk.par_workers = 1L)
    dir <- withr::local_tempdir()
    file.copy(list.files(av6, full.names = TRUE), dir, recursive = TRUE)

    # the wrapper passes the 10x reader arguments (bounds as a named list,
    # load_images = "focus", load_aligned_images, path overrides)
    g <- suppressWarnings(Giotto::createGiottoAteraObject(atera_dir = dir,
        backend = withr::local_tempdir(), verbose = FALSE))
    expect_s4_class(g, "giotto")
    expect_setequal(
        GiottoClass::list_spatial_info_names(g), c("cell", "nucleus")
    )
    expect_length(GiottoClass::list_images_names(g), 4L)
    pe <- GiottoClass::getExpression(g, spat_unit = "cell",
        feat_type = "rna", values = "raw")[]
    expect_s4_class(pe, "parquetExprStore")
})
