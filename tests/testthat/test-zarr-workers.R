# Converter workers (cells / boundaries / transcripts) against the
# fixture's ground-truth tables.

skip_if_no_zarr_deps()

fx <- make_zarr_fixture()

test_that(".encode_xenium_id matches the documented a-p encoding", {
    # ffffffff -> pppppppp; 0 -> aaaaaaaa; and a value above 2^31
    expect_identical(
        GiottoDisk:::.encode_xenium_id(4294967295, 1L), "pppppppp-1"
    )
    expect_identical(GiottoDisk:::.encode_xenium_id(0, 1L), "aaaaaaaa-1")
    expect_identical(
        GiottoDisk:::.encode_xenium_id(3000000000, 2L),
        .zf_barcode(3000000000, 2L)
    )
})

test_that("cells worker writes the 10x cellmeta schema", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cells)
    on.exit(GiottoDisk:::.zarr_close(src))
    out <- withr::local_tempfile(fileext = ".parquet")
    r <- GiottoDisk:::.zarr_cells_to_parquet(src, out, verbose = FALSE)
    expect_equal(r$rows, 12)
    dt <- as.data.frame(arrow::read_parquet(out))
    truth <- as.data.frame(fx$truth$cells)
    expect_identical(colnames(dt), colnames(truth))
    expect_identical(dt$cell_id, truth$cell_id)
    expect_equal(dt$x_centroid, truth$x_centroid)
    expect_identical(dt$nucleus_count, truth$nucleus_count)
    expect_identical(dt$segmentation_method, truth$segmentation_method)
})

test_that("boundaries worker is exact, cell-ordered, single-file", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cells)
    on.exit(GiottoDisk:::.zarr_close(src))
    cid <- GiottoDisk:::.zarr_array(src, "cell_id")
    oc <- withr::local_tempfile(fileext = ".parquet")
    on_ <- withr::local_tempfile(fileext = ".parquet")
    # block_size smaller than the polygon count so the block loop and the
    # range reads actually iterate
    GiottoDisk:::.zarr_boundaries_to_parquet(
        src, oc, on_,
        cell_id_lookup = list(prefix = cid[, 1L], suffix = cid[, 2L]),
        block_size = 5L, verbose = FALSE
    )
    bc <- as.data.frame(arrow::read_parquet(oc))
    bn <- as.data.frame(arrow::read_parquet(on_))
    truth_c <- as.data.frame(fx$truth$cell_boundaries)
    truth_n <- as.data.frame(fx$truth$nucleus_boundaries)
    # exact, in order -- includes the zero-vertex nucleus polygon (absent
    # from the output but not shifting later polygons)
    expect_identical(bc$cell_id, truth_c$cell_id)
    expect_equal(bc$vertex_x, truth_c$vertex_x)
    expect_equal(bc$vertex_y, truth_c$vertex_y)
    expect_identical(bn$cell_id, truth_n$cell_id)
    expect_equal(bn$vertex_x, truth_n$vertex_x)
    expect_identical(nrow(bn), sum(fx$truth$nv_nucleus))
})

test_that("transcripts worker matches truth; all-invalid tile drops out", {
    src <- GiottoDisk:::.zarr_open(fx$paths$transcripts)
    on.exit(GiottoDisk:::.zarr_close(src))
    gl <- GiottoDisk:::.load_gene_panel_lookup(fx$dir)
    expect_identical(gl, fx$truth$panel)
    out <- withr::local_tempfile(fileext = ".parquet")
    r <- GiottoDisk:::.zarr_transcripts_to_parquet(
        src, out, gene_lookup = gl, verbose = FALSE
    )
    truth <- as.data.frame(fx$truth$transcripts)
    expect_equal(r$rows, nrow(truth)) # tile "1,1" contributed nothing
    dt <- as.data.frame(arrow::read_parquet(out))
    o1 <- order(suppressWarnings(as.numeric(dt$transcript_id)))
    o2 <- order(suppressWarnings(as.numeric(truth$transcript_id)))
    expect_identical(dt$feature_name[o1], truth$feature_name[o2])
    expect_equal(dt$x_location[o1], truth$x_location[o2])
    expect_equal(dt$qv[o1], truth$qv[o2])
    expect_identical(dt$fov_name[o1], truth$fov_name[o2])
    expect_identical(dt$codeword_index[o1], truth$codeword_index[o2])
    # not derivable from zarr: placeholders
    expect_true(all(dt$cell_id == "UNASSIGNED"))
    expect_true(all(is.na(dt$nucleus_distance)))
})

test_that("qv filter at conversion time is optional and exact", {
    src <- GiottoDisk:::.zarr_open(fx$paths$transcripts)
    on.exit(GiottoDisk:::.zarr_close(src))
    gl <- GiottoDisk:::.load_gene_panel_lookup(fx$dir)
    out <- withr::local_tempfile(fileext = ".parquet")
    GiottoDisk:::.zarr_transcripts_to_parquet(
        src, out, gene_lookup = gl, qv_threshold = 20, verbose = FALSE
    )
    dt <- as.data.frame(arrow::read_parquet(out))
    truth <- fx$truth$transcripts
    expect_identical(nrow(dt), sum(truth$qv > 20))
    expect_true(all(dt$qv > 20))
})

test_that("parallel transcripts equal serial (sorted comparison)", {
    skip_on_os("windows") # fork path; lapply_flex needs a future plan
    src <- GiottoDisk:::.zarr_open(fx$paths$transcripts)
    gl <- GiottoDisk:::.load_gene_panel_lookup(fx$dir)
    out1 <- withr::local_tempfile(fileext = ".parquet")
    GiottoDisk:::.zarr_transcripts_to_parquet(
        src, out1, gene_lookup = gl, verbose = FALSE
    )
    GiottoDisk:::.zarr_close(src)
    src <- GiottoDisk:::.zarr_open(fx$paths$transcripts)
    on.exit(GiottoDisk:::.zarr_close(src))
    outd <- withr::local_tempdir()
    r <- GiottoDisk:::.zarr_transcripts_to_parquet(
        src, outd, gene_lookup = gl, workers = 2L,
        zarr_path = fx$paths$transcripts, verbose = FALSE
    )
    expect_identical(r$n_parts, 2L)
    d1 <- as.data.frame(arrow::read_parquet(out1))
    d2 <- as.data.frame(dplyr::collect(arrow::open_dataset(outd)))
    d1 <- d1[order(suppressWarnings(as.numeric(d1$transcript_id))), ]
    d2 <- d2[order(suppressWarnings(as.numeric(d2$transcript_id))), ]
    rownames(d1) <- rownames(d2) <- NULL
    expect_equal(d1, d2)
})
