# Zarr source + array layer over synthetic fixtures built at test time
# (helper-zarr-fixture.R). Ground truth is the values the fixture was
# written from.

skip_if_no_zarr_deps()

fx <- make_zarr_fixture()

test_that("dtype parsing covers the Xenium dtypes", {
    expect_identical(
        GiottoDisk:::.parse_zarr_dtype("<u4"),
        list(kind = "u", size = 4L, endian = "little")
    )
    expect_identical(GiottoDisk:::.parse_zarr_dtype(">f8")$endian, "big")
    expect_identical(GiottoDisk:::.parse_zarr_dtype("|u1")$size, 1L)
})

test_that("zip source reads uint32 above 2^31 exactly", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cells)
    on.exit(GiottoDisk:::.zarr_close(src))
    cid <- GiottoDisk:::.zarr_array(src, "cell_id")
    expect_identical(cid[6L, 1L], 3000000000)
})

test_that("blosc-compressed chunks decode (Rarr internal path)", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cells)
    on.exit(GiottoDisk:::.zarr_close(src))
    cs <- GiottoDisk:::.zarr_array(src, "cell_summary")
    expect_equal(cs[, 1L], seq_len(12L) * 10)
    expect_equal(cs[, 8L], rep(c(1, 2), 6))
})

test_that("zip and unzipped-dir sources read identically", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cells)
    td <- withr::local_tempdir()
    utils::unzip(fx$paths$cells, exdir = td)
    src2 <- GiottoDisk:::.zarr_open(td)
    on.exit({
        GiottoDisk:::.zarr_close(src)
        GiottoDisk:::.zarr_close(src2)
    })
    for (pfx in c("cell_id", "polygon_sets/1/vertices",
        "polygon_sets/0/num_vertices")) {
        expect_identical(
            GiottoDisk:::.zarr_array(src, pfx),
            GiottoDisk:::.zarr_array(src2, pfx)
        )
    }
    expect_setequal(
        GiottoDisk:::.zarr_list(src, "polygon_sets", dirs_only = TRUE),
        GiottoDisk:::.zarr_list(src2, "polygon_sets", dirs_only = TRUE)
    )
})

test_that("range read equals full-read-then-subset across chunk bounds", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cells)
    on.exit(GiottoDisk:::.zarr_close(src))
    # 2D, chunked (5, 4) on a (12, 10) array: ranges crossing row-chunk
    # boundaries with column chunks always assembled in full
    v <- GiottoDisk:::.zarr_array(src, "polygon_sets/1/vertices")
    for (r in list(c(1L, 12L), c(4L, 9L), c(5L, 6L), c(11L, 12L),
        c(1L, 1L))) {
        expect_identical(
            GiottoDisk:::.zarr_array(src, "polygon_sets/1/vertices",
                range = r),
            v[r[1L]:r[2L], , drop = FALSE]
        )
    }
    # 1D multi-chunk
    nv <- GiottoDisk:::.zarr_array(src, "polygon_sets/1/num_vertices")
    expect_identical(
        GiottoDisk:::.zarr_array(src, "polygon_sets/1/num_vertices",
            range = c(4L, 11L)),
        nv[4:11]
    )
    expect_error(
        GiottoDisk:::.zarr_array(src, "polygon_sets/1/num_vertices",
            range = c(0L, 5L)),
        "bad range"
    )
})

test_that("uint64 arrays decode as exact doubles", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cell_feature_matrix)
    on.exit(GiottoDisk:::.zarr_close(src))
    meta <- GiottoDisk:::.zarr_meta(src, "cell_features/indptr")
    expect_identical(meta$dtype, "<u8")
    indptr <- GiottoDisk:::.zarr_array(src, "cell_features/indptr")
    expect_type(indptr, "double")
    expect_identical(indptr[1L], 0)
    expect_false(is.unsorted(indptr))
})

test_that("missing chunks fill with fill_value", {
    td <- withr::local_tempdir()
    x <- seq_len(20L) + 0.5
    .zf_write_array(td, "arr", x, "<f8", chunks = 6L)
    unlink(file.path(td, "arr", "1")) # drop values 7..12
    src <- GiottoDisk:::.zarr_open(td)
    on.exit(GiottoDisk:::.zarr_close(src))
    out <- GiottoDisk:::.zarr_array(src, "arr")
    expect_identical(out[7:12], rep(0, 6L))
    expect_identical(out[c(1:6, 13:20)], x[c(1:6, 13:20)])
})

test_that("chunk reader streams a multi-chunk 1D array", {
    src <- GiottoDisk:::.zarr_open(fx$paths$cell_feature_matrix)
    on.exit(GiottoDisk:::.zarr_close(src))
    rdr <- GiottoDisk:::.zarr_chunk_reader(src, "cell_features/indices",
        u4_as_integer = TRUE)
    expect_gte(rdr$n_chunks, 3L) # fixture guarantees >= 3 chunks
    got <- unlist(lapply(seq_len(rdr$n_chunks) - 1L, rdr$read_chunk))
    expect_identical(
        got,
        GiottoDisk:::.zarr_array(src, "cell_features/indices",
            u4_as_integer = TRUE)
    )
    expect_length(got, rdr$total)
})

test_that(".zarray metadata is memoized per handle", {
    td <- withr::local_tempdir()
    .zf_write_array(td, "arr", 1:10, "<i4", chunks = 10L)
    src <- GiottoDisk:::.zarr_open(td)
    on.exit(GiottoDisk:::.zarr_close(src))
    m1 <- GiottoDisk:::.zarr_meta(src, "arr")
    # rewrite the on-disk metadata; the cached parse must still be served
    writeLines("{}", file.path(td, "arr", ".zarray"))
    expect_identical(GiottoDisk:::.zarr_meta(src, "arr"), m1)
})

test_that("deflated zip entries are refused", {
    td <- withr::local_tempdir()
    # compressible payload so deflate actually kicks in
    .zf_write_array(td, "arr", rep(1L, 5000L), "<i4", chunks = 5000L)
    zf <- withr::local_tempfile(fileext = ".zip")
    zip::zip(zf, files = list.files(td, recursive = TRUE, all.files = TRUE,
        no.. = TRUE), root = td, compression_level = 9,
        include_directories = FALSE, mode = "mirror")
    src <- GiottoDisk:::.zarr_open(zf)
    on.exit(GiottoDisk:::.zarr_close(src))
    expect_error(GiottoDisk:::.zarr_array(src, "arr"), "stored")
})

test_that("detectZarrLayout returns a validated descriptor", {
    lay <- detectZarrLayout(fx$dir)
    expect_identical(lay$layout_version, "xenium-zarr-v1")
    expect_identical(lay$products$cells, fx$paths$cells)
    expect_identical(lay$products$transcripts, fx$paths$transcripts)
    expect_identical(
        lay$products$cell_feature_matrix, fx$paths$cell_feature_matrix
    )
    expect_true(file.exists(lay$gene_panel_json))
    # a directory without zarr products
    empty <- withr::local_tempdir()
    expect_identical(detectZarrLayout(empty)$layout_version, NA_character_)
})
