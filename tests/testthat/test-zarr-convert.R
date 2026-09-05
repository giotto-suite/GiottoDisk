# Standalone converter (xeniumZarrToParquet) and the fingerprint-keyed
# conversion cache used by the disk readers (.zarr_ensure_parquet).

skip_if_no_zarr_deps()

fx <- make_zarr_fixture()

test_that("xeniumZarrToParquet converts all products with a manifest", {
    out <- withr::local_tempdir()
    res <- xeniumZarrToParquet(fx$dir, out, verbose = FALSE)
    expect_true(file.exists(res$cells))
    expect_true(file.exists(res$cell_boundaries))
    expect_true(file.exists(res$nucleus_boundaries))
    expect_true(file.exists(res$transcripts) || dir.exists(res$transcripts))
    expect_true(file.exists(res$cell_feature_matrix))
    man <- jsonlite::fromJSON(file.path(out, "conversion_manifest.json"))
    expect_identical(man$layout_version, "xenium-zarr-v1")
    expect_true("transcripts" %in% names(man$products))

    # cell_feature_matrix drains the tenxZarrInput iterator: triplet
    # layout, sorted by row_id, matching the dense truth matrix
    tri <- as.data.frame(arrow::read_parquet(res$cell_feature_matrix))
    expect_identical(colnames(tri), c("row_id", "col_id", "value"))
    expect_false(is.unsorted(tri$row_id))
    m <- matrix(0, nrow = 12L, ncol = 6L)
    m[cbind(tri$row_id, tri$col_id)] <- tri$value
    expect_identical(m, t(unname(fx$truth$cfm)))
})

test_that("existing outputs are protected unless overwrite = TRUE", {
    out <- withr::local_tempdir()
    xeniumZarrToParquet(fx$dir, out, what = "cells", verbose = FALSE)
    expect_error(
        xeniumZarrToParquet(fx$dir, out, what = "cells", verbose = FALSE),
        "overwrite"
    )
    expect_silent(
        xeniumZarrToParquet(fx$dir, out, what = "cells",
            overwrite = TRUE, verbose = FALSE)
    )
})

test_that(".zarr_ensure_parquet converts once and reuses the cache", {
    dump <- withr::local_tempdir()
    withr::local_options(giottodisk.artifact_dump = dump)
    p1 <- GiottoDisk:::.zarr_ensure_parquet(fx$paths$cells, "cells",
        verbose = FALSE)
    mt1 <- file.info(p1)$mtime
    Sys.sleep(1)
    # cache hit: same file, untouched
    p2 <- GiottoDisk:::.zarr_ensure_parquet(fx$paths$cells, "cells",
        verbose = FALSE)
    expect_identical(p1, p2)
    expect_identical(mt1, file.info(p2)$mtime)
    # sibling products of the same archive came from the one pass
    p3 <- GiottoDisk:::.zarr_ensure_parquet(fx$paths$cells,
        "nucleus_boundaries", verbose = FALSE)
    expect_identical(dirname(p3), dirname(p1))
    expect_identical(mt1, file.info(p1)$mtime) # still not reconverted
})

test_that("touching the source invalidates the fingerprint", {
    dump <- withr::local_tempdir()
    withr::local_options(giottodisk.artifact_dump = dump)
    p1 <- GiottoDisk:::.zarr_ensure_parquet(fx$paths$cells, "cells",
        verbose = FALSE)
    Sys.setFileTime(fx$paths$cells, Sys.time() + 10)
    p2 <- GiottoDisk:::.zarr_ensure_parquet(fx$paths$cells, "cells",
        verbose = FALSE)
    expect_false(identical(dirname(p1), dirname(p2)))
})
