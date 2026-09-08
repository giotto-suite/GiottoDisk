# End-to-end: importXeniumDisk on a zarr-only directory (fixture), plus
# env-gated golden tests against a real Xenium dataset.
#
# Golden tests run only when GIOTTODISK_XENIUM_DATA points at a full
# Xenium output directory that ships BOTH the zarr archives and the 10x
# parquet files (the parquet is the independent ground truth). Known,
# accepted divergences: transcripts `cell_id` / `overlaps_nucleus` /
# `nucleus_distance` (not derivable from zarr) and `fov_name` (synthetic
# "FOV%03d"; the shipped alphanumeric codes come from instrument config).

skip_if_no_zarr_deps()

test_that("importXeniumDisk builds a gobject from a zarr-only dir", {
    skip_if_not_installed("Giotto")
    skip_if_not_installed("tilework")
    fx <- make_zarr_fixture()
    dump <- withr::local_tempdir()
    withr::local_options(giottodisk.artifact_dump = dump)
    withr::local_options(tilework.warn_sequential = FALSE)
    withr::local_options(giotto.warn_sequential = FALSE)

    rdr <- suppressWarnings(importXeniumDisk(
        xenium_dir = fx$dir,
        backend = withr::local_tempdir(),
        qv_threshold = 20
    ))
    # detection resolved every slot to the zarr archives
    expect_identical(rdr@paths$tx_path, fx$paths$transcripts)
    expect_identical(rdr@paths$cell_bound_path, fx$paths$cells)
    expect_identical(rdr@paths$nuc_bound_path, fx$paths$cells)
    expect_identical(rdr@paths$cell_meta_path, fx$paths$cells)
    expect_identical(rdr@paths$expr_path, fx$paths$cell_feature_matrix)

    g <- suppressWarnings(
        rdr$create_gobject(load_images = NULL, verbose = FALSE)
    )
    expect_s4_class(g, "giotto")
    # polygons: both sets, correct cell ids
    expect_setequal(
        GiottoClass::list_spatial_info_names(g), c("cell", "nucleus")
    )
    gpoly <- GiottoClass::getPolygonInfo(g, name = "cell",
        return_giottoPolygon = TRUE)
    expect_setequal(GiottoClass::spatIDs(gpoly), fx$truth$cell_ids)
    # transcripts: qv-filtered rna points
    tt <- fx$truth$transcripts
    keep <- tt$qv >= 20 & tt$feature_name %in% fx$truth$panel[1:4]
    gpoints <- GiottoClass::getFeatureInfo(g, feat_type = "rna",
        return_giottoPoints = TRUE)
    expect_setequal(
        GiottoClass::featIDs(gpoints), unique(tt$feature_name[keep])
    )
    # expression: rna store matches the truth matrix for the gene rows
    ex <- GiottoClass::getExpression(g, spat_unit = "cell",
        feat_type = "rna", values = "raw")
    pe <- ex[]
    expect_s4_class(pe, "parquetExprStore")
    expect_identical(pe@cell_ids, fx$truth$cell_ids)
    tri <- dplyr::collect(storeRead(pe))
    gi <- pe@gene_idx
    if (!length(gi)) gi <- seq_len(pe@n_genes)
    m <- matrix(0, nrow = 12L, ncol = length(gi))
    m[cbind(tri$row_id, match(tri$col_id, gi))] <- tri$value
    truth_rna <- t(unname(fx$truth$cfm))[, seq_along(gi)]
    expect_identical(m, truth_rna)
})

test_that("importAteraDisk inherits the zarr pathway unchanged", {
    skip_if_not_installed("Giotto")
    fx <- make_zarr_fixture()
    dump <- withr::local_tempdir()
    withr::local_options(giottodisk.artifact_dump = dump)
    rdr <- suppressWarnings(importAteraDisk(
        atera_dir = fx$dir,
        backend = withr::local_tempdir(),
        qv_threshold = 20
    ))
    expect_s4_class(rdr, "AteraDiskReader")
    expect_identical(rdr@platform, "Atera")
    expect_identical(rdr@paths$tx_path, fx$paths$transcripts)
    expect_identical(rdr@paths$expr_path, fx$paths$cell_feature_matrix)
    # one modality end to end through the inherited closures
    ex <- rdr$load_expression(verbose = FALSE)
    expect_true(any(vapply(ex, function(e) e@feat_type == "rna",
        logical(1L))))
})

# ---- golden tests against a real dataset (env-gated) ----

.xen_data_dir <- Sys.getenv("GIOTTODISK_XENIUM_DATA", "")
skip_if_no_xenium_data <- function() {
    testthat::skip_if(
        !nzchar(.xen_data_dir) || !dir.exists(.xen_data_dir),
        "GIOTTODISK_XENIUM_DATA not set"
    )
}

test_that("golden: converted cells/boundaries match 10x-shipped parquet", {
    skip_if_no_xenium_data()
    shipped_cells <- file.path(.xen_data_dir, "cells.parquet")
    shipped_cb <- file.path(.xen_data_dir, "cell_boundaries.parquet")
    skip_if(!file.exists(shipped_cells) || !file.exists(shipped_cb),
        "shipped parquet not present")
    out <- withr::local_tempdir()
    res <- xeniumZarrToParquet(.xen_data_dir, out,
        what = c("cells", "cell_boundaries"), verbose = FALSE)

    got <- as.data.frame(arrow::read_parquet(res$cells))
    ref <- as.data.frame(arrow::read_parquet(shipped_cells))
    expect_identical(nrow(got), nrow(ref))
    expect_identical(got$cell_id, ref$cell_id)
    expect_equal(got$x_centroid, ref$x_centroid, tolerance = 1e-6)
    expect_equal(got$cell_area, ref$cell_area, tolerance = 1e-6)

    gotb <- as.data.frame(arrow::read_parquet(res$cell_boundaries))
    refb <- as.data.frame(arrow::read_parquet(shipped_cb))
    expect_identical(nrow(gotb), nrow(refb))
    expect_identical(gotb$cell_id, refb$cell_id)
    expect_equal(gotb$vertex_x, refb$vertex_x, tolerance = 1e-6)
    expect_equal(gotb$vertex_y, refb$vertex_y, tolerance = 1e-6)
})

test_that("golden: converted transcripts match 10x-shipped parquet", {
    skip_if_no_xenium_data()
    shipped <- file.path(.xen_data_dir, "transcripts.parquet")
    skip_if(!file.exists(shipped), "shipped parquet not present")
    out <- withr::local_tempdir()
    res <- xeniumZarrToParquet(.xen_data_dir, out, what = "transcripts",
        verbose = FALSE)
    got <- arrow::open_dataset(res$transcripts) |>
        dplyr::select(transcript_id, feature_name, x_location,
            y_location, z_location, qv, fov_name, codeword_index) |>
        dplyr::collect() |>
        as.data.frame()
    ref <- arrow::open_dataset(shipped) |>
        dplyr::select(transcript_id, feature_name, x_location,
            y_location, z_location, qv, fov_name, codeword_index) |>
        dplyr::collect() |>
        as.data.frame()
    expect_identical(nrow(got), nrow(ref))
    o1 <- order(suppressWarnings(as.numeric(got$transcript_id)))
    o2 <- order(suppressWarnings(as.numeric(ref$transcript_id)))
    expect_identical(got$feature_name[o1], ref$feature_name[o2])
    # format v4 zarr carries the instrument fov names on .zattrs
    expect_identical(got$fov_name[o1], ref$fov_name[o2])
    expect_equal(got$x_location[o1], ref$x_location[o2], tolerance = 1e-6)
    expect_equal(got$y_location[o1], ref$y_location[o2], tolerance = 1e-6)
    expect_equal(got$qv[o1], ref$qv[o2], tolerance = 1e-6)
    expect_identical(got$codeword_index[o1], ref$codeword_index[o2])
})

test_that("golden: tenxZarrInput triplets equal tenxH5Input triplets", {
    skip_if_no_xenium_data()
    skip_if_not_installed("hdf5r")
    h5 <- file.path(.xen_data_dir, "cell_feature_matrix.h5")
    zarr <- file.path(.xen_data_dir, "cell_feature_matrix.zarr.zip")
    skip_if(!file.exists(h5) || !file.exists(zarr),
        "h5 + zarr pair not present")
    drain <- function(itr) {
        parts <- list()
        repeat {
            b <- itr$next_batch()
            if (is.null(b)) break
            if (nrow(b)) parts[[length(parts) + 1L]] <- b
        }
        out <- data.table::rbindlist(parts)
        data.table::setorder(out, row_id, col_id)
        out
    }
    tz <- tenxZarrInput(zarr)
    th <- tenxH5Input(h5)
    expect_identical(tz@cell_ids, th@cell_ids)
    expect_identical(tz@feat_ids, th@feat_ids)
    tri_z <- drain(storeRead(tz))
    tri_h <- drain(storeRead(th))
    expect_identical(nrow(tri_z), nrow(tri_h))
    expect_identical(tri_z$row_id, tri_h$row_id)
    expect_identical(tri_z$col_id, tri_h$col_id)
    expect_identical(tri_z$value, tri_h$value)
})
