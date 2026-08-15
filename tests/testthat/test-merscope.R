# MERSCOPE reader unit tests.
#
# All fixtures are synthetic, so these run without the vendor dataset. That
# also lets the true-3D branch be exercised: real MERSCOPE exports seen so far
# are replicated 2D, and there is no way to reach the 3D error with them.

# helpers ####

# Build a boundary parquet shaped like cell_boundaries.parquet.
# `mode`:
#   "replicated" - identical geometry on every z-plane (the common vendor case)
#   "3d"         - geometry genuinely differs per plane
#   "single"     - one z-plane only
.mk_bounds <- function(path, n_cells = 20L, n_z = 3L,
                       mode = c("replicated", "3d", "single")) {
    mode <- match.arg(mode)
    if (mode == "single") n_z <- 1L
    sq <- function(cx, cy, r) {
        terra::vect(sprintf(
            "POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
            cx - r, cy - r, cx + r, cy - r, cx + r, cy + r,
            cx - r, cy + r, cx - r, cy - r))
    }
    rows <- list()
    for (i in seq_len(n_cells)) {
        cx <- (i %% 5L) * 10; cy <- (i %/% 5L) * 10
        for (z in seq_len(n_z) - 1L) {
            # 3d: radius grows with z, so vertices differ per plane
            r <- if (mode == "3d") 2 + z * 0.5 else 2
            rows[[length(rows) + 1L]] <- data.frame(
                EntityID = i, ZIndex = z, ZLevel = 1.5 * (z + 1L),
                Type = "cell", stringsAsFactors = FALSE
            ) |> within(Geometry <- I(list(terra::geom(sq(cx, cy, r), wkb = TRUE)[[1L]])))
        }
    }
    df <- do.call(rbind, rows)
    arrow::write_parquet(df, path)
    path
}

# transcript parquet carrying an UNNAMED leading column, as the vendor writes it
.mk_tx <- function(path, n = 100L) {
    tb <- arrow::arrow_table(
        idx = seq_len(n),
        global_x = as.numeric(seq_len(n)),
        global_y = as.numeric(rev(seq_len(n))),
        fov = rep(c(0L, 1L), length.out = n),
        transcript_score = seq(0.95, 1.0, length.out = n),
        gene = rep(c("A", "B", "Blank-1"), length.out = n),
        cell_id = as.integer(seq_len(n))
    )
    names(tb)[1] <- ""          # <- the whole point
    arrow::write_parquet(tb, path)
    path
}


# feature classification ####

test_that("Blank features are separated from rna by keyword", {
    ids <- c("PDK4", "CCL26", "Blank-1", "Blank-84", "CD4")
    cl <- GiottoDisk:::.merscope_feat_classes(ids, "Blank")
    expect_equal(cl, c("rna", "rna", "Blank", "Blank", "rna"))
})

test_that("feature classification accepts multiple keywords", {
    ids <- c("GeneA", "Blank-1", "NegCtrl-3")
    cl <- GiottoDisk:::.merscope_feat_classes(ids, list("Blank", "NegCtrl"))
    expect_equal(cl, c("rna", "Blank", "NegCtrl"))
})


# transcript schema patch ####

test_that("unnamed index column is renamed and dropped", {
    skip_if_not_installed("arrow")
    p <- .mk_tx(tempfile(fileext = ".parquet"))

    # the raw dataset cannot be used with dplyr at all
    raw <- arrow::open_dataset(p)
    expect_true(any(!nzchar(names(raw))))
    expect_error(dplyr::collect(dplyr::filter(raw, fov == 0L)),
                 "zero-length variable name")

    # after the patch it works and the column is gone
    q <- GiottoDisk:::.merscope_tx_query(p)
    expect_false(any(!nzchar(names(q))))
    expect_false("idx" %in% names(q))
    expect_silent(dplyr::collect(dplyr::filter(q, fov == 0L)))
})

test_that("score and FOV filters are applied", {
    p <- .mk_tx(tempfile(fileext = ".parquet"), n = 100L)
    all_n <- nrow(dplyr::collect(GiottoDisk:::.merscope_tx_query(p)))
    expect_equal(all_n, 100L)

    fov0 <- dplyr::collect(GiottoDisk:::.merscope_tx_query(p, FOVs = 0L))
    expect_true(all(fov0$fov == 0L))

    hi <- dplyr::collect(
        GiottoDisk:::.merscope_tx_query(p, score_threshold = 0.99))
    expect_true(all(hi$transcript_score >= 0.99))
    expect_lt(nrow(hi), all_n)
})

test_that("cell_id is cast to character so it joins against poly_ID", {
    p <- .mk_tx(tempfile(fileext = ".parquet"))
    q <- dplyr::collect(GiottoDisk:::.merscope_tx_query(p))
    expect_type(q$cell_id, "character")
})


# z-plane architecture ####

test_that("replicated 2D is detected and redundant planes dropped", {
    p <- .mk_bounds(tempfile(fileext = ".parquet"), n_z = 4L, mode = "replicated")
    arch <- GiottoDisk:::.merscope_zplane_architecture(p, verbose = FALSE)
    expect_equal(arch$kind, "replicated_2d")
    expect_equal(arch$use_z, 0L)
    expect_length(arch$z_indices, 4L)
})

test_that("a single z-plane is reported as plain 2D", {
    p <- .mk_bounds(tempfile(fileext = ".parquet"), mode = "single")
    arch <- GiottoDisk:::.merscope_zplane_architecture(p, verbose = FALSE)
    expect_equal(arch$kind, "2d")
})

test_that("true 3D raises an explicit error naming aggregateStacks", {
    p <- .mk_bounds(tempfile(fileext = ".parquet"), n_z = 3L, mode = "3d")
    expect_error(
        GiottoDisk:::.merscope_zplane_architecture(p, verbose = FALSE),
        "true 3D segmentation detected"
    )
    expect_error(
        GiottoDisk:::.merscope_zplane_architecture(p, verbose = FALSE),
        "aggregateStacks"
    )
})

test_that("poly_z_indices selects a subset, and a single index escapes the 3D error", {
    p <- .mk_bounds(tempfile(fileext = ".parquet"), n_z = 3L, mode = "3d")
    # picking one plane explicitly is the documented escape hatch
    arch <- GiottoDisk:::.merscope_zplane_architecture(
        p, poly_z_indices = 1L, verbose = FALSE)
    expect_equal(arch$kind, "2d")
    expect_equal(arch$use_z, 1L)
})

test_that("poly_z_indices rejects indices that are not present", {
    p <- .mk_bounds(tempfile(fileext = ".parquet"), n_z = 2L)
    expect_error(
        GiottoDisk:::.merscope_zplane_architecture(p, poly_z_indices = 99L),
        "none of `poly_z_indices`"
    )
})


# path detection ####

test_that("region files are found and experiment.json is read from the parent", {
    root <- file.path(tempfile(), "EXP_1")
    region <- file.path(root, "region_R1")
    dir.create(region, recursive = TRUE)
    for (f in c("detected_transcripts.parquet", "cell_boundaries.parquet",
                "cell_by_gene.csv", "cell_metadata.csv"))
        file.create(file.path(region, f))
    file.create(file.path(root, "experiment.json"))   # parent level

    p <- GiottoDisk:::.merscope_detect_paths(region)
    expect_match(p$tx_path, "detected_transcripts\\.parquet$")
    expect_match(p$cell_bound_path, "cell_boundaries\\.parquet$")
    expect_match(p$expr_path, "cell_by_gene\\.csv$")
    expect_match(p$experiment_path, "experiment\\.json$")
    # absent modalities come back as NA, not NULL
    expect_true(is.na(p$image_dir))
})

test_that("a directory with no MERSCOPE files is rejected", {
    d <- tempfile(); dir.create(d)
    expect_error(GiottoDisk:::.merscope_detect_paths(d),
                 "Is this a MERSCOPE region directory")
})


# reader construction ####

test_that("importMerscopeDisk requires a backend", {
    expect_error(importMerscopeDisk("somewhere"), "`backend` is required")
})

test_that("hdf5 boundaries are refused with an explicit message", {
    skip_if_not_installed("terra")
    d <- tempfile(); dir.create(d)
    expect_error(
        GiottoDisk:::.merscope_poly_disk(
            path = tempfile(), gsource = structure(list(), class = "gsource"),
            polygon_format = "hdf5"),
        "not implemented yet"
    )
})
