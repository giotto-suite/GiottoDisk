# Tests for snapshotSave/snapshotLoad on giottoMulti.
#
# A multi save writes one .rds at the multi vault and one per-child .rds in
# each child's vault, named "<multi_name>_<child_name>". The multi .rds
# carries the entire reconstituted topology — load is one readRDS().
#
# These tests require a build of GiottoClass that carries @source on
# giottoMulti (gmulti branch).

skip_if_not_installed("GiottoClass")
skip_if_not_installed("GiottoData")

skip_if(
    !"source" %in% methods::slotNames("giottoMulti"),
    "GiottoClass lacks @source on giottoMulti (pre-gmulti build)"
)

.mk_test_multi <- function(td) {
    g1 <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    g2 <- GiottoData::loadGiottoMini("visium", verbose = FALSE)

    src_multi <- gDirSource(file.path(td, "multi_proj"))
    src1 <- gDirSource(file.path(td, "sample1_proj"))
    src2 <- gDirSource(file.path(td, "sample2_proj"))
    g1@source <- src1
    g2@source <- src2

    GiottoClass::createGiottoMulti(
        list(s1 = g1, s2 = g2),
        source = src_multi
    )
}


test_that("snapshotSave(gDirSource, giottoMulti) writes multi + per-child .rds", {
    td <- tempfile("multi_snap_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    mg <- .mk_test_multi(td)
    snapshotSave(mg@source, mg, verbose = FALSE)

    # Multi vault gets <name>.rds
    multi_gsd <- file.path(mg@source@path, "giottosave")
    expect_true(dir.exists(multi_gsd))
    multi_rds <- list.files(multi_gsd, pattern = "\\.rds$")
    expect_length(multi_rds, 1L)

    # Each child vault gets <multi_name>_<child_name>.rds
    multi_name <- tools::file_path_sans_ext(multi_rds)
    for (child_name in names(mg@objects)) {
        child_gsd <- file.path(mg@objects[[child_name]]@source@path, "giottosave")
        expect_true(dir.exists(child_gsd))
        child_rds <- list.files(child_gsd, pattern = "\\.rds$")
        expect_length(child_rds, 1L)
        expect_identical(
            tools::file_path_sans_ext(child_rds),
            paste0(multi_name, "_", child_name)
        )
    }
})


test_that("snapshotLoad reconstitutes the giottoMulti with children intact", {
    td <- tempfile("multi_load_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    mg <- .mk_test_multi(td)
    snapshotSave(mg@source, mg, verbose = FALSE)

    mg2 <- snapshotLoad(mg@source)
    expect_s4_class(mg2, "giottoMulti")
    expect_identical(names(mg2@objects), c("s1", "s2"))

    # children retain their own sources (compare basenames — macOS
    # /var/folders is a symlink to /private/var/folders so absolute
    # paths differ post-normalize).
    expect_identical(basename(mg2@objects$s1@source@path), "sample1_proj")
    expect_identical(basename(mg2@objects$s2@source@path), "sample2_proj")
    expect_identical(basename(mg2@source@path), "multi_proj")
})


test_that("snapshotSave on a multi with in-memory child writes only multi .rds", {
    td <- tempfile("multi_inmem_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g1 <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    src_multi <- gDirSource(file.path(td, "multi_proj"))
    # No child source — child is in-memory only
    mg <- GiottoClass::createGiottoMulti(
        list(only = g1),
        source = src_multi
    )

    snapshotSave(mg@source, mg, verbose = FALSE)

    multi_gsd <- file.path(mg@source@path, "giottosave")
    expect_true(dir.exists(multi_gsd))
    expect_length(list.files(multi_gsd, pattern = "\\.rds$"), 1L)

    # No per-child vault created since child has no @source
    expect_false(dir.exists(file.path(td, "only_proj")))
})


test_that("snapshotDelete on a giottoMulti cascades to per-child snapshots", {
    td <- tempfile("multi_delete_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    mg <- .mk_test_multi(td)
    snapshotSave(mg@source, mg, name = "snap_v1", verbose = FALSE)

    # Pre-condition: multi + per-child snapshots all exist on disk
    multi_rds <- file.path(mg@source@path, "giottosave", "snap_v1.rds")
    s1_rds <- file.path(mg@objects$s1@source@path, "giottosave", "snap_v1_s1.rds")
    s2_rds <- file.path(mg@objects$s2@source@path, "giottosave", "snap_v1_s2.rds")
    expect_true(all(file.exists(c(multi_rds, s1_rds, s2_rds))))

    # Delete the multi → per-child snapshots are removed too
    snapshotDelete(mg@source, "snap_v1")
    expect_false(file.exists(multi_rds))
    expect_false(file.exists(s1_rds))
    expect_false(file.exists(s2_rds))
})


test_that("snapshotDelete cascade is best-effort: parent removed even if a child fails", {
    td <- tempfile("multi_delete_be_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    mg <- .mk_test_multi(td)
    snapshotSave(mg@source, mg, name = "snap_v1", verbose = FALSE)

    # Pre-emptively delete one child snapshot to force a cascade failure
    s1_rds <- file.path(mg@objects$s1@source@path, "giottosave", "snap_v1_s1.rds")
    unlink(s1_rds)
    multi_rds <- file.path(mg@source@path, "giottosave", "snap_v1.rds")
    s2_rds <- file.path(mg@objects$s2@source@path, "giottosave", "snap_v1_s2.rds")
    expect_true(file.exists(multi_rds))
    expect_true(file.exists(s2_rds))

    # Best-effort: warn on missing s1 cascade, still delete multi + s2
    expect_warning(snapshotDelete(mg@source, "snap_v1"),
                   "child snapshot 'snap_v1_s1'.*failed")
    expect_false(file.exists(multi_rds))
    expect_false(file.exists(s2_rds))
})


test_that("multi snapshotSave overwrite gate fires on name collision", {
    td <- tempfile("multi_overwrite_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    mg <- .mk_test_multi(td)
    snapshotSave(mg@source, mg, name = "snap_v1", verbose = FALSE)

    expect_error(
        snapshotSave(mg@source, mg, name = "snap_v1", verbose = FALSE),
        "already exists"
    )

    expect_no_error(
        snapshotSave(mg@source, mg, name = "snap_v1",
            overwrite = TRUE, verbose = FALSE)
    )
})
