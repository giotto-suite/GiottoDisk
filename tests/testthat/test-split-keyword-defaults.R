# Atera WTA panels name their genomic control probes `Intergenic_Region_*`
# instead of `GenomicControl*`, and ship `DeprecatedCodeword_*` features. Every
# disk reader's defaults have to split both out, or they load as `rna` and
# aggregation adds them to the gene expression.

skip_if_no_zarr_deps()

.default <- function(f, arg) eval(formals(f)[[arg]])

test_that("Atera and Xenium disk readers share the same defaults", {
    skip_if_not_installed("Giotto")
    withr::local_options(giottodisk.artifact_dump = withr::local_tempdir())

    atera <- importAteraDisk(test_path("fixtures", "atera_v6"),
        backend = withr::local_tempdir())
    xenium <- suppressWarnings(importXeniumDisk(make_zarr_fixture()$dir,
        backend = withr::local_tempdir()))

    ref <- GiottoDisk:::.xenium_transcript_disk
    expect_true("Intergenic_Region" %in% unlist(.default(ref, "split_keyword")))
    expect_true("DeprecatedCodeword" %in% .default(ref, "feat_type"))
    for (f in list(atera$load_transcripts, atera$create_gobject,
        xenium$load_transcripts, xenium$create_gobject)) {
        for (arg in c("feat_type", "split_keyword")) {
            expect_identical(.default(f, arg), .default(ref, arg))
        }
    }
})

test_that("the default split keeps control probes out of rna", {
    skip_if_not_installed("Giotto")
    ref <- GiottoDisk:::.xenium_transcript_disk
    tx <- data.table::data.table(
        feat_ID = c("GENE_A", "NegControlProbe_00002", "GenomicControl_0001",
            "Intergenic_Region_100", "DeprecatedCodeword_11440"),
        x = 1:5, y = 1:5
    )
    gp <- GiottoClass::createGiottoPoints(tx,
        feat_type = .default(ref, "feat_type"),
        split_keyword = .default(ref, "split_keyword"),
        verbose = FALSE
    )
    names(gp) <- vapply(gp, GiottoClass::objName, character(1L))

    expect_identical(GiottoClass::featIDs(gp$rna), "GENE_A")
    expect_setequal(GiottoClass::featIDs(gp$GenomicControl),
        c("GenomicControl_0001", "Intergenic_Region_100"))
    expect_identical(GiottoClass::featIDs(gp$DeprecatedCodeword),
        "DeprecatedCodeword_11440")
})
