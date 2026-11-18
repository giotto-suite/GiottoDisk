# 10x feature names on disk follow Giotto's in-memory readers: every copy of a
# repeated name gets `--N`, counted within each feature type after the split.

skip_if_not_installed("hdf5r")

# A gene repeated within its type (G1) and a gene that shares its name with an
# antibody (CD4), as an h5 and an mtx folder (helper-visiumhd-fixture.R).
.write_tenx_types <- function(dir) {
    m <- Matrix::sparseMatrix(i = c(1, 2, 3, 4, 5, 1), j = c(1, 1, 2, 3, 3, 2),
                              x = c(5, 2, 7, 1, 3, 4), dims = c(5, 3),
                              dimnames = list(c("G1", "G1", "CD4", "CD4", "G2"),
                                              paste0("c", 1:3)))
    .vhd_fx_write_matrix(m, dir, "raw_feature_bc_matrix", types = c(
        rep("Gene Expression", 3), "Antibody Capture", "Gene Expression"))
}

test_that(".disambiguate_feat_ids counts within types when given", {
    ids <- c("G1", "G1", "CD4", "CD4", "G2")
    types <- c("g", "g", "g", "a", "g")
    expect_identical(.disambiguate_feat_ids(ids),
                     c("G1--1", "G1--2", "CD4--1", "CD4--2", "G2"))
    expect_identical(.disambiguate_feat_ids(ids, types),
                     c("G1--1", "G1--2", "CD4", "CD4", "G2"))
})

test_that("split 10x stores are named within each type, from h5 and mtx", {
    d <- withr::local_tempdir()
    .write_tenx_types(d)
    for (p in c("raw_feature_bc_matrix", "raw_feature_bc_matrix.h5")) {
        st <- .tenx_expression_disk(file.path(d, p),
            gsource = gDirSource(withr::local_tempdir()), output = "store")
        expect_identical(rownames(st[["rna"]]), c("G1--1", "G1--2", "CD4", "G2"))
        expect_identical(rownames(st[["Antibody_Capture"]]), "CD4")
    }
})

test_that("disk names match Giotto's in-memory readers", {
    skip_if_not_installed("Giotto")
    skip_if_not(exists(".tenx_unique_feat_ids", asNamespace("Giotto")))
    d <- withr::local_tempdir()
    .write_tenx_types(d)
    mem <- Giotto::get10Xmatrix_h5(file.path(d, "raw_feature_bc_matrix.h5"))
    st <- .tenx_expression_disk(file.path(d, "raw_feature_bc_matrix.h5"),
        gsource = gDirSource(withr::local_tempdir()), output = "store")
    expect_identical(rownames(st[["rna"]]), rownames(mem[["Gene Expression"]]))
    expect_identical(rownames(st[["Antibody_Capture"]]),
                     rownames(mem[["Antibody Capture"]]))
})
