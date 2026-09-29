# Technology-agnostic reader pieces (R/reader-shared.R), on synthetic data.
# Ground truth is the input matrix / geometry itself, not a second run of the
# same path.

# 3 genes x 6 bins
.grid_fixture <- function() {
    mat <- Matrix::sparseMatrix(
        i = c(1L, 2L, 1L, 3L, 2L, 1L, 3L, 2L),
        j = c(1L, 1L, 2L, 3L, 4L, 5L, 5L, 6L),
        x = c(4, 1, 2, 7, 3, 5, 1, 6),
        dims = c(3L, 6L),
        dimnames = list(c("g1", "g2", "g3"), sprintf("b%d", 1:6))
    )
    coords <- data.table::data.table(
        cell_ID = colnames(mat),
        sdimx = rep(c(0.5, 1.5, 2.5), 2L),
        sdimy = rep(c(0.5, 1.5), each = 3L)
    )
    list(mat = mat, coords = coords)
}

test_that(".add_parent_units keeps only one-to-one parents", {
    fx <- .grid_fixture()
    g <- GiottoClass::createGiottoObject(expression = fx$mat, verbose = FALSE)

    # two rows per unit (as a finer layer would give): `blk` is a function of
    # the unit, `cellx` is not for b2, the flag is aggregated with any()
    mapping <- data.table::data.table(
        cell_ID = rep(colnames(fx$mat), each = 2L),
        blk     = rep(c("L", "R", "R", "L", "R", "R"), each = 2L),
        cellx   = c("c1", "c1", "c1", "c2", NA, NA, "c3", "c3", NA, NA, NA, NA),
        in_x    = c(TRUE, FALSE, rep(FALSE, 10L))
    )
    g2 <- .add_parent_units(g, spat_unit = "cell", mapping = mapping)
    cm <- GiottoClass::pDataDT(g2)
    expect_false("cellx" %in% names(cm))
    expect_identical(cm[match(colnames(fx$mat), cell_ID)]$blk,
                     c("L", "R", "R", "L", "R", "R"))
    expect_identical(cm[match(colnames(fx$mat), cell_ID)]$in_x,
                     c(TRUE, rep(FALSE, 5L)))

    # one row per unit: taken as-is
    one <- mapping[seq(1L, .N, by = 2L)]
    cm1 <- GiottoClass::pDataDT(.add_parent_units(g, "cell", one))
    expect_identical(cm1[match(one$cell_ID, cell_ID)]$cellx, one$cellx)
})
