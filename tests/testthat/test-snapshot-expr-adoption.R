# snapshotSave must adopt this package's own parquet expression stores, not just
# BPCells IterableMatrix. Before the fix, `.ss_gdsrc_register_external_expr()`
# skipped anything that was not an IterableMatrix, so a parquetExprStore written
# to the default tempdir() dump was never moved into the vault -- a reloaded
# snapshot then pointed at a directory the OS had deleted.
#
# The subtlety: lazy ops mean `raw` and `normalized` are two handles onto ONE
# on-disk path, and adoption MOVES files. The second handle must be remapped via
# the adoption session map rather than re-adopted or rejected.

.mk_pe <- function(nf = 8L, nc = 12L, seed = 1L) {
    set.seed(seed)
    m <- Matrix::Matrix(rpois(nf * nc, 3), nrow = nf, sparse = TRUE)
    dimnames(m) <- list(paste0("g", seq_len(nf)), paste0("c", seq_len(nc)))
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
}

.mk_src <- function(tag) {
    td <- tempfile(tag); dir.create(td, recursive = TRUE)
    list(src = gDirSource(td), td = td)
}


test_that("a parquetExprStore is a fileStore and not an IterableMatrix", {
    pe <- .mk_pe()
    # this is precisely why the old class test skipped it
    expect_true(inherits(pe, "fileStore"))
    expect_false(inherits(pe, "IterableMatrix"))
})


test_that("sourceAdopt moves an external parquet expression store into the vault", {
    x <- .mk_src("adoptExpr_"); on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    GiottoDisk:::.adopt_session_reset()

    pe <- .mk_pe()
    before <- as.data.frame(storeRead(pe, output = "tibble"))
    old_path <- pe@path
    expect_false(startsWith(normalizePath(old_path), normalizePath(x$td)))

    pe2 <- sourceAdopt(x$src, pe)

    expect_true(startsWith(normalizePath(pe2@path), normalizePath(x$td)))
    expect_true(storeExists(pe2))
    expect_false(file.exists(old_path))          # moved, not copied
    after <- as.data.frame(storeRead(pe2, output = "tibble"))
    expect_equal(after[order(after$row_id, after$col_id), ],
                 before[order(before$row_id, before$col_id), ],
                 ignore_attr = TRUE)              # values survive the move
})


test_that("two handles on one path: the second is remapped, not rejected", {
    x <- .mk_src("adoptShared_"); on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    GiottoDisk:::.adopt_session_reset()

    raw <- .mk_pe()
    # what normalizeGiotto() produces: same @path, one extra lazy op
    norm <- processData(raw, new("libraryNormParam", param = list(scalefactor = 1e4)))
    expect_identical(raw@path, norm@path)
    expect_gt(length(norm@ops), length(raw@ops))

    raw2 <- sourceAdopt(x$src, raw)
    # the files have moved out from under `norm`; adoption must remap it
    norm2 <- sourceAdopt(x$src, norm)

    expect_identical(raw2@path, norm2@path)
    expect_true(storeExists(norm2))
    expect_gt(length(norm2@ops), 0L)             # ops preserved through remap
    expect_true(startsWith(normalizePath(norm2@path), normalizePath(x$td)))
})


test_that("an adopted store reads as vault-resident, so snapshotSave skips it", {
    x <- .mk_src("adoptIdem_"); on.exit(unlink(x$td, recursive = TRUE), add = TRUE)
    GiottoDisk:::.adopt_session_reset()

    pe <- .mk_pe()
    expect_false(sourceContains(x$src, pe))      # external before adoption
    pe2 <- sourceAdopt(x$src, pe)
    # this is the guard snapshotSave uses to avoid re-adopting on a second save
    expect_true(sourceContains(x$src, pe2))
})
