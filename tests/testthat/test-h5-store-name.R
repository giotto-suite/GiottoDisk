# A dense matrix is written to an h5ArrayStore whose dataset is "data", and
# read back after .store_nostate() has cleared the store's params (#157).

describe("h5ArrayStore", {
    skip_if_not_installed("HDF5Array")

    it("names its dataset data", {
        expect_identical(h5ArrayStore()@params$name, "data")
    })

    it("reads back after its params are cleared", {
        m <- matrix(as.numeric(1:12), 3, dimnames = list(letters[1:3], LETTERS[1:4]))
        s <- storeWrite(h5ArrayStore(), m)
        on.exit(unlink(s@path), add = TRUE)
        expect_equal(as.matrix(storeRead(.store_nostate(s))), m, ignore_attr = TRUE)
    })

    it("takes a dense matrix through sourceWrite() into a project", {
        td <- tempfile("h5name_")
        on.exit(unlink(td, recursive = TRUE), add = TRUE)
        src <- suppressMessages(gDirSource(td))
        m <- Matrix::Matrix(matrix(as.numeric(1:12), 3), sparse = FALSE)
        s <- sourceWrite(src, m)
        expect_s4_class(s, "h5ArrayStore")
        expect_true(sourceContains(src, s))
        expect_equal(as.matrix(storeRead(s)), as.matrix(m), ignore_attr = TRUE)
    })
})
