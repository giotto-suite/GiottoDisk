# fingerprint_data() on stores: identity from the handle, never the data

.fp_pe <- function(seed = 1L, prefix = "c") {
    storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
        .rand_counts(seed = seed, cell_prefix = prefix))
}

.fp_geom <- function() {
    sv <- terra::vect(
        data.frame(x = 1:5, y = 1:5, cell_ID = paste0("cell_", 1:5)),
        geom = c("x", "y"), crs = ""
    )
    storeWrite(parquetGeomStore(), sv)
}

.fp_edges <- function() {
    dt <- data.table::data.table(
        from = c("a", "a", "b", "c"), to = c("b", "c", "c", "d"),
        weight = c(1, 2, 3, 4)
    )
    storeWrite(storeCreate(type = "parquetEdgeStore"), dt,
        type = "sNN", directed = FALSE)
}

.fp_src <- function() {
    gdir <- file.path(tempdir(), paste0("fp_", basename(tempfile())))
    gDirSource(path = gdir)
}

describe("fingerprint_data on stores", {
    it("returns the same hash for the same handle", {
        pe <- .fp_pe()
        expect_type(fingerprint_data(pe), "character")
        expect_identical(fingerprint_data(pe), fingerprint_data(pe))
    })

    it("tells two artifacts apart even with the same content", {
        expect_false(identical(fingerprint_data(.fp_pe(1L)),
            fingerprint_data(.fp_pe(1L))))
    })

    it("changes when the view changes", {
        pe <- .fp_pe()
        fp <- fingerprint_data(pe)
        expect_false(identical(fp, fingerprint_data(pe[1:3, ])))
        expect_false(identical(fp, fingerprint_data(pe[, 1:3])))
        expect_false(identical(fp, fingerprint_data(pe > 2)))

        g <- .fp_geom()
        fg <- fingerprint_data(g)
        expect_false(identical(fg,
            fingerprint_data(subset(g, subset = cell_ID == "cell_1"))))
        expect_false(identical(fg,
            fingerprint_data(terra::crop(g, terra::ext(0, 2, 0, 2)))))
    })

    it("does not change when cached marginals are filled in", {
        pe <- .fp_pe()
        pe2 <- pe
        pe2@stats <- list(cell_sum = stats::setNames(1:30, NULL))
        expect_identical(fingerprint_data(pe), fingerprint_data(pe2))
    })

    it("ignores the CRS", {
        g <- .fp_geom()
        g2 <- g
        g2@params$crs <- "EPSG:4326"
        expect_identical(fingerprint_data(g), fingerprint_data(g2))
    })

    it("survives adoption into a project vault", {
        src <- .fp_src()
        on.exit(unlink(src@path, recursive = TRUE), add = TRUE)

        pe <- .fp_pe()
        expect_identical(fingerprint_data(sourceAdopt(src, pe)),
            fingerprint_data(pe))

        g <- .fp_geom()
        expect_identical(fingerprint_data(sourceAdopt(src, g)),
            fingerprint_data(g))

        e <- .fp_edges()
        expect_identical(fingerprint_data(sourceAdopt(src, e)),
            fingerprint_data(e))
    })

    it("hashes a queued arrow id table by its values", {
        g <- .fp_geom()
        keep <- function() {
            list(arrow = arrow::arrow_table(cell_ID = c("cell_1", "cell_2")))
        }
        a <- GiottoDisk:::.queue_keep(g, keep())
        b <- GiottoDisk:::.queue_keep(g, keep())
        expect_identical(fingerprint_data(a), fingerprint_data(b))
    })

    it("fingerprints a union store through its substores", {
        a <- .fp_pe(1L, "a")
        u <- unionParquetExprStore(list(a, .fp_pe(2L, "b")))
        u2 <- unionParquetExprStore(list(a, .fp_pe(2L, "b")))
        expect_type(fingerprint_data(u), "character")
        expect_false(identical(fingerprint_data(u), fingerprint_data(u2)))
        expect_false(identical(fingerprint_data(u), fingerprint_data(u[1:3, ])))
    })

    it("counts a function in an op by its code", {
        pe <- .fp_pe()
        op <- function(f) {
            pe@ops <- c(pe@ops, list(list(type = "custom", fun = f)))
            pe
        }
        expect_identical(fingerprint_data(op(function(x) x)),
            fingerprint_data(op(function(x) x)))
        expect_false(identical(fingerprint_data(op(function(x) x)),
            fingerprint_data(op(function(x) x + 1))))
    })

    it("errors on a handle value it has no canonical form for", {
        pe <- .fp_pe()
        pe@ops <- c(pe@ops, list(list(type = "custom", env = new.env())))
        expect_error(fingerprint_data(pe), "no canonical form for environment")
    })
})
