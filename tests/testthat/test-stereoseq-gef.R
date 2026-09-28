# Tests for the Stereo-seq GEF inputs: cellbinGefInput / binGefInput and
# their ingest into a parquetExprStore.
#
# The fixtures below write real HDF5 compound datasets in the layouts the
# readers address, so the group and field names are under test too -- that is
# the part a mocked iterator would not catch. Ground truth is the matrix the
# fixture was built from, never a second run of the same read path.

skip_if_not_installed("rhdf5")

# helper: gene-major triplet records for a genes x cells matrix. Mirrors how
# a GEF stores expression -- all records for gene 1, then gene 2, and so on.
.gef_records <- function(mat) {
    sm <- Matrix::summary(methods::as(mat, "CsparseMatrix"))
    sm <- sm[order(sm$i, sm$j), , drop = FALSE]
    data.frame(
        gene_row = as.integer(sm$i),
        cell_col = as.integer(sm$j),
        count    = as.integer(sm$x)
    )
}

# helper: write a cellbin .gef. `gene_names` may repeat or name a gene with no
# records, both of which occur in real files. Like SAW, it writes the matrix
# twice: gene-major `geneExp` and cell-major `cellExp`, the latter indexed by
# the cell table's `offset` / `geneCount`. `cell_exp = FALSE` leaves both out.
.write_cellbin_gef <- function(mat, path, gene_names = rownames(mat),
                               cell_exp = TRUE) {
    rec <- .gef_records(mat)
    n_genes <- nrow(mat)
    n_cells <- ncol(mat)
    cell_id <- seq_len(n_cells)          # on-disk ids; cell_IDs are "cell_<id>"

    per_gene <- tabulate(rec$gene_row, nbins = n_genes)
    per_cell <- tabulate(rec$cell_col, nbins = n_cells)

    cell_tab <- data.frame(id = as.integer(cell_id),
                           x  = as.integer(seq_len(n_cells) * 10L),
                           y  = as.integer(seq_len(n_cells) * 20L))
    if (cell_exp) {
        cell_tab$offset    <- as.integer(c(0L, cumsum(per_cell)[-n_cells]))
        cell_tab$geneCount <- as.integer(per_cell)
    }

    rhdf5::h5createFile(path)
    rhdf5::h5createGroup(path, "cellBin")
    rhdf5::h5write(cell_tab, path, "cellBin/cell")
    if (cell_exp) {
        cm <- rec[order(rec$cell_col, rec$gene_row), , drop = FALSE]
        rhdf5::h5write(
            data.frame(geneID = as.integer(cm$gene_row - 1L),   # 0-based row
                       count  = as.integer(cm$count)),
            path, "cellBin/cellExp"
        )
    }
    rhdf5::h5write(
        data.frame(geneName  = as.character(gene_names),
                   geneID    = paste0("ENSG", sprintf("%05d", seq_len(n_genes))),
                   cellCount = as.integer(per_gene)),
        path, "cellBin/gene"
    )
    rhdf5::h5write(
        data.frame(cellID = as.integer(cell_id[rec$cell_col]),
                   count  = as.integer(rec$count)),
        path, "cellBin/geneExp"
    )
    rhdf5::h5closeAll()
    path
}

# helper: write a bin .gef under `geneExp/<bin_key>/`. `coords` supplies the
# (x, y) for each cell column, so the test controls first-appearance order.
.write_bin_gef <- function(mat, path, bin_key = "bin100", coords = NULL,
                           y_attrs = TRUE) {
    rec <- .gef_records(mat)
    n_genes <- nrow(mat)
    n_cells <- ncol(mat)
    if (is.null(coords)) {
        coords <- data.frame(x = seq_len(n_cells) * 5L,
                             y = rev(seq_len(n_cells)) * 7L)
    }
    per_gene <- tabulate(rec$gene_row, nbins = n_genes)

    rhdf5::h5createFile(path)
    rhdf5::h5createGroup(path, "geneExp")
    rhdf5::h5createGroup(path, file.path("geneExp", bin_key))
    rhdf5::h5write(
        data.frame(geneName = as.character(rownames(mat)),
                   geneID   = paste0("ENSG", sprintf("%05d", seq_len(n_genes))),
                   count    = as.integer(per_gene)),
        path, paste0("geneExp/", bin_key, "/gene")
    )
    ex_name <- paste0("geneExp/", bin_key, "/expression")
    rhdf5::h5write(
        data.frame(x     = as.integer(coords$x[rec$cell_col]),
                   y     = as.integer(coords$y[rec$cell_col]),
                   count = as.integer(rec$count)),
        path, ex_name
    )
    # SAW records the coordinate range on the expression dataset
    fid <- rhdf5::H5Fopen(path)
    did <- rhdf5::H5Dopen(fid, ex_name)
    for (nm in if (y_attrs) c("minY", "maxY") else character(0L)) {
        v <- coords$y
        rhdf5::h5writeAttribute(as.integer(if (startsWith(nm, "min")) min(v) else max(v)),
                                did, nm)
    }
    rhdf5::H5Dclose(did); rhdf5::H5Fclose(fid)
    rhdf5::h5closeAll()
    list(path = path, records = rec, coords = coords)
}

# helper: per part file, the row_id range it holds, in part order
.part_ranges <- function(out) {
    parts <- list.files(out, pattern = "[.]parquet$", recursive = TRUE, full.names = TRUE)
    parts <- parts[order(as.integer(sub(".*part-([0-9]+)[.]parquet$", "\\1", parts)))]
    do.call(rbind, lapply(parts, function(f) {
        r <- range(arrow::read_parquet(f, col_select = "row_id")$row_id)
        data.frame(lo = r[1L], hi = r[2L])
    }))
}

# cell-major: each part holds one ascending, non-overlapping range of cells
.expect_cell_major <- function(out) {
    pr <- .part_ranges(out)
    if (nrow(pr) > 1L) expect_true(all(pr$lo[-1L] > pr$hi[-nrow(pr)]))
    invisible(pr)
}

# helper: pull a store back into a genes x cells sparse matrix
.pe_as_matrix <- function(pe) {
    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    m <- Matrix::sparseMatrix(
        i = df$col_id, j = df$row_id, x = df$value,
        dims = c(nrow(pe), ncol(pe))
    )
    dimnames(m) <- list(pe@feat_ids, pe@cell_ids)
    m
}


# ---- cellbin ---------------------------------------------------------------

test_that("cellbinGefInput + storeWrite round-trips a cellbin gef losslessly", {
    # Gene names deliberately out of alphabetical order: the input sorts
    # feat_ids, so a correct col_id mapping has to survive the reorder.
    mat <- Matrix::sparseMatrix(
        i = c(1L, 1L, 2L, 3L, 3L, 3L),
        j = c(1L, 3L, 2L, 1L, 2L, 4L),
        x = c(5, 7, 2, 1, 9, 4),
        dims = c(3L, 4L)
    )
    rownames(mat) <- c("gzmb", "actb", "cd8a")
    colnames(mat) <- paste0("cell_", 1:4)

    gef <- .write_cellbin_gef(mat, file.path(tempdir(), "cellbin_rt.gef"))
    out <- file.path(tempdir(), "cellbin_rt_out")
    on.exit(unlink(c(gef, out), recursive = TRUE), add = TRUE)

    inp <- cellbinGefInput(gef)
    pe  <- storeWrite(parquetExprStore(path = out), inp)

    expect_s4_class(pe, "parquetExprStore")
    expect_equal(pe@cell_ids, paste0("cell_", 1:4))
    expect_equal(pe@feat_ids, sort(rownames(mat)))
    expect_equal(pe@n_cells, 4)
    expect_equal(pe@n_genes, 3)

    expected <- mat[sort(rownames(mat)), , drop = FALSE]
    expect_equal(as.matrix(.pe_as_matrix(pe)), as.matrix(expected))
})

test_that("cellbinGefInput sums duplicate gene names and drops empty genes", {
    # Rows 1 and 3 share the name "actb" and must collapse to one feature;
    # row 4 has no records at all and must not appear as a feature.
    mat <- Matrix::sparseMatrix(
        i = c(1L, 2L, 3L),
        j = c(1L, 1L, 1L),
        x = c(4, 6, 5),
        dims = c(4L, 2L)
    )
    gef <- .write_cellbin_gef(
        mat, file.path(tempdir(), "cellbin_dup.gef"),
        gene_names = c("actb", "cd8a", "actb", "empty")
    )
    out <- file.path(tempdir(), "cellbin_dup_out")
    on.exit(unlink(c(gef, out), recursive = TRUE), add = TRUE)

    pe <- storeWrite(parquetExprStore(path = out), cellbinGefInput(gef))

    expect_equal(pe@feat_ids, c("actb", "cd8a"))
    got <- .pe_as_matrix(pe)
    expect_equal(unname(got["actb", 1L]), 4 + 5)   # both rows summed
    expect_equal(unname(got["cd8a", 1L]), 6)
})


test_that("cellbinGefInput sums duplicate gene names split across chunks", {
    # The case a real gene table produces: the table is ordered by geneID, so
    # two rows sharing a geneName sit far apart and fall in different gene
    # chunks. Only the geneExp path reads gene chunks, so the file here has no
    # cellExp; the batch aggregate must still leave one row per (cell, gene).
    n_genes <- 10L
    mat <- Matrix::sparseMatrix(
        i = c(1L, 10L, 5L), j = c(1L, 1L, 2L), x = c(4, 5, 7),
        dims = c(n_genes, 3L)
    )
    names_v <- paste0("g", seq_len(n_genes))
    names_v[c(1L, 10L)] <- "shared"          # duplicated, maximally far apart

    gef <- .write_cellbin_gef(mat, file.path(tempdir(), "cellbin_split.gef"),
                              gene_names = names_v, cell_exp = FALSE)
    out <- file.path(tempdir(), "cellbin_split_out")
    on.exit(unlink(c(gef, out), recursive = TRUE), add = TRUE)

    # batch_genes = 1 guarantees rows 1 and 10 land in different chunks
    inp <- cellbinGefInput(gef, batch_genes = 1L)
    pe  <- storeWrite(parquetExprStore(path = out), inp)

    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    # one row per (cell, gene) -- the invariant that duplicate names break
    expect_false(any(duplicated(df[, c("row_id", "col_id")])))
    expect_equal(unname(.pe_as_matrix(pe)["shared", 1L]), 4 + 5)
})


test_that("cellbinGefInput writes cell-major parts and needs no cellExp to", {
    set.seed(11)
    mat <- Matrix::rsparsematrix(8L, 60L, density = 0.4,
        rand.x = function(n) as.double(rpois(n, 3) + 1))
    names_v <- paste0("g", 1:8)
    names_v[c(2L, 7L)] <- "dup"
    withr::local_options(list(giottodisk.gef_batch_rows = 30L))

    got <- list()
    for (ce in c(TRUE, FALSE)) {
        gef <- .write_cellbin_gef(mat, tempfile(fileext = ".gef"), gene_names = names_v,
                                  cell_exp = ce)
        out <- tempfile("cellbin_cm_")
        inp <- cellbinGefInput(gef, batch_genes = 1L)
        expect_identical(inp@params$has_cell_exp, ce)
        pe <- storeWrite(parquetExprStore(path = out), inp)
        pr <- .expect_cell_major(out)
        expect_gt(nrow(pr), 1L)
        got[[as.character(ce)]] <- .pe_as_matrix(pe)
        unlink(c(gef, out), recursive = TRUE)
    }
    expect_equal(got[["TRUE"]], got[["FALSE"]])
    ref <- rowsum(as.matrix(mat), names_v)[got[["TRUE"]]@Dimnames[[1L]], ]
    expect_equal(unname(as.matrix(got[["TRUE"]])), unname(ref))
})


# ---- bin -------------------------------------------------------------------

test_that("binGefInput + storeWrite round-trips a bin gef losslessly", {
    mat <- Matrix::sparseMatrix(
        i = c(1L, 1L, 2L, 2L, 3L),
        j = c(2L, 3L, 1L, 3L, 4L),
        x = c(3, 8, 1, 6, 2),
        dims = c(3L, 4L)
    )
    rownames(mat) <- c("aaa", "bbb", "ccc")

    fx  <- .write_bin_gef(mat, file.path(tempdir(), "bin_rt.gef"))
    out <- file.path(tempdir(), "bin_rt_out")
    on.exit(unlink(c(fx$path, out), recursive = TRUE), add = TRUE)

    inp <- binGefInput(fx$path, bin_size = "bin100")
    # cell identity is unknown until the stream has been driven
    expect_equal(inp@n_cells, 0L)
    expect_length(inp@cell_ids, 0L)

    pe <- storeWrite(parquetExprStore(path = out), inp)

    # bin_IDs are named in first-appearance order over the record stream (the
    # in-memory reader's numbering); the store positions them in grid order.
    appearance <- unique(fx$records$cell_col)
    expect_equal(pe@n_cells, length(appearance))
    expect_setequal(pe@cell_ids, paste0("bin_", seq_along(appearance)))

    src_col <- appearance[as.integer(sub("bin_", "", pe@cell_ids))]
    expected <- mat[, src_col, drop = FALSE]
    expect_equal(unname(as.matrix(.pe_as_matrix(pe))),
                 unname(as.matrix(expected)))
})

test_that("binGefInput positions bins in (y, x) grid order", {
    mat <- Matrix::sparseMatrix(
        i = c(1L, 2L, 1L, 2L, 1L, 2L), j = c(1L, 2L, 3L, 4L, 5L, 6L),
        x = as.double(1:6), dims = c(2L, 6L)
    )
    rownames(mat) <- c("aaa", "bbb")
    coords <- data.frame(x = c(30L, 10L, 20L, 10L, 30L, 20L),
                         y = c(5L, 5L, 1L, 9L, 1L, 9L))
    fx  <- .write_bin_gef(mat, file.path(tempdir(), "bin_grid.gef"), coords = coords)
    out <- file.path(tempdir(), "bin_grid_out")
    on.exit(unlink(c(fx$path, out), recursive = TRUE), add = TRUE)

    inp <- binGefInput(fx$path, bin_size = "bin100")
    pe  <- storeWrite(parquetExprStore(path = out), inp)

    bc <- inp@params$coord_env$bin_coords
    lab <- as.integer(sub("bin_", "", pe@cell_ids))
    at <- bc[match(lab, bc$bin_ID)]
    expect_equal(order(at$y, at$x), seq_len(nrow(at)))
    expect_equal(bc$pos[match(lab, bc$bin_ID)], seq_along(lab))
})

test_that("binGefInput addresses the prefixed geneExp/<bin_size> group", {
    # Real GEFs key on "bin100", not "100" -- Giotto's in-memory reader
    # hardcodes "geneExp/bin1/expression". A reader that strips the prefix
    # finds nothing on a real file, so pin the group name here.
    mat <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 3, dims = c(1L, 1L))
    rownames(mat) <- "aaa"
    fx <- .write_bin_gef(mat, file.path(tempdir(), "bin_key.gef"),
                         bin_key = "bin100")
    on.exit(unlink(fx$path), add = TRUE)

    expect_s4_class(binGefInput(fx$path, bin_size = "bin100"), "binGefInput")
    expect_error(binGefInput(fx$path, bin_size = "100"))
})

test_that("binGefInput publishes bin coordinates in first-appearance order", {
    # The parity that lets spatial locations come from the ingest stream:
    # the streamed map must equal what the in-memory reader derives with
    # unique(exprDT[, c("x","y")])[, bin_ID := .I] over the whole table.
    mat <- Matrix::sparseMatrix(
        i = c(1L, 1L, 2L, 2L, 3L, 3L),
        j = c(4L, 2L, 1L, 4L, 3L, 2L),
        x = c(1, 2, 3, 4, 5, 6),
        dims = c(3L, 4L)
    )
    rownames(mat) <- c("aaa", "bbb", "ccc")

    fx  <- .write_bin_gef(mat, file.path(tempdir(), "bin_coord.gef"))
    out <- file.path(tempdir(), "bin_coord_out")
    on.exit(unlink(c(fx$path, out), recursive = TRUE), add = TRUE)

    inp <- binGefInput(fx$path, bin_size = "bin100")
    storeWrite(parquetExprStore(path = out), inp)

    got <- inp@params$coord_env$bin_coords
    expect_s3_class(got, "data.table")

    # independent ground truth: the full on-disk record table, deduped
    full <- data.table::setDT(
        rhdf5::h5read(fx$path, "geneExp/bin100/expression")
    )
    rhdf5::h5closeAll()
    want <- unique(full[, c("x", "y")], by = c("x", "y"))
    want[, bin_ID := seq_len(nrow(want))]

    expect_equal(nrow(got), nrow(want))
    expect_equal(got$bin_ID, want$bin_ID)
    expect_equal(got$x, want$x)
    expect_equal(got$y, want$y)
})

test_that("binGefInput does not publish coordinates from a partial pass", {
    # A half-consumed iterator has only some of the bins. Publishing that
    # would yield spatial locations for a subset of the store's columns.
    mat <- Matrix::sparseMatrix(
        i = c(1L, 2L, 3L, 4L), j = c(1L, 2L, 3L, 4L),
        x = c(1, 2, 3, 4), dims = c(4L, 4L)
    )
    rownames(mat) <- paste0("g", 1:4)
    fx <- .write_bin_gef(mat, file.path(tempdir(), "bin_partial.gef"))
    on.exit(unlink(fx$path), add = TRUE)

    inp <- binGefInput(fx$path, bin_size = "bin100", batch_genes = 1L)
    itr <- storeRead(inp)
    itr$next_batch()      # one chunk only, then abandon
    itr$close()

    expect_null(inp@params$coord_env$bin_coords)
})

# ---- reader wiring ---------------------------------------------------------

# helper: a minimal Stereo-seq output directory holding one bin gef, laid out
# the way .stereoseq_detect_paths() expects to find it.
.write_ss_dir <- function(mat, dir, coords = NULL) {
    dir.create(file.path(dir, "feature_expression"), recursive = TRUE,
               showWarnings = FALSE)
    .write_bin_gef(
        mat, file.path(dir, "feature_expression", "SAMPLE.tissue.gef"),
        bin_key = "bin100", coords = coords
    )
}

test_that("createGiottoStereoSeqObjectBin(backend =) builds a disk-backed object", {
    skip_if_not_installed("Giotto")

    mat <- Matrix::sparseMatrix(
        i = c(1L, 1L, 2L, 2L, 3L), j = c(2L, 3L, 1L, 3L, 4L),
        x = c(3, 8, 1, 6, 2), dims = c(3L, 4L)
    )
    rownames(mat) <- c("aaa", "bbb", "ccc")
    coords <- data.frame(x = (1:4) * 5L, y = rev(1:4) * 7L)

    ss   <- file.path(tempdir(), "ss_reader")
    proj <- file.path(tempdir(), "ss_reader_proj")
    unlink(c(ss, proj), recursive = TRUE)
    dir.create(proj, recursive = TRUE)
    fx <- .write_ss_dir(mat, ss, coords = coords)
    on.exit(unlink(c(ss, proj), recursive = TRUE), add = TRUE)

    g <- Giotto::createGiottoStereoSeqObjectBin(
        stereoseq_dir = ss,
        bin_size      = "bin100",
        backend       = gDirSource(proj),
        load_image    = FALSE,
        load_mask     = FALSE,
        load_binpoints = FALSE,
        verbose       = FALSE
    )

    ex <- GiottoClass::getExpression(g, spat_unit = "bin100",
                                     feat_type = "rna", values = "raw",
                                     output = "exprObj")
    # the point of `backend =`: expression is on disk, not a dgCMatrix
    expect_s4_class(ex[], "parquetExprStore")
    expect_equal(dim(ex[]), c(3, 4))

    # spatial locations came off the ingest stream, so their IDs must be the
    # store's own columns -- not a separately-derived numbering
    sl <- GiottoClass::getSpatialLocations(g, spat_unit = "bin100",
                                           output = "data.table")
    expect_setequal(sl$cell_ID, ex[]@cell_ids)
    # rows in store column order, which is not bin_ID order
    expect_identical(sl$cell_ID, ex[]@cell_ids)

    appearance <- unique(fx$records$cell_col)
    expected <- data.frame(
        cell_ID = paste0("bin_", seq_along(appearance)),
        sdimx   = as.integer(coords$x[appearance]),
        sdimy   = as.integer(-coords$y[appearance])
    )
    got <- as.data.frame(sl)[order(as.data.frame(sl)$cell_ID), ]
    expected <- expected[order(expected$cell_ID), ]
    expect_equal(got$sdimx, expected$sdimx)
    expect_equal(got$sdimy, expected$sdimy)   # negative_y default applied
})

test_that("the disk reader's load_expression matches the in-memory shape", {
    skip_if_not_installed("Giotto")
    # A reader used piecewise must be substitutable with `backend =` set or
    # unset, so load_expression() returns list(exprObj) either way.
    mat <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(3, 4),
                                dims = c(2L, 2L))
    rownames(mat) <- c("aaa", "bbb")
    ss   <- file.path(tempdir(), "ss_shape")
    proj <- file.path(tempdir(), "ss_shape_proj")
    unlink(c(ss, proj), recursive = TRUE); dir.create(proj, recursive = TRUE)
    .write_ss_dir(mat, ss)
    on.exit(unlink(c(ss, proj), recursive = TRUE), add = TRUE)

    rd <- Giotto::importStereoSeq(ss, type = "bin", bin_size = "bin100",
                                  backend = gDirSource(proj))
    rm_ <- Giotto::importStereoSeq(ss, type = "bin", bin_size = "bin100")
    ed <- rd@calls$load_expression(verbose = FALSE)
    em <- rm_@calls$load_expression(verbose = FALSE)

    expect_type(ed, "list")
    expect_length(ed, length(em))
    expect_s4_class(ed[[1L]], "exprObj")
    expect_s4_class(ed[[1L]][], "parquetExprStore")
})

test_that("binGefInput writes several cell-major parts on a small budget", {
    set.seed(7)
    mat <- Matrix::rsparsematrix(6L, 40L, density = 0.5,
        rand.x = function(n) as.double(rpois(n, 4) + 1))
    rownames(mat) <- paste0("g", 1:6)
    coords <- data.frame(x = rep(1:8, 5) * 10L, y = rep(1:5, each = 8) * 10L)

    fx  <- .write_bin_gef(mat, file.path(tempdir(), "bin_chunk.gef"), coords = coords)
    out <- file.path(tempdir(), "bin_chunk_out")
    on.exit(unlink(c(fx$path, out), recursive = TRUE), add = TRUE)

    withr::local_options(list(giottodisk.gef_batch_rows = 25L))
    inp <- binGefInput(fx$path, bin_size = "bin100", batch_genes = 1L)
    pe  <- storeWrite(parquetExprStore(path = out), inp)

    pr <- .expect_cell_major(out)
    expect_gt(nrow(pr), 1L)
    expect_equal(max(pr$hi), pe@n_cells)

    appearance <- unique(fx$records$cell_col)
    src_col <- appearance[as.integer(sub("bin_", "", pe@cell_ids))]
    expect_equal(unname(as.matrix(.pe_as_matrix(pe))),
                 unname(as.matrix(mat[, src_col, drop = FALSE])))
    # grid order holds across stripes, not only within one
    at <- coords[src_col, ]
    expect_equal(order(at$y, at$x), seq_len(nrow(at)))
})

test_that("binGefInput derives the y range when the file does not record it", {
    set.seed(8)
    mat <- Matrix::rsparsematrix(5L, 30L, density = 0.5,
        rand.x = function(n) as.double(rpois(n, 4) + 1))
    rownames(mat) <- paste0("g", 1:5)
    coords <- data.frame(x = rep(1:6, 5) * 10L, y = rep(1:5, each = 6) * 10L)
    fx  <- .write_bin_gef(mat, file.path(tempdir(), "bin_noattr.gef"),
                          coords = coords, y_attrs = FALSE)
    out <- file.path(tempdir(), "bin_noattr_out")
    on.exit(unlink(c(fx$path, out), recursive = TRUE), add = TRUE)

    withr::local_options(list(giottodisk.gef_batch_rows = 20L))
    inp <- binGefInput(fx$path, bin_size = "bin100", batch_genes = 1L)
    expect_null(inp@params$y_range)
    pe <- storeWrite(parquetExprStore(path = out), inp)

    expect_gt(nrow(.expect_cell_major(out)), 1L)
    appearance <- unique(fx$records$cell_col)
    src_col <- appearance[as.integer(sub("bin_", "", pe@cell_ids))]
    expect_equal(unname(as.matrix(.pe_as_matrix(pe))),
                 unname(as.matrix(mat[, src_col, drop = FALSE])))
})

test_that("binGefInput sums duplicate gene names split across chunks", {
    mat <- Matrix::sparseMatrix(i = c(1L, 4L, 2L), j = c(1L, 1L, 2L),
        x = c(4, 5, 7), dims = c(4L, 2L))
    rownames(mat) <- c("shared", "g2", "g3", "shared")
    fx  <- .write_bin_gef(mat, file.path(tempdir(), "bin_dup.gef"))
    out <- file.path(tempdir(), "bin_dup_out")
    on.exit(unlink(c(fx$path, out), recursive = TRUE), add = TRUE)

    pe <- storeWrite(parquetExprStore(path = out),
                     binGefInput(fx$path, bin_size = "bin100", batch_genes = 1L))
    df <- as.data.frame(dplyr::collect(storeRead(pe)))
    expect_false(any(duplicated(df[, c("row_id", "col_id")])))
    got <- .pe_as_matrix(pe)
    first <- paste0("bin_", 1L)   # cell column 1 appears first
    expect_equal(unname(got["shared", first]), 4 + 5)
})

