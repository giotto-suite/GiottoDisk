# CosMx multiomic exports (protein_dir) and fovs on every input. Fixtures:
# helper-cosmx-fixture.R. Cells of the n_fov1 = 300 fixture:
#   RNA:      FOV 1 cells 1-300, FOV 2 cells 1-3       (303)
#   protein:  RNA without c_1_1_5, plus c_1_2_4        (303)
#   polygons: FOV 1 cells 1-300, FOV 2 cells 1, 2, 4   (303)
# so RNA, protein and a polygon: 301 cells.

.cosmx_p_skip <- function() {
    skip_if_not_installed("arrow")
    skip_if_not_installed("Giotto")
}

.cosmx_reader <- function(d, ...) {
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    suppressWarnings(importCosMxDisk(d, backend = p, version = "v6", ...))
}

# value plus every warning message (line breaks of wrap_txt() undone), so
# that one warning can be checked among the others (tilework, image folders)
.with_warnings <- function(expr) {
    w <- character(0L)
    value <- withCallingHandlers(expr, warning = function(cnd) {
        w <<- c(w, gsub("\\s+", " ", conditionMessage(cnd)))
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = w)
}

.cosmx_mo_obj <- function(d, pd, ..., fovs = NULL) {
    r <- .cosmx_reader(d, protein_dir = pd, fovs = fovs)
    .with_warnings(suppressMessages(
        r$create_gobject(load_images = NULL, verbose = FALSE, ...)
    ))
}

.expr_ids <- function(g, feat_type) {
    colnames(GiottoClass::getExpression(g, feat_type = feat_type,
                                        values = "raw")[])
}

.expr_dense <- function(g, feat_type) {
    pe <- GiottoClass::getExpression(g, feat_type = feat_type,
                                     values = "raw")[]
    as.matrix(storeRead(pe, output = "dgcmatrix", max_rows = Inf))
}

# protein matrix ####

test_that("load_protein splits targets and isotypes and drops Channel-*", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    r <- .cosmx_reader(d, protein_dir = pd)
    st <- r$load_protein(output = "store", verbose = FALSE)

    expect_named(st, c("protein", "protein_neg"))
    expect_identical(st$protein@feat_ids, c("CD3", "GENEA", "NF-kB p65"))
    expect_identical(st$protein_neg@feat_ids, c("Ms IgG1", "Rb IgG"))
    vals <- attr(pd, "values")
    expect_identical(st$protein@cell_ids, rownames(vals))
})

test_that("protein values are the file's decimals, unchanged", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    r <- .cosmx_reader(d, protein_dir = pd)
    st <- r$load_protein(output = "store", verbose = FALSE)

    ref <- data.table::fread(file.path(pd, "Fixture_exprMat_file.csv"))
    ref <- t(as.matrix(ref[, c("CD3", "GENEA", "NF-kB p65", "Ms IgG1",
                               "Rb IgG")]))
    got <- rbind(
        as.matrix(storeRead(st$protein, output = "dgcmatrix")),
        as.matrix(storeRead(st$protein_neg, output = "dgcmatrix"))
    )
    expect_identical(unname(got), unname(ref))
})

test_that("isotypes fall back to the IgG pattern without a plex file", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d, plex = FALSE)
    r <- .cosmx_reader(d, protein_dir = pd)
    st <- r$load_protein(output = "store", verbose = FALSE)
    expect_identical(st$protein_neg@feat_ids, c("Ms IgG1", "Rb IgG"))
})

test_that("a decimal .csv.gz matrix is read by the CSV route", {
    .cosmx_p_skip()
    skip_if_not_installed("cosmxscan")
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d, gz = TRUE)
    path <- file.path(pd, "Fixture_exprMat_file.csv.gz")
    expect_false(GiottoDisk:::.cosmx_values_are_integer(path))

    r <- .cosmx_reader(d, protein_dir = pd)
    st <- r$load_protein(output = "store", verbose = FALSE)
    expect_length(st$protein@cell_ids, 303L)

    # scanner = "auto" with nothing to drop: cosmxscan would stop at the
    # first `.` (UPSTREAM_ISSUES #16)
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    whole <- GiottoDisk:::.cosmx_expression_disk(path, gDirSource(path = p),
        feat_type = "x", split_keyword = list(), output = "store",
        verbose = FALSE)
    expect_true("Channel-DNA" %in% whole$x@feat_ids)
})

# create_gobject ####

test_that("create_gobject keeps cells with RNA, protein and a polygon", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    res <- .cosmx_mo_obj(d, pd)
    g <- res$value

    ids <- GiottoClass::spatIDs(g)
    expect_length(ids, 301L)
    expect_false(any(c("c_1_1_5", "c_1_2_3", "c_1_2_4") %in% ids))
    expect_true(any(grepl(paste(
        "3 cells are not loaded: 1 with protein but no RNA expression",
        "\\(e.g. c_1_2_4\\); 1 with RNA but no protein expression",
        "\\(e.g. c_1_1_5\\); 1 with expression but no polygon",
        "\\(e.g. c_1_2_3\\)"), res$warnings)))
    # c_1_2_4 has protein: its polygon is not an orphan
    expect_false(any(grepl("polygons belong to cells", res$warnings)))
    expect_setequal(GiottoClass::list_expression(g)$feat_type,
        c("rna", "negprobes", "falsecode", "protein", "protein_neg"))
    poly_ids <- GiottoClass::spatIDs(GiottoClass::getPolygonInfo(g,
        return_giottoPolygon = TRUE))
    expect_setequal(poly_ids, ids)
})

test_that("a cell with NA protein values is left out as missing from it", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    f <- file.path(pd, "Fixture_exprMat_file.csv")
    x <- data.table::fread(f)
    x[["CD3"]][x$fov == 1L & x$cell_ID == 7L] <- NA
    data.table::fwrite(x, f)
    res <- .cosmx_mo_obj(d, pd)
    ids <- GiottoClass::spatIDs(res$value)
    expect_length(ids, 300L)
    expect_false("c_1_1_7" %in% ids)
    expect_true(any(grepl(paste("1 cell in Fixture_exprMat_file.csv has NA",
        "values and is left out, as cells missing from it \\(e.g. c_1_1_7\\)"),
        res$warnings)))
    expect_true(any(grepl("2 with RNA but no protein expression",
                          res$warnings)))
    expect_identical(.expr_ids(res$value, "protein"), ids)
    pe <- GiottoClass::getExpression(res$value, feat_type = "protein",
                                     values = "raw")[]
    got <- .expr_dense(res$value, "protein")
    dimnames(got) <- list(pe@feat_ids, pe@cell_ids)
    want <- attr(pd, "values")
    expect_equal(unname(got["CD3", ids]), unname(want[ids, "CD3"]))
})

test_that("every feature type has the cells in one order", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    # protein rows in another order than RNA: the reader must reorder
    f <- file.path(pd, "Fixture_exprMat_file.csv")
    x <- data.table::fread(f)
    data.table::fwrite(x[rev(seq_len(nrow(x))), ], f)
    g <- .cosmx_mo_obj(d, pd)$value

    ids <- GiottoClass::spatIDs(g)
    for (ft in c("rna", "negprobes", "falsecode", "protein",
                 "protein_neg")) {
        expect_identical(.expr_ids(g, ft), ids, label = ft)
    }
    vals <- attr(pd, "values")
    m <- .expr_dense(g, "protein")
    expect_identical(unname(m["CD3", ]), unname(vals[ids, "CD3"]))
    m_neg <- .expr_dense(g, "protein_neg")
    expect_identical(unname(m_neg["Ms IgG1", ]), unname(vals[ids, "Ms IgG1"]))
})

test_that("protein metadata is renamed and kept apart from the RNA's", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    g <- .cosmx_mo_obj(d, pd)$value

    rm <- GiottoClass::pDataDT(g, feat_type = "rna")
    pm <- GiottoClass::pDataDT(g, feat_type = "protein")
    rfile <- data.table::fread(file.path(d, "Fixture_metadata_file.csv"))
    pfile <- data.table::fread(file.path(pd, "Fixture_metadata_file.csv"))

    expect_identical(rm$nCount_RNA,
                     rfile$nCount_RNA[match(rm$cell_ID, rfile$cell)])
    expect_false(any(c("nCount_RNA", "nFeature_RNA", "nCount_negprobes",
                       "Panel") %in% names(pm)))
    expect_true(all(c("fov", "Area.um2", "nCount_protein",
                      "nFeature_protein", "nCount_protein_neg",
                      "nFeature_protein_neg", "Panel_protein",
                      "assay_type_protein") %in% names(pm)))
    i <- match(pm$cell_ID, pfile$cell)
    expect_identical(pm$nCount_protein, pfile$nCount_RNA[i])
    expect_identical(pm$nCount_protein_neg, pfile$nCount_negprobes[i])
    expect_type(pm$fov, "integer")
    expect_identical(pm$fov, rm$fov[match(pm$cell_ID, rm$cell_ID)])
})

test_that("a protein-only object keeps the protein cells with a polygon", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    res <- .cosmx_mo_obj(d, pd, load_expression = FALSE)
    g <- res$value

    ids <- GiottoClass::spatIDs(g)
    expect_length(ids, 302L)
    expect_true("c_1_2_4" %in% ids)
    expect_setequal(GiottoClass::list_expression(g)$feat_type,
                    c("protein", "protein_neg"))
    pm <- GiottoClass::pDataDT(g, feat_type = "protein")
    expect_true(all(c("fov", "Area.um2", "nCount_protein") %in% names(pm)))
    expect_true(any(grepl("1 with expression but no polygon",
                          res$warnings)))
})

test_that("load_protein = TRUE without a protein matrix is an error", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    r <- .cosmx_reader(d)
    expect_error(suppressWarnings(r$create_gobject(load_protein = TRUE,
        load_images = NULL, verbose = FALSE)), "no protein exprMat_file")
})

test_that("one protein feat_type keeps the isotypes with the targets", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    r <- .cosmx_reader(d, protein_dir = pd)
    st <- r$load_protein(feat_type = "protein", output = "store",
                         verbose = FALSE)
    expect_named(st, "protein")
    expect_setequal(st$protein@feat_ids,
                    c("CD3", "GENEA", "NF-kB p65", "Ms IgG1", "Rb IgG"))
    expect_error(r$load_protein(feat_type = "protein",
                                split_keyword = list("IgG"),
                                output = "store", verbose = FALSE),
                 "1 feat_type\\(s\\) given for 2")
})

test_that("a protein matrix that will not split as expected warns", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    r <- .cosmx_reader(d, protein_dir = pd)
    expect_warning(r$load_protein(split_keyword = list("^NoSuchIsotype$"),
                                  output = "store", verbose = FALSE),
                   "no isotype control")
    expect_warning(r$load_protein(path = file.path(d, "Fixture_exprMat_file.csv"),
                                  output = "store", verbose = FALSE),
                   "is `protein_dir` the RNA export")
    expect_no_warning(r$load_protein(output = "store", verbose = FALSE))
})

test_that("protein_drop = NULL keeps the Channel columns as targets", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    res <- .cosmx_mo_obj(d, pd, protein_drop = NULL)
    expect_true("Channel-DNA" %in% GiottoClass::featIDs(res$value,
                                                        feat_type = "protein"))
})

test_that("several matching files in protein_dir warn and take the first", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    file.copy(file.path(pd, "Fixture_exprMat_file.csv"),
              file.path(pd, "Other_exprMat_file.csv"))
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    expect_warning(
        importCosMxDisk(d, backend = p, version = "v6", protein_dir = pd),
        "2 exprMat_file files in protein_dir; using Fixture"
    )
})

test_that("an empty split group is dropped by name, not shifted", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    r <- .cosmx_reader(d, protein_dir = pd)
    st <- r$load_protein(feat_type = c("protein", "protein_neg", "protein_morph"),
                         split_keyword = list("^NoSuch", "^Channel-"),
                         drop_cols = NULL, output = "store", verbose = FALSE)
    expect_named(st, c("protein", "protein_morph"))
    expect_true(length(st$protein_morph@feat_ids) > 0L)
    expect_true(all(grepl("^Channel-", st$protein_morph@feat_ids)))
    expect_false(any(grepl("^Channel-", st$protein@feat_ids)))

    rl <- r$load_expression(feat_type = c("rna", "negprobes", "falsecode"),
                            split_keyword = list("^NoSuch", "^SystemControl"),
                            output = "store", verbose = FALSE)
    expect_named(rl, c("rna", "falsecode"))
    expect_true(all(grepl("^SystemControl", rl$falsecode@feat_ids)))

    expect_error(r$load_protein(feat_type = c("protein", "protein_neg", "x"),
                                output = "store", verbose = FALSE),
                 "3 feat_type\\(s\\) given for 2")
})

test_that("RNA loaded only as transcripts takes its cells from the RNA metadata", {
    .cosmx_p_skip()
    skip_if_not("protein_dir" %in% names(formals(Giotto::importCosMx)))
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    .mk_cosmx_tx(d, n_blank = 10L)
    args <- list(d, protein_dir = pd, version = "v6", poly_pref = "csv",
        feat_type = c("rna", "negprobes", "falsecode"),
        split_keyword = list("^Negative", "^SystemControl"),
        load_expression = FALSE, load_transcripts = TRUE, load_images = NULL,
        cores = 1L, verbose = FALSE)
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    for (backend in list(NULL, p)) {
        res <- .with_warnings(suppressMessages(do.call(
            Giotto::createGiottoCosMxObject, c(args, list(backend = backend)))))
        ids <- GiottoClass::spatIDs(res$value)
        expect_length(ids, 301L)
        expect_false(any(c("c_1_1_5", "c_1_2_3", "c_1_2_4") %in% ids))
        expect_true(any(grepl("1 with protein but no RNA expression",
                              res$warnings)))
        if (is.null(backend)) {
            rm <- GiottoClass::pDataDT(res$value, feat_type = "rna")
            expect_identical(nrow(rm), length(ids))
        }
    }
})

test_that("the in-memory and disk readers build the same multiomic object", {
    .cosmx_p_skip()
    skip_if_not("protein_dir" %in% names(formals(Giotto::importCosMx)))
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    args <- list(d, protein_dir = pd, version = "v6", poly_pref = "csv",
        feat_type = c("rna", "negprobes", "falsecode"),
        split_keyword = list("^Negative", "^SystemControl"),
        load_expression = TRUE, load_transcripts = FALSE, load_images = NULL,
        cores = 1L, verbose = FALSE)
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    mem <- .with_warnings(suppressMessages(
        do.call(Giotto::createGiottoCosMxObject, args)))
    disk <- .with_warnings(suppressMessages(
        do.call(Giotto::createGiottoCosMxObject, c(args, list(backend = p)))))
    gm <- mem$value
    gd <- disk$value
    expect_false(is.null(gd@source))
    expect_identical(GiottoClass::spatIDs(gd), GiottoClass::spatIDs(gm))
    fts <- GiottoClass::list_expression(gm)$feat_type
    expect_setequal(GiottoClass::list_expression(gd)$feat_type, fts)
    for (ft in fts) {
        m <- as.matrix(GiottoClass::getExpression(gm, feat_type = ft,
            values = "raw", output = "matrix"))
        pe <- GiottoClass::getExpression(gd, feat_type = ft, values = "raw")[]
        x <- as.matrix(storeRead(pe, output = "dgcmatrix", max_rows = Inf))
        dimnames(x) <- list(pe@feat_ids, pe@cell_ids)
        expect_equal(x[rownames(m), colnames(m), drop = FALSE], m, label = ft)
    }
    not_loaded <- function(r) {
        sub("^\\[(CosMxDiskReader|CosmxReader)\\] ", "",
            grep("not loaded", r$warnings, value = TRUE))
    }
    expect_identical(not_loaded(disk), not_loaded(mem))
    pm <- GiottoClass::pDataDT(gm, feat_type = "protein")
    pd_ <- GiottoClass::pDataDT(gd, feat_type = "protein")
    shared <- setdiff(intersect(names(pd_), names(pm)), c("cell_ID", "fov"))
    shared <- shared[vapply(shared, function(col) is.numeric(pm[[col]]),
                            logical(1L))]
    expect_true("nCount_protein" %in% shared)
    for (col in shared) {
        expect_equal(as.numeric(pd_[[col]][match(pm$cell_ID, pd_$cell_ID)]),
                     as.numeric(pm[[col]]), label = col)
    }
})

test_that("an empty cell selection stops with the reason", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    expect_error(.cosmx_mo_obj(d, pd, fovs = 99),
                 "no cells in Fixture_exprMat_file.csv for fovs 99")

    f <- file.path(pd, "Fixture_exprMat_file.csv")
    x <- data.table::fread(f)
    x$cell_ID <- x$cell_ID + 1000L
    data.table::fwrite(x, f)
    expect_error(.cosmx_mo_obj(d, pd), "share no cell")
})

# fovs ####

test_that("fovs restrict expression, protein, polygons and metadata", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    res <- .cosmx_mo_obj(d, pd, fovs = 2)
    g <- res$value

    ids <- c("c_1_2_1", "c_1_2_2")
    expect_identical(GiottoClass::spatIDs(g), ids)
    for (ft in c("rna", "protein", "protein_neg")) {
        expect_identical(.expr_ids(g, ft), ids, label = ft)
    }
    poly_ids <- GiottoClass::spatIDs(GiottoClass::getPolygonInfo(g,
        return_giottoPolygon = TRUE))
    expect_setequal(poly_ids, ids)
    expect_identical(GiottoClass::pDataDT(g, feat_type = "rna")$cell_ID, ids)
    expect_identical(GiottoClass::pDataDT(g, feat_type = "protein")$cell_ID,
                     ids)
    # c_1_2_4 (protein only) and c_1_2_3 (no polygon) are in FOV 2
    expect_true(any(grepl(paste("2 cells are not loaded: 1 with protein but",
        "no RNA expression \\(e.g. c_1_2_4\\); 1 with expression but no",
        "polygon"), res$warnings)))
})

test_that("the polygon CSV is filtered by FOV in file order before conversion", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    f <- file.path(d, "Fixture-polygons.csv")
    ref <- data.table::fread(f)
    sub <- GiottoDisk:::.cosmx_csv_fov_subset(f, 2L, chunk_lines = 7L)
    expect_identical(data.table::fread(sub), ref[ref$fov == 2L])
    nofov <- tempfile(fileext = ".csv")
    data.table::fwrite(data.frame(a = 1:3, b = 4:6), nofov)
    expect_null(GiottoDisk:::.cosmx_csv_fov_subset(nofov, 2L))
})

test_that("fovs restrict the polygons before the orphan count", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    r <- .cosmx_reader(d, fovs = 2)
    res <- .with_warnings(suppressMessages(r$create_gobject(
        load_images = NULL, verbose = FALSE)))
    expect_true(any(grepl("1 of 3 polygons belong to cells", res$warnings)))
    expect_identical(GiottoClass::spatIDs(res$value),
                     c("c_1_2_1", "c_1_2_2"))
})

test_that("a decimal after the first row stops cosmxscan with advice", {
    .cosmx_p_skip()
    skip_if_not_installed("cosmxscan")
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    f <- file.path(d, "Fixture_exprMat_file.csv")
    x <- data.table::fread(f)
    x[["GENEA"]] <- as.numeric(x[["GENEA"]])
    x[["GENEA"]][50L] <- 2.5
    data.table::fwrite(x, paste0(f, ".gz"))
    unlink(f)
    r <- .cosmx_reader(d)
    expect_error(r$load_expression(feat_type = "x", split_keyword = list(),
                                   output = "store", verbose = FALSE),
                 "non-integer values after its first row.*scanner = \"csv\"")
    st <- r$load_expression(feat_type = "x", split_keyword = list(),
                            scanner = "csv", output = "store", verbose = FALSE)
    got <- as.matrix(storeRead(st$x, output = "dgcmatrix"))
    expect_identical(unname(got[match("GENEA", st$x@feat_ids), 50L]), 2.5)
    p <- tempfile("cosmx_proj_")
    dir.create(p)
    expect_error(GiottoDisk:::.cosmx_expression_disk(paste0(f, ".gz"),
        gDirSource(path = p), feat_type = "x", split_keyword = list(),
        drop_cols = "GENEB", scanner = "cosmxscan", output = "store",
        verbose = FALSE), "cannot skip the columns")
})

test_that("the cosmxscan header is read without quotes", {
    .cosmx_p_skip()
    skip_if_not_installed("cosmxscan")
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    f <- file.path(d, "Fixture_exprMat_file.csv")
    x <- data.table::fread(f)
    data.table::fwrite(x, paste0(f, ".gz"), quote = TRUE)
    unlink(f)
    inp <- GiottoDisk:::.cosmx_scan_input(paste0(f, ".gz"))
    expect_identical(inp@feat_ids, c("GENEA", "GENEB", "Negative1",
                                     "SystemControl1"))
})

test_that("the cosmxscan route honours fovs and drops background rows", {
    .cosmx_p_skip()
    skip_if_not_installed("cosmxscan")
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    ref <- data.table::fread(file.path(d, "Fixture_exprMat_file.csv"))
    .cosmx_gzip_expr(d, background = TRUE)
    feats <- c("GENEA", "GENEB", "Negative1", "SystemControl1")

    r <- .cosmx_reader(d, fovs = 2)
    expect_message(
        st <- r$load_expression(feat_type = "x", split_keyword = list(),
                                output = "store", verbose = TRUE),
        "scanning .gz")
    pe <- st$x
    expect_identical(pe@cell_ids, c("c_1_2_1", "c_1_2_2", "c_1_2_3"))
    want <- t(as.matrix(ref[ref$fov == 2L, ..feats]))
    got <- as.matrix(storeRead(pe, output = "dgcmatrix"))
    expect_identical(unname(got), unname(want) + 0)

    r_all <- .cosmx_reader(d)
    pe_all <- r_all$load_expression(feat_type = "x", split_keyword = list(),
                                    output = "store", verbose = FALSE)$x
    expect_length(pe_all@cell_ids, 303L)
    expect_false(any(grepl("_0$", pe_all@cell_ids)))
    got_all <- as.matrix(storeRead(pe_all, output = "dgcmatrix"))
    expect_identical(unname(got_all),
                     unname(t(as.matrix(ref[, ..feats]))) + 0)

    r_none <- .cosmx_reader(d, fovs = 99)
    expect_error(r_none$load_expression(output = "store", verbose = FALSE),
                 "no cells in Fixture_exprMat_file.csv.gz for fovs 99")
})

test_that("the scan filter renumbers rows across chunks", {
    skip_if_not_installed("cosmxscan")
    d <- .mk_cosmx_disk_dir(n_fov1 = 30L)
    ref <- data.table::fread(file.path(d, "Fixture_exprMat_file.csv"))
    path <- .cosmx_gzip_expr(d, background = TRUE)
    feats <- c("GENEA", "GENEB", "Negative1", "SystemControl1")

    # 3-row chunks: some hold only FOV 1 rows (all dropped), one holds the
    # FOV 2 background row and FOV 1 rows
    for (fv in list(2L, c(1L, 2L))) {
        inp <- GiottoDisk:::.cosmx_scan_input(path, fovs = fv,
                                              batch_rows = 3L)
        p <- tempfile("cosmx_proj_")
        dir.create(p)
        pe <- sourceWrite(gDirSource(path = p), inp,
                          store_type = "parquetExpr")
        sel <- ref$fov %in% fv
        expect_identical(pe@cell_ids,
            sprintf("c_1_%d_%d", ref$fov[sel], ref$cell_ID[sel]))
        got <- as.matrix(storeRead(pe, output = "dgcmatrix"))
        expect_identical(unname(got),
                         unname(t(as.matrix(ref[sel, ..feats]))) + 0)
    }
})

# reader ####

test_that("protein files are found by pattern and protein_dir survives $<-", {
    .cosmx_p_skip()
    d <- .mk_cosmx_disk_dir(n_fov1 = 300L)
    pd <- .mk_cosmx_protein_dir(d)
    r <- .cosmx_reader(d, protein_dir = pd)
    expect_identical(r$protein_dir, as.character(pd))
    expect_identical(basename(r@paths$protein_plex_path), "plex-fixture.txt")
    expect_identical(basename(r@paths$protein_meta_path),
                     "Fixture_metadata_file.csv")

    suppressWarnings(r$fovs <- 2)
    expect_identical(r$protein_dir, as.character(pd))
    st <- r$load_protein(output = "store", verbose = FALSE)
    expect_identical(st$protein@cell_ids,
                     c("c_1_2_1", "c_1_2_2", "c_1_2_3", "c_1_2_4"))

    r2 <- .cosmx_reader(d)
    expect_null(r2@paths$protein_expr_path)
    suppressWarnings(r2$protein_dir <- pd)
    expect_identical(basename(r2@paths$protein_expr_path),
                     "Fixture_exprMat_file.csv")
    suppressWarnings(r2$protein_dir <- NULL)
    expect_length(r2$protein_dir, 0L)
    expect_null(r2@paths$protein_expr_path)
})
