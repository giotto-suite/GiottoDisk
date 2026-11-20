# AnnData .h5ad inputs and loaders (R/h5ad.R), against the synthetic file of
# helper-h5ad-fixture.R. Every read is compared with the source R objects,
# never with a second read through the same path.

skip_if_not_installed("hdf5r")

proj <- withr::local_tempdir()
gsrc <- suppressMessages(gDirSource(proj))


# expression ####

test_that("CSR X streams into a cell-major store equal to the source", {
    fx <- make_h5ad(format = "csr")
    inp <- h5adInput(fx$path, batch_records = 40)
    expect_identical(inp@cell_ids, fx$cells)
    expect_identical(inp@feat_ids, .disambiguate_feat_ids(fx$genes))
    pe <- sourceWrite(gsrc, inp, store_type = "parquetExpr")
    expect_equal(dim(pe), c(ncol(fx$M), nrow(fx$M)))
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
    expect_true(.parts_cell_major(pe))
})

test_that("CSC X is reordered to the same cell-major store", {
    fx <- make_h5ad(format = "csc", seed = 2L)
    pe <- sourceWrite(gsrc, h5adInput(fx$path, batch_records = 40),
                      store_type = "parquetExpr")
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
    expect_true(.parts_cell_major(pe))
})

test_that("the session option pins the batch size; the argument wins", {
    fx <- make_h5ad(seed = 3L)
    withr::local_options(giottodisk.h5_batch_rows = 40)
    pe <- sourceWrite(gsrc, h5adInput(fx$path), store_type = "parquetExpr")
    expect_true(.parts_cell_major(pe))                      # several parts
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
    expect_identical(.h5ad_max_rows(NULL), 40)
    expect_identical(.h5ad_max_rows(1e3), 1e3)
})

test_that("layers and raw resolve; a dense layer streams its non-zeros", {
    fx <- make_h5ad()
    pe <- sourceWrite(gsrc, h5adInput(fx$path, layer = "counts"),
                      store_type = "parquetExpr")
    expect_equal(.store_triplets(pe), .source_triplets(fx$M * 2))
    raw <- h5adInput(fx$path, layer = "raw")
    expect_identical(raw@feat_ids, fx$genes[1:10])
    expect_identical(raw@n_genes, 10L)
    dn <- sourceWrite(gsrc, h5adInput(fx$path, layer = "dense", batch_records = 100),
                      store_type = "parquetExpr")
    expect_true(.parts_cell_major(dn))                      # several blocks
    expect_equal(.store_triplets(dn), .source_triplets(fx$M))
    expect_error(h5adInput(fx$path, layer = "nope"), "Available: X, raw, counts, dense")
})


# metadata, spatial, dimreds ####

test_that("obs decodes categoricals, nullables, strings and booleans in order", {
    fx <- make_h5ad()
    cm <- h5adCellMeta(fx$path)
    dt <- cm[]
    expect_identical(names(dt),
        c("cell_ID", "n_counts", "group", "label", "flag", "n_extra"))
    expect_identical(dt$cell_ID, fx$cells)
    expect_identical(dt$group, fx$grp)                      # -1 -> NA kept
    expect_identical(dt$label, paste0("L", seq_along(fx$cells)))
    expect_identical(dt$flag, seq_along(fx$cells) %% 2 == 0)
    expect_identical(is.na(dt$n_extra), fx$cnt_mask)
    expect_equal(dt$n_extra[!fx$cnt_mask], fx$cnt[!fx$cnt_mask])
})

test_that("var ids match the expression input, also for raw", {
    fx <- make_h5ad()
    fm <- h5adFeatMeta(fx$path)
    expect_identical(fm[]$feat_ID, h5adInput(fx$path)@feat_ids)
    expect_identical(fm[]$highly_variable, seq_along(fx$genes) %% 3 == 0)
    expect_identical(h5adFeatMeta(fx$path, layer = "raw")[]$feat_ID, fx$genes[1:10])
})

test_that("spatial locations keep orientation, with optional y flip and z", {
    fx <- make_h5ad()
    sl <- h5adSpatLocs(fx$path)[]
    expect_equal(sl$sdimx, fx$xy[, 1]); expect_equal(sl$sdimy, fx$xy[, 2])
    expect_identical(sl$cell_ID, fx$cells)
    expect_equal(h5adSpatLocs(fx$path, flip_y = TRUE)[]$sdimy, -fx$xy[, 2])
    expect_true("sdimz" %in% names(h5adSpatLocs(fx$path, key = "spatial3d")[]))
})

test_that("dim reductions carry PCA loadings and eigenvalues", {
    fx <- make_h5ad()
    dr <- h5adDimReds(fx$path, verbose = FALSE)
    expect_setequal(names(dr), c("pca", "umap"))
    expect_equal(unname(dr$pca[]), fx$pca)
    expect_identical(rownames(dr$pca[]), fx$cells)
    expect_identical(dr$pca@reduction_method, "pca")
    expect_equal(unname(dr$pca@misc$loadings), fx$load)
    expect_identical(rownames(dr$pca@misc$loadings),
                     .disambiguate_feat_ids(fx$genes))
    expect_equal(dr$pca@misc$eigenvalues, c(4, 3, 2, 1))
    expect_equal(unname(dr$umap[]), fx$umap)
    expect_identical(dr$umap@reduction_method, "umap")
})


# networks ####

test_that("the kNN streams without self-loops, weighted by connectivity", {
    fx <- make_h5ad()
    nets <- h5adNetworks(fx$path, gsrc, batch_records = 20, verbose = FALSE)
    expect_setequal(names(nets), c("kNN.pca", "spatial_neighbors"))
    knn <- nets[["kNN.pca"]]
    expect_s4_class(knn, "nnNetObj")
    st <- knn[]
    expect_true(st@directed)
    expect_identical(st@n_edges, as.numeric(length(fx$nb)))
    expect_gt(length(list.files(file.path(st@path, "edges"))), 1L)

    got <- data.table::as.data.table(storeRead(st, output = "tibble"))
    data.table::setnames(got, c("from_id", "to_id"), c("from", "to"))
    from <- rep(seq_along(fx$cells), each = ncol(fx$nb)); to <- as.vector(t(fx$nb))
    want <- data.table::data.table(from = fx$cells[from], to = fx$cells[to],
        weight = fx$conn[cbind(from, to)], distance = fx$D[cbind(from, to)])
    data.table::setorder(got, from, to); data.table::setorder(want, from, to)
    expect_equal(got[, .(from, to, weight, distance)], want)
    expect_true(all(got[, .N, by = from]$N == ncol(fx$nb)))  # uniform degree

    # same edges through the in-memory writer give the same store contents
    ref <- sourceWrite(gsrc, edgeDTInput(want, weight_col = "weight",
        distance_col = "distance"), store_type = "parquetEdge", type = "kNN")
    ref_dt <- data.table::as.data.table(storeRead(ref, output = "tibble"))
    data.table::setnames(ref_dt, c("from_id", "to_id"), c("from", "to"))
    data.table::setorder(ref_dt, from, to)
    expect_equal(got[, .(from, to, weight, distance)],
                 ref_dt[, .(from, to, weight, distance)])
})

test_that("a connectivities-only graph is undirected, upper triangle", {
    fx <- make_h5ad()
    sp <- h5adNetworks(fx$path, gsrc, keys = "spatial_neighbors",
                       verbose = FALSE)[[1L]]
    expect_s4_class(sp, "spatialNetworkObj")
    st <- sp[]
    expect_false(st@directed)
    up <- methods::as(Matrix::triu(fx$conn, 1), "TsparseMatrix")
    expect_identical(st@n_edges, as.numeric(length(up@x)))
    e <- data.table::setDT(dplyr::collect(arrow::open_dataset(file.path(st@path, "edges"))))
    expect_true(all(e$from_id < e$to_id))
    expect_equal(sort(e$weight), sort(up@x))
})


# reader ####

test_that("importH5adDisk() builds a backed object from every slot", {
    withr::local_options(giotto.use_conda = FALSE, giotto.no_python_warn = TRUE)
    fx <- make_h5ad(seed = 4L)
    rdr <- importH5adDisk(fx$path, backend = gsrc, layer = "counts")
    expect_s4_class(rdr, "H5adDiskReader")
    expect_identical(rdr$layer, "counts")
    g <- rdr$create_gobject()

    ex <- GiottoClass::getExpression(g, values = "raw", output = "exprObj")
    expect_equal(.store_triplets(ex[]), .source_triplets(fx$M * 2))
    expect_identical(GiottoClass::spatIDs(g), rownames(fx$M))
    expect_identical(nrow(GiottoClass::getSpatialLocations(g, output = "data.table")),
                     nrow(fx$M))
    expect_true("group" %in% names(GiottoClass::pDataDT(g)))
    expect_true("highly_variable" %in% names(GiottoClass::fDataDT(g)))
    dr <- GiottoClass::list_dim_reductions(g)$name
    expect_true(all(c("pca", "umap") %in% dr))
    expect_true("kNN.pca" %in% GiottoClass::list_nearest_networks(g)$name)
    expect_true("spatial_neighbors" %in% GiottoClass::list_spatial_networks(g)$name)
})

test_that("create_gobject() skips what it is told to and adds into a gobject", {
    withr::local_options(giotto.use_conda = FALSE, giotto.no_python_warn = TRUE)
    fx <- make_h5ad(seed = 5L)
    rdr <- importH5adDisk(fx$path, backend = gsrc)
    g <- rdr$create_gobject(load_dimreds = FALSE, load_networks = FALSE)
    expect_null(GiottoClass::list_dim_reductions(g))
    expect_null(GiottoClass::list_nearest_networks(g))
    g <- rdr$create_gobject(load_expression = FALSE, load_spatlocs = FALSE,
        load_cellmeta = FALSE, load_featmeta = FALSE, load_networks = FALSE,
        gobject = g)
    expect_true(all(c("pca", "umap") %in% GiottoClass::list_dim_reductions(g)$name))
})


test_that("a dense X streams, and the sparse layers are named", {
    fx <- make_h5ad(seed = 6L)
    h5 <- hdf5r::H5File$new(fx$path, mode = "r+")
    h5$link_delete("X")
    .w_array(h5, "X", as.matrix(fx$M))
    h5$close_all()
    expect_message(inp <- h5adInput(fx$path, batch_records = 100), "'counts'")
    pe <- sourceWrite(gsrc, inp, store_type = "parquetExpr")
    expect_equal(.store_triplets(pe), .source_triplets(fx$M))
})


# spatial network names ####

.with_params <- function(path, params) {
    h5 <- hdf5r::H5File$new(path, mode = "r+")
    on.exit(h5$close_all())
    pg <- h5[["uns/spatial_neighbors"]]$create_group("params")
    for (n in names(params)) pg[[n]] <- params[[n]]
    invisible(path)
}
.spat_net_name <- function(path) {
    names(h5adNetworks(path, gsrc, keys = "spatial_neighbors", verbose = FALSE))
}

test_that("squidpy params name the spatial network the Giotto way", {
    expect_identical(.spat_net_name(.with_params(make_h5ad()$path,
        list(coord_type = "generic", n_neighbors = 6L, transform = "none"))),
        "kNN_network")
    expect_identical(.spat_net_name(.with_params(make_h5ad()$path,
        list(coord_type = "generic", transform = "none"))),
        "Delaunay_network")
    expect_identical(.spat_net_name(.with_params(make_h5ad()$path,
        list(coord_type = "grid", n_neighbors = 6L, n_rings = 1L,
             delaunay = FALSE))),
        "spatial_neighbors")                                # grid: no Giotto kind
    expect_identical(.spat_net_name(.with_params(make_h5ad()$path,
        list(coord_type = "generic", radius = 50))),
        "spatial_neighbors")                                # radius graph
    expect_identical(.spat_net_name(.with_params(make_h5ad()$path,
        list(n_neighbors = 6L, coord_type = "generic", radius = 50,
             transform = "none"))),
        "spatial_neighbors")                                # older squidpy
    expect_identical(.spat_net_name(make_h5ad()$path), "spatial_neighbors")
})
