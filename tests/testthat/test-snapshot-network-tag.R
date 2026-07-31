# Tests for snapshotSave's network adoption + uid-tagging path.
#
# Backed gobjects with parquetEdgeStore networks (from setter auto-write or
# direct sourceWrite) must have their network artifacts tagged with the
# snapshot name in the source manifest. Otherwise sourcePrune sees the
# manifest entry untagged by any snapshot and deletes the artifact, leaving
# the snapshot's serialized handle dangling.

skip_if_not_installed("GiottoClass")
skip_if_not_installed("GiottoData")


.mk_backed_g <- function(td) {
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    src <- gDirSource(td)
    g@source <- src
    g
}

.attach_pes_to_nn_network <- function(g, slot_name = "sNN.pca", type = "sNN") {
    nn <- GiottoClass::getNearestNetwork(g, name = slot_name,
        output = "nnNetObj")
    stopifnot(inherits(nn@network, "igraph"))
    pes <- sourceWrite(g@source, nn@network, type = type)
    nn@network <- pes
    g@nn_network[["cell"]][["rna"]][[slot_name]] <- nn
    list(g = g, pes = pes)
}


# ---- detect path: network already vault-resident gets tagged --------------

test_that("snapshotSave tags parquetEdgeStore (nn_network) with snapshot name", {
    td <- tempfile("netTag_nn_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    res <- .attach_pes_to_nn_network(g)
    g <- res$g
    pes_uid <- res$pes@uid

    snapshotSave(g@source, g, name = "snap_v1", verbose = FALSE)

    # Read consolidated manifest from disk
    mfst <- GiottoDisk:::.gdsrc_json_read(g@source@path,
        consolidate = TRUE)$content
    entry <- mfst[[pes_uid]]
    expect_false(is.null(entry))
    expect_identical(entry$store, "parquetEdgeStore")
    expect_identical(entry$giottosave, "snap_v1")
})

test_that("snapshotSave tags parquetEdgeStore (spatial_network) with snapshot name", {
    td <- tempfile("netTag_sn_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    sn <- GiottoClass::getSpatialNetwork(g, name = "Delaunay_network",
        output = "spatialNetworkObj")
    skip_if_not(inherits(sn@network, "igraph"),
        "Spatial network on mini fixture is not igraph-backed")
    pes <- sourceWrite(g@source, sn@network, type = "spatial")
    sn@network <- pes
    g@spatial_network[["cell"]][["Delaunay_network"]] <- sn

    snapshotSave(g@source, g, name = "snap_v1", verbose = FALSE)

    mfst <- GiottoDisk:::.gdsrc_json_read(g@source@path,
        consolidate = TRUE)$content
    entry <- mfst[[pes@uid]]
    expect_false(is.null(entry))
    expect_identical(entry$store, "parquetEdgeStore")
    expect_identical(entry$giottosave, "snap_v1")
})

test_that("multiple snapshots accumulate in the network manifest entry", {
    td <- tempfile("netTag_multi_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    res <- .attach_pes_to_nn_network(g)
    g <- res$g
    pes_uid <- res$pes@uid

    snapshotSave(g@source, g, name = "snap_v1", verbose = FALSE)
    snapshotSave(g@source, g, name = "snap_v2", verbose = FALSE)

    mfst <- GiottoDisk:::.gdsrc_json_read(g@source@path,
        consolidate = TRUE)$content
    entry <- mfst[[pes_uid]]
    expect_setequal(entry$giottosave, c("snap_v1", "snap_v2"))
})


# ---- detect_uid wiring ----------------------------------------------------

test_that(".ss_gdsrc_detect_uid surfaces network uids alongside other types", {
    td <- tempfile("netTag_detect_"); dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)

    g <- .mk_backed_g(td)
    res <- .attach_pes_to_nn_network(g)
    g <- res$g

    uids <- GiottoDisk:::.ss_gdsrc_detect_uid(g)
    expect_true(res$pes@uid %in% uids)
})
