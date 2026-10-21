# A small MERSCOPE region directory: 6 square cells over 2 FOVs (cells 1-3 in
# FOV 0, 4-6 in FOV 1), genes A, B and the control Blank-1, three transcripts
# per cell. Cell IDs are 19-digit, as in real exports, so anything that reads
# them as doubles loses digits and fails the comparisons.
#
# Mirrors Giotto's tests/testthat/helper-insitu-fixture.R.

# built as text: a double cannot hold 19 digits
.mer_ids <- function() sprintf("24510672000001000%02d", 1:6)

.mk_merscope_disk_dir <- function(tx_format = c("csv", "parquet"),
    bounds = c("parquet", "hdf5"), z_planes = 1L) {
    tx_format <- match.arg(tx_format)
    bounds <- match.arg(bounds)
    d <- tempfile("merscope_")
    dir.create(d)
    ids <- .mer_ids()
    fov <- rep(0:1, each = 3L)
    cx <- 10 * seq_along(ids)
    cy <- rep(c(10, 30), length.out = length(ids))
    sq <- function(x, y, r = 2) {
        cbind(c(x - r, x + r, x + r, x - r, x - r),
              c(y - r, y - r, y + r, y + r, y - r))
    }

    if (bounds == "parquet") {
        wkb <- lapply(seq_along(ids), function(i) {
            terra::geom(terra::vect(list(sq(cx[i], cy[i])), type = "polygons"),
                wkb = TRUE)[[1L]]
        })
        n <- length(ids) * z_planes
        tb <- arrow::arrow_table(
            EntityID = arrow::Array$create(rep(ids, z_planes))$cast(arrow::int64()),
            ZIndex = rep(seq_len(z_planes) - 1L, each = length(ids)),
            Type = rep("cell", n)
        )
        tb$Geometry <- arrow::Array$create(rep(wkb, z_planes),
            type = arrow::binary())
        arrow::write_parquet(tb, file.path(d, "cell_boundaries.parquet"))
    } else {
        # Vizgen per-FOV HDF5: /featuredata/<cell>/zIndex_<k>/p_0/coordinates
        # holds a 2 x n x 1 array; z_coordinates lists the 7 plane depths
        bdir <- file.path(d, "cell_boundaries")
        dir.create(bdir)
        for (f in unique(fov)) {
            h5 <- file.path(bdir, sprintf("feature_data_%d.hdf5", f))
            rhdf5::h5createFile(h5)
            rhdf5::h5createGroup(h5, "featuredata")
            for (i in which(fov == f)) {
                cell <- file.path("featuredata", ids[i])
                rhdf5::h5createGroup(h5, cell)
                rhdf5::h5write(seq(1.5, 10.5, by = 1.5),
                    h5, file.path(cell, "z_coordinates"))
                for (z in 0:6) {
                    zg <- file.path(cell, sprintf("zIndex_%d", z))
                    rhdf5::h5createGroup(h5, zg)
                    rhdf5::h5createGroup(h5, file.path(zg, "p_0"))
                    xy <- t(sq(cx[i], cy[i]))
                    rhdf5::h5write(array(xy, dim = c(2L, ncol(xy), 1L)),
                        h5, file.path(zg, "p_0", "coordinates"))
                }
            }
            rhdf5::h5closeAll()
        }
    }

    genes <- c("A", "B", "Blank-1")
    tx <- data.table::data.table(
        V1 = seq_len(3L * length(ids)) - 1L,
        barcode_id = 1L,
        global_x = rep(cx + 0.5, each = 3L),
        global_y = rep(cy - 0.5, each = 3L),
        global_z = 0,
        x = rep(cx, each = 3L) * 9, y = rep(cy, each = 3L) * 9,
        fov = rep(fov, each = 3L),
        gene = rep(genes, times = length(ids)),
        transcript_id = "t",
        cell_id = rep(ids, each = 3L)
    )
    data.table::setnames(tx, "V1", "")
    if (tx_format == "csv") {
        data.table::fwrite(tx, file.path(d, "detected_transcripts.csv"))
    } else {
        arrow::write_parquet(tx, file.path(d, "detected_transcripts.parquet"))
    }

    data.table::fwrite(data.table::data.table(
        cell = ids, A = seq_along(ids), B = 2L, `Blank-1` = 1L
    ), file.path(d, "cell_by_gene.csv"))
    data.table::fwrite(data.table::data.table(
        EntityID = ids, fov = fov, volume = 10, center_x = cx, center_y = cy
    ), file.path(d, "cell_metadata.csv"))
    d
}
