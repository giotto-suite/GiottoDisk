# A small CosMx flat-file directory: FOV 1 holds a grid of `n_fov1` square
# cells, FOV 2 holds four. It carries the export quirks the disk reader has
# to absorb:
#   - one vertex written as `1e+05`, as AtoMx writes round coordinates. It
#     sits in FOV 2, after > 1 MB of integer vertices, because arrow infers
#     column types from the first block of the file and only an int64 guess
#     fails on it
#   - one polygon (c_1_2_4) for a cell with no expression or metadata row,
#     as for the cells the vendor's QC removes
#   - one cell (c_1_2_3) with expression and metadata but no polygon
# The polygon file is written as text so the `1e+05` reaches the reader
# verbatim.

.mk_cosmx_disk_dir <- function(n_fov1 = 15000L) {
    d <- tempfile("cosmx_")
    dir.create(d)
    origin <- data.frame(FOV = 1:2, x_global_px = c(0L, 5000L),
                         y_global_px = c(50000L, 99000L))

    # cell centers, FOV-local; FOV 2 cell 1 has its top edge at y = 100000
    grid <- expand.grid(i = 0:149, j = 0:(ceiling(n_fov1 / 150) - 1L))
    grid <- grid[seq_len(n_fov1), ]
    cells <- data.frame(
        fov    = c(rep(1L, n_fov1), rep(2L, 4L)),
        cellID = c(seq_len(n_fov1), 1:4),
        lx     = c(30L * grid$i + 20L, c(100L, 300L, 500L, 700L)),
        ly     = c(30L * grid$j + 20L, c(990L, 100L, 100L, 100L))
    )
    cells$cell <- sprintf("c_1_%d_%d", cells$fov, cells$cellID)
    cells$gx <- origin$x_global_px[cells$fov] + cells$lx
    cells$gy <- origin$y_global_px[cells$fov] + cells$ly

    # polygons: every cell except c_1_2_3; c_1_2_4 is the orphan
    r <- 10L
    has_poly <- cells$cell != "c_1_2_3"
    pc <- cells[rep(which(has_poly), each = 4L), ]
    dx <- rep(c(-r, r, r, -r), sum(has_poly))
    dy <- rep(c(-r, -r, r, r), sum(has_poly))
    lines <- sprintf("%d,%d,%s,%d,%d,%d,%s",
                     pc$fov, pc$cellID, pc$cell, pc$lx + dx, pc$ly + dy,
                     pc$gx + dx, format(pc$gy + dy, scientific = FALSE))
    # write.csv-style formatting of round numbers: 100000 -> "1e+05"
    lines <- sub(",100000$", ",1e+05", lines)
    hdr <- "fov,cellID,cell,x_local_px,y_local_px,x_global_px,y_global_px"
    writeLines(c(hdr, lines),
               file.path(d, "Fixture-polygons.csv"))

    # expression + metadata: the cells the vendor kept
    kept <- cells[cells$cell != "c_1_2_4", ]
    n <- nrow(kept)
    counts <- cbind(GENEA = (seq_len(n) %% 5L) + 1L,
                    GENEB = seq_len(n) %% 3L,
                    Negative1 = seq_len(n) %% 2L,
                    SystemControl1 = as.integer(seq_len(n) %% 7L == 0L))
    data.table::fwrite(data.frame(fov = kept$fov, cell_ID = kept$cellID,
                                  counts),
                       file.path(d, "Fixture_exprMat_file.csv"))
    data.table::fwrite(data.frame(
        fov = kept$fov, cell_ID = kept$cellID, cell = kept$cell,
        CenterX_global_px = kept$gx, CenterY_global_px = kept$gy,
        Area.um2 = 5.8, nCount_RNA = rowSums(counts[, 1:2]),
        nFeature_RNA = rowSums(counts[, 1:2] > 0)),
        file.path(d, "Fixture_metadata_file.csv"))
    data.table::fwrite(cbind(origin, x_global_mm = origin$x_global_px * 1.2e-4,
                             y_global_mm = origin$y_global_px * 1.2e-4),
                       file.path(d, "Fixture_fov_positions_file.csv"))
    d
}
