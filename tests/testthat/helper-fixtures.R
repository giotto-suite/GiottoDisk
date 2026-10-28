# Synthetic fixtures shared across test files. File-specific shapes stay in
# their test file as thin wrappers that pin their own defaults.

# Sparse count matrix: Poisson(lambda) + 1 nonzeros, so every stored value is
# >= 1 and no explicit zeros reach the store. `pad` zero-pads the numeric part
# of the gene / cell names (recycled to length 2); 0 means unpadded.
.rand_counts <- function(n_genes = 12L, n_cells = 30L, density = 0.5,
                         seed = 1L, lambda = 5L,
                         gene_prefix = "g", cell_prefix = "c", pad = 0L) {
    pad <- rep_len(pad, 2L)
    ids <- function(prefix, n, width) {
        i <- seq_len(n)
        if (width > 0L) i <- formatC(i, width = width, flag = "0")
        paste0(prefix, i)
    }
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = density,
        rand.x = function(n) as.double(rpois(n, lambda) + 1L))
    rownames(m) <- ids(gene_prefix, n_genes, pad[[1L]])
    colnames(m) <- ids(cell_prefix, n_cells, pad[[2L]])
    m
}

# n points on the diagonal (i, i), keyed cell_1..cell_n.
make_pts <- function(n = 5) {
    terra::vect(
        data.frame(
            x = seq_len(n),
            y = seq_len(n),
            cell_ID = paste0("cell_", seq_len(n))
        ),
        geom = c("x", "y"),
        crs = ""
    )
}

# 5x5 grid of points at tile centres (5, 15, ..., 45) over a 50x50 extent.
make_tile_pts <- function() {
    coords <- data.frame(
        x = rep(seq(5, 45, by = 10), 5),
        y = rep(seq(5, 45, by = 10), each = 5)
    )
    terra::vect(
        data.frame(coords, id = seq_len(nrow(coords))),
        geom = c("x", "y"),
        crs = ""
    )
}
