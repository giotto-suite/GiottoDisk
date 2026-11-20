# A synthetic AnnData .h5ad, written with hdf5r in anndata's on-disk encoding
# (encoding-type / encoding-version / shape / _index / column-order
# attributes), small enough to write in every test run. It carries:
#   X (CSR or CSC), layers/counts (CSR), layers/dense (refused), raw/X + raw/var
#   obs: numeric, categorical with a missing code, string, boolean, and a
#        nullable integer; var with a duplicate name
#   obsm: spatial, spatial3d, X_pca, X_umap; varm/PCs; uns/pca/variance
#   obsp: a k = 5 kNN whose distances keep each cell as its own neighbour at
#         distance 0, umap-style symmetric connectivities with one kNN pair
#         dropped, and a connectivities-only spatial graph; uns names both
# make_h5ad() returns the path and the source objects every read is checked
# against.

.w_attr <- function(node, ...) {
    a <- list(...)
    for (n in names(a)) hdf5r::h5attr(node, n) <- a[[n]]
}
.w_sparse <- function(h5, p, M, format = c("csr", "csc")) {
    format <- match.arg(format)
    g <- h5$create_group(p)
    S <- if (format == "csr") as(t(M), "CsparseMatrix") else as(M, "CsparseMatrix")
    g[["data"]] <- S@x; g[["indices"]] <- S@i; g[["indptr"]] <- S@p
    .w_attr(g, `encoding-type` = paste0(format, "_matrix"),
            `encoding-version` = "0.1.0", shape = dim(M))
}
.w_array <- function(h5, p, m) {
    h5[[p]] <- if (is.matrix(m)) t(m) else m   # stored (n_obs, k)
    .w_attr(h5[[p]], `encoding-type` = "array", `encoding-version` = "0.2.0")
}
.w_df <- function(h5, p, index, cols) {
    g <- h5$create_group(p)
    g[["_index"]] <- index
    for (n in names(cols)) {
        v <- cols[[n]]
        if (is.factor(v)) {
            cg <- g$create_group(n)
            codes <- as.integer(v) - 1L; codes[is.na(codes)] <- -1L
            cg[["codes"]] <- codes; cg[["categories"]] <- levels(v)
            .w_attr(cg, `encoding-type` = "categorical",
                    `encoding-version` = "0.2.0", ordered = FALSE)
        } else if (is.list(v)) {          # list(values, mask): nullable
            ng <- g$create_group(n)
            ng[["values"]] <- v$values; ng[["mask"]] <- v$mask
            .w_attr(ng, `encoding-type` = v$enc, `encoding-version` = "0.1.0")
        } else {
            g[[n]] <- v
            .w_attr(g[[n]], `encoding-type` =
                if (is.character(v)) "string-array" else "array",
                `encoding-version` = "0.2.0")
        }
    }
    .w_attr(g, `encoding-type` = "dataframe", `encoding-version` = "0.2.0",
            `_index` = "_index", `column-order` = names(cols))
}

make_h5ad <- function(n_obs = 60L, n_var = 25L, format = "csr", seed = 1L) {
    set.seed(seed)
    M <- Matrix::rsparsematrix(n_obs, n_var, density = 0.2,
                       rand.x = function(n) as.double(rpois(n, 4) + 1))
    cells <- sprintf("cell_%03d", seq_len(n_obs))
    genes <- sprintf("g%02d", seq_len(n_var))
    genes[n_var] <- genes[1L]                         # a duplicate var name
    dimnames(M) <- list(cells, genes)

    # kNN with the cell itself first at distance 0, as some writers store it
    k <- 5L
    xy <- matrix(runif(n_obs * 2) * 100, ncol = 2)
    D <- as.matrix(dist(xy))
    nb <- t(apply(D, 1, function(d) order(d)[seq_len(k + 1L)]))   # incl. self
    dist_m <- Matrix::sparseMatrix(i = rep(seq_len(n_obs), each = k + 1L),
        j = as.vector(t(nb)), x = D[cbind(rep(seq_len(n_obs), each = k + 1L),
        as.vector(t(nb)))], dims = c(n_obs, n_obs))
    nb_off <- nb[, -1L]
    A <- Matrix::sparseMatrix(i = rep(seq_len(n_obs), each = k), j = as.vector(t(nb_off)),
        x = 1 / (1 + D[cbind(rep(seq_len(n_obs), each = k), as.vector(t(nb_off)))]),
        dims = c(n_obs, n_obs))
    conn <- pmax(A, t(A))                              # symmetric, like umap
    # umap's fuzzy graph can drop a kNN pair; that edge must keep weight 0
    p1 <- c(1L, nb_off[1L, 1L])
    conn[p1[1], p1[2]] <- 0; conn[p1[2], p1[1]] <- 0; conn <- Matrix::drop0(conn)

    path <- tempfile(fileext = ".h5ad")
    h5 <- hdf5r::H5File$new(path, mode = "w")
    .w_attr(h5, `encoding-type` = "anndata", `encoding-version` = "0.1.0")
    .w_sparse(h5, "X", M, format)
    lg <- h5$create_group("layers")
    .w_sparse(lg, "counts", M * 2, "csr")
    .w_array(lg, "dense", as.matrix(M))
    rg <- h5$create_group("raw")
    .w_sparse(rg, "X", M[, 1:10], "csr")
    .w_df(rg, "var", genes[1:10], list(gene_len = as.numeric(1:10)))

    grp <- factor(sample(c("a", "b", NA), n_obs, TRUE), levels = c("a", "b"))
    cnt <- rpois(n_obs, 10); cnt_mask <- seq_len(n_obs) %% 7 == 0
    .w_df(h5, "obs", cells, list(
        n_counts = as.numeric(Matrix::rowSums(M)),
        group    = grp,
        label    = paste0("L", seq_len(n_obs)),
        flag     = seq_len(n_obs) %% 2 == 0,
        n_extra  = list(values = cnt, mask = cnt_mask, enc = "nullable-integer")))
    .w_df(h5, "var", genes, list(highly_variable = seq_len(n_var) %% 3 == 0))

    om <- h5$create_group("obsm")
    pca <- matrix(rnorm(n_obs * 4), ncol = 4)
    umap <- matrix(rnorm(n_obs * 2), ncol = 2)
    .w_array(om, "spatial", xy)
    .w_array(om, "spatial3d", cbind(xy, 1))
    .w_array(om, "X_pca", pca)
    .w_array(om, "X_umap", umap)
    vm <- h5$create_group("varm")
    load <- matrix(rnorm(n_var * 4), ncol = 4)
    .w_array(vm, "PCs", load)

    op <- h5$create_group("obsp")
    .w_sparse(op, "distances", dist_m, "csr")
    .w_sparse(op, "connectivities", conn, "csr")
    .w_sparse(op, "spatial_connectivities", conn, "csr")
    un <- h5$create_group("uns")
    nbg <- un$create_group("neighbors")
    nbg[["connectivities_key"]] <- "connectivities"
    nbg[["distances_key"]] <- "distances"
    pg <- un$create_group("pca")
    pg[["variance"]] <- c(4, 3, 2, 1)
    sp <- un$create_group("spatial_neighbors")
    sp[["connectivities_key"]] <- "spatial_connectivities"
    h5$close_all()

    list(path = path, M = M, cells = cells, genes = genes, grp = grp,
         cnt = cnt, cnt_mask = cnt_mask, xy = xy, pca = pca, umap = umap,
         load = load, D = D, nb = nb_off, conn = conn)
}

# comparison helpers ####

.store_triplets <- function(pe) {
    dt <- data.table::setDT(dplyr::collect(arrow::open_dataset(pe@path)))
    dt <- dt[, .(row_id = as.integer(row_id), col_id = as.integer(col_id), value)]
    data.table::setorder(dt, row_id, col_id)[]
}
.source_triplets <- function(M) {
    s <- methods::as(M, "TsparseMatrix")
    data.table::setorder(data.table::data.table(
        row_id = s@i + 1L, col_id = s@j + 1L, value = s@x),
        row_id, col_id)[]
}
.parts_cell_major <- function(pe) {
    fs <- sort(list.files(pe@path, "[.]parquet$", recursive = TRUE, full.names = TRUE))
    fs <- fs[order(as.integer(sub(".*part-(\\d+)[.]parquet$", "\\1", fs)))]
    r <- t(vapply(fs, function(f) {
        x <- arrow::read_parquet(f, col_select = "row_id")$row_id
        c(min(x), max(x), !is.unsorted(x))
    }, numeric(3)))
    length(fs) > 1L && all(r[, 3] == 1) && all(r[-1, 1] > r[-nrow(r), 2])
}
