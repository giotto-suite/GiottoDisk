#' @include pkg_imports.R class-parquetExprStore.R class-parquetEdgeStore.R
NULL

# processData(parquetExprBase, spatialSmoothParam): each cell's values
# averaged with its neighbours', S = X M, written to a new store (GiottoClass
# defines the param and the in-memory method). Composed from existing blocks:
#
#   pairs in view       .edge_to_cell_map + .edge_view_pairs
#   operator M          GiottoClass::spatial_smooth_operator
#   local cell order    .pairs_local, else .pairs_bfs_order +
#                       .pe_rewrite_cell_order (which also bakes the chain)
#   windowed product    per window of target cells, the row range their
#                       neighbourhoods span, read with .pe_read_chunk_sub
#   write               .ingest_spill_iterator + .pestore_write_itr
#
# The result keeps the input's cell and feature order. It is denser than the
# input -- a cell gains every feature its neighbours have -- so smooth a
# feature subset (`x[feats, ]`) rather than a whole panel.

#' @rdname processData
#' @param network for a `spatialSmoothParam`: the spatial network as a
#'   `parquetEdgeStore`.
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "spatialSmoothParam"),
    function(x, param, ..., network) {
        if (!inherits(network, "parquetEdgeStore")) {
            stop("[processData(parquetExprBase, spatialSmoothParam)] ",
                "`network` must be a parquetEdgeStore.", call. = FALSE)
        }
        n_cells <- as.integer(x@n_cells)
        pairs <- .edge_view_pairs(network, .edge_to_cell_map(network, x),
            isTRUE(param$weighted))
        M <- spatial_smooth_operator(pairs$a, pairs$b, pairs$w, n = n_cells,
            smoothing = param$smoothing)
        M1 <- M
        for (h in seq_len(param$hops - 1L)) M <- M %*% M1

        # a window's product is several times its input: the source range,
        # the operator block and the denser output
        budget <- max(1L, as.integer(.exprbase_chunk_size(x) %/% 4L))
        ord <- if (budget >= n_cells || .pairs_local(pairs, budget)) {
            seq_len(n_cells)
        } else {
            .pairs_bfs_order(pairs, n_cells)
        }
        tmp <- .dump_tempfile()
        on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
        src <- if (identical(ord, seq_len(n_cells))) {
            storeWrite(parquetExprStore(path = tmp), x)
        } else {
            .pe_rewrite_cell_order(x, ord, tmp)
        }
        .pe_smooth_write(src, M[ord, ord], ord, x@cell_ids, x@feat_ids,
            budget)
    }
)

# S = X M for a materialized single store `src` whose cells are the view's
# reordered by `ord`; `Mo` is M in that order. Triplets are written under the
# cells' view positions, so the new store has the view's order.
.pe_smooth_write <- function(src, Mo, ord, cell_ids, feat_ids, budget) {
    n_cells <- length(cell_ids)
    info <- list(sub = .exprbase_substores(src)[[1L]]$store)
    windows <- .pe_chunk_ranges(1L, n_cells, budget)
    fpos <- stats::setNames(seq_along(feat_ids), feat_ids)
    Mo <- methods::as(Mo, "CsparseMatrix")

    max_rows <- .ingest_batch_rows(96)
    dens <- tryCatch(.pe_view_density(src), error = function(e) 1)
    n_rec <- as.numeric(n_cells) * length(feat_ids) * min(1, 4 * dens)
    n_buckets <- .ingest_n_buckets(n_rec, max_rows)
    per_bucket <- max(1L, as.integer(ceiling(n_cells / n_buckets)))

    produce <- function(i) {
        tgt <- windows[[i]][[1L]]:windows[[i]][[2L]]
        blk <- Mo[, tgt, drop = FALSE]
        rng <- range(blk@i + 1L)
        X <- .pe_read_chunk_sub(info, rng[1L], rng[2L], list(),
            length(feat_ids))
        if (is.null(X)) return(NULL)
        # columns by cell ID onto the full range: a cell without values in
        # the range has no column in the read
        Xt <- methods::as(X, "TsparseMatrix")
        cj <- match(colnames(X), src@cell_ids[rng[1L]:rng[2L]])
        X <- Matrix::sparseMatrix(i = Xt@i + 1L, j = cj[Xt@j + 1L],
            x = Xt@x, dims = c(nrow(X), rng[2L] - rng[1L] + 1L),
            dimnames = list(rownames(X), NULL))
        S <- methods::as(X %*% blk[rng[1L]:rng[2L], , drop = FALSE],
            "TsparseMatrix")
        if (!length(S@x)) return(NULL)
        vp <- ord[tgt[S@j + 1L]]
        data.table::data.table(row_id = vp,
            col_id = unname(fpos[rownames(X)[S@i + 1L]]),
            value = S@x, b = (vp - 1L) %/% per_bucket + 1L)
    }
    finish <- function(dt) {
        dt <- dt[, list(row_id = as.integer(row_id), col_id = as.integer(col_id),
            value = as.numeric(value))]
        data.table::setkeyv(dt, c("row_id", "col_id"))
        dt
    }
    row_id <- col_id <- value <- NULL  # NSE
    itr <- .ingest_spill_iterator(produce, length(windows), n_buckets,
        finish, max_rows)
    itr$cell_ids <- function() cell_ids
    itr$feat_ids <- function() feat_ids
    itr$n_cells <- function() n_cells
    itr$n_genes <- function() length(feat_ids)

    path <- .dump_tempfile()
    out <- parquetExprStore(path = path)
    partition_dir <- .idpath(path, out@uid)
    dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)
    .pestore_write_itr(out, itr, partition_dir)
}
