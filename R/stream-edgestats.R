#' @include pkg_imports.R class-parquetExprStore.R class-parquetEdgeStore.R
NULL

# analyzeData(parquetExprBase, edgeStatsParam): the per-feature sums over a
# spatial network that edge statistics are computed from (GiottoClass defines
# the contract and the in-memory method; Giotto computes Moran's I, Geary's C
# and join counts from it). Composed from existing blocks:
#
#   pairs in view      .edge_to_cell_map + .edge_view_pairs (edges read once)
#   degree, weights    .pairs_degree and sums over the pairs
#   power and degree-  .pe_accum_raw: one pass grouped by degree class, or a
#   weighted sums      multiply pass per weighted sum when degrees are many
#   edge products      .pe_edge_products
#
# Only edges with both ends among the view's cells count, as in the contract.
#
# With `groups =`, the "grouped" attribute adds per-group power sums (the same
# accumulator pass grouped by cell group), neighbour-group sums (a weighted
# `by_cell` pass, each cell listed under the groups its edges reach) and the
# group edge tables (.edge_group_tables), all from the same pairs.

#' @rdname analyzeData
#' @param network for an `edgeStatsParam`: the spatial network as a
#'   `parquetEdgeStore`.
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "edgeStatsParam"),
    function(x, param, ..., network, groups = NULL) {
        if (!inherits(network, "parquetEdgeStore")) {
            stop("[analyzeData(parquetExprBase, edgeStatsParam)] `network` ",
                "must be a parquetEdgeStore.", call. = FALSE)
        }
        weighted <- isTRUE(param$weighted)
        map <- .edge_to_cell_map(network, x)
        pairs <- .edge_view_pairs(network, map, weighted)

        # per view cell degree, 0 where a cell has no edge in view
        n_cells <- as.integer(x@n_cells)
        deg <- .pairs_degree(pairs, n_cells)

        sums <- .pe_degree_sums(x, deg)
        out <- data.table::data.table(
            feats = x@feat_ids,
            sum = sums$sum, sumsq = sums$sumsq,
            sum3 = sums$sum3, sum4 = sums$sum4,
            deg_sum = sums$deg_sum, deg_sumsq = sums$deg_sumsq,
            edge_prod = .pe_edge_products(x, network, pairs = pairs)
        )
        data.table::setattr(out, "graph", list(
            n_cells = n_cells, n_edges = nrow(pairs),
            sum_w = sum(pairs$w), sum_w2 = sum(pairs$w^2), sum_deg2 = sum(deg^2)
        ))
        if (!is.null(groups)) {
            data.table::setattr(out, "grouped",
                .pe_edge_stats_grouped(x, pairs, groups))
        }
        out
    }
)

# The "grouped" attribute of the contract on a store.
.pe_edge_stats_grouped <- function(x, pairs, groups) {
    pos <- k <- NULL  # NSE
    g <- .pe_cell_groups(x, groups)
    if (anyNA(g)) {
        stop("[edgeStatsParam] every cell needs a group; ", sum(is.na(g)),
            " have none.", call. = FALSE)
    }
    g <- droplevels(if (is.factor(g)) g else factor(g))
    lv <- levels(g)
    code <- as.integer(g)
    K <- length(lv)
    G <- as.integer(x@n_genes)

    cells <- data.table::copy(.pe_axis_pos_map(x, "cell"))
    data.table::setnames(cells, "key_id", "row_id")
    key <- c(intersect("source_id", names(cells)), "row_id")
    to_mat <- function(agg, col) {
        out <- matrix(0, G, K, dimnames = list(x@feat_ids, lv))
        if (!is.null(agg)) out[cbind(agg$pos, agg$k)] <- agg[[col]]
        out
    }

    pw <- c("sum", "sumsq", "sum3", "sum4")
    by_group <- cells[, c(key, "pos"), with = FALSE][, k := code[pos]][
        , c(key, "k"), with = FALSE]
    agg <- .pe_accum_raw(x, axis = "feat", stats = pw, by_cell = by_group)

    gt <- .edge_group_tables(pairs, code, K)
    nb <- merge(gt$a, cells[, c(key, "pos"), with = FALSE], by = "pos")
    nb_cell <- nb[, c(key, "l", "a"), with = FALSE]
    data.table::setnames(nb_cell, c("l", "a"), c("k", "w"))
    nbr <- .pe_accum_raw(x, axis = "feat", stats = "wsum", by_cell = nb_cell)

    dimnames(gt$edges) <- list(lv, lv)
    dimnames(gt$paths) <- list(lv, lv, lv)
    list(
        levels = lv,
        n_cells = stats::setNames(as.numeric(tabulate(code, K)), lv),
        sum = to_mat(agg, "sum"), sumsq = to_mat(agg, "sumsq"),
        sum3 = to_mat(agg, "sum3"), sum4 = to_mat(agg, "sum4"),
        nbr_sum = to_mat(nbr, "wsum"),
        edges = gt$edges,
        paths = gt$paths
    )
}

# Per-feature power sums and degree-weighted sums, `deg` given per view cell.
#
# Few distinct degrees (any unweighted network): one pass, grouped by degree
# class through `by_cell`, and the weighted sums are the class sums times the
# class degree. Otherwise one ungrouped pass for the power sums and a
# `multiply(axis = "cell")` pass for each weighted sum.
.pe_degree_sums <- function(x, deg, max_classes = 256L) {
    k <- pos <- NULL  # NSE
    G <- as.integer(x@n_genes)
    pw <- c("sum", "sumsq", "sum3", "sum4")
    classes <- sort(unique(deg))
    if (length(classes) <= max_classes) {
        map <- data.table::copy(.pe_axis_pos_map(x, "cell"))
        map[, k := match(deg[pos], classes)]
        data.table::setnames(map, "key_id", "row_id")
        by_cell <- map[, c(intersect("source_id", names(map)), "row_id", "k"),
            with = FALSE]
        agg <- .pe_accum_raw(x, axis = "feat", stats = pw, by_cell = by_cell)
        out <- lapply(c(pw, "deg_sum", "deg_sumsq"), function(i) numeric(G))
        names(out) <- c(pw, "deg_sum", "deg_sumsq")
        if (is.null(agg)) return(out)
        d <- classes[agg$k]
        data.table::set(agg, j = "deg_sum", value = d * agg$sum)
        data.table::set(agg, j = "deg_sumsq", value = d * agg$sumsq)
        cols <- c(pw, "deg_sum", "deg_sumsq")
        tot <- agg[, lapply(.SD, sum), by = pos, .SDcols = cols]
        for (nm in cols) out[[nm]][tot$pos] <- tot[[nm]]
        return(out)
    }
    out <- .stream_expr_accum(x, axis = "feat", stats = pw)
    by_deg <- function(f) .pe_push_op_last(x, list(type = "multiply",
        axis = "cell", factors = .pe_axis_factor_payload(x, f, "cell")))
    out$deg_sum <- .stream_expr_accum(by_deg(deg), stats = "sum")$sum
    out$deg_sumsq <- .stream_expr_accum(by_deg(sqrt(deg)), stats = "sumsq")$sumsq
    out
}
