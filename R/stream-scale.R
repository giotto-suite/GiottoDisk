#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-scale ####
# Z-score scaling for parquetExprBase, as Giotto's in-memory zscoreScaleParam
# (`base::scale()` per row or column of the matrix): per position of one axis,
# divide by the standard deviation (`scale`) after subtracting the mean
# (`center`). The statistics are computed when the verb runs, over the values
# the chain gives at that point and the population in view, and frozen into
# the records (adr/0006):
#
#   scale   -> a `multiply` record, 1 / sd   (sparse; lowers to Acero)
#   center  -> an `add` record, -mean / sd  (dense; bounded reads only)
#
# A centred store is read only as a slice (adr/0027): the intended consumers
# are displays (heatmaps, the "scaled" slot), and every analysis that needs
# centring folds it into its own algebra instead (PCA's `scale = TRUE`).

# Sums and sums of squares of the chain's values per position of `axis`, and
# `n`, the count of positions on the other axis. The values are dense where the
# chain ends in an `add`, so the accumulators cannot see them; one trailing
# `add` (a single zscoreScaleParam before this one) is resolved algebraically
# from the chain before it, `y`, as v = y + b:
#
#   add on `axis`:   sum = S + n b,     sumsq = Q + 2 b S + n b^2
#   add across it:   sum = S + sum(b),  sumsq = Q + 2 sum_k(y_k b_k) + sum(b^2)
#
# with sum_k(y_k b_k) the sum of `y` after a `multiply` by b. A chain with an
# `add` anywhere else is refused: storeWrite() it first.
.pe_axis_moments <- function(x, axis) {
    other <- if (identical(axis, "feat")) "cell" else "feat"
    n <- as.numeric(if (identical(other, "feat")) x@n_genes else x@n_cells)
    chain <- c(x@ops, x@post_ops)
    is_add <- vapply(chain, function(o) identical(o$type, "add"), NA)
    if (!any(is_add)) {
        st <- .stream_expr_accum(x, axis = axis, stats = c("sum", "sumsq"))
        return(list(sum = st$sum, sumsq = st$sumsq, n = n))
    }
    if (sum(is_add) > 1L || !is_add[[length(is_add)]] || !length(x@post_ops)) {
        stop("[processData(parquetExprBase, zscoreScaleParam)] statistics over ",
             "a chain with a centring step other than the last cannot be ",
             "streamed. storeWrite() the store first.", call. = FALSE)
    }
    add <- x@post_ops[[length(x@post_ops)]]
    y <- x
    y@post_ops <- y@post_ops[-length(y@post_ops)]
    st <- .stream_expr_accum(y, axis = axis, stats = c("sum", "sumsq"))
    b <- .pe_payload_view(y, add$terms, add$axis %||% "cell")
    if (identical(add$axis, axis)) {
        return(list(sum = st$sum + n * b,
            sumsq = st$sumsq + 2 * b * st$sum + n * b^2, n = n))
    }
    yb <- .pe_push_op_last(y, list(type = "multiply", axis = other, factors = add$terms))
    cross <- .stream_expr_accum(yb, axis = axis, stats = "sum")$sum
    list(sum = st$sum + sum(b), sumsq = st$sumsq + 2 * cross + sum(b^2), n = n)
}

# A payload keyed by on-disk id, as a vector in view order along `axis`
.pe_payload_view <- function(pe, payload, axis) {
    if (!is.list(payload)) {
        n <- if (identical(axis, "feat")) pe@n_genes else pe@n_cells
        return(rep_len(as.numeric(payload), n))
    }
    map <- .pe_axis_pos_map(pe, axis)
    uid <- if ("source_id" %in% names(map)) map$source_id
        else rep(as.character(.exprbase_substores(pe)[[1L]]$store@uid), nrow(map))
    out <- numeric(max(map$pos))
    for (u in unique(uid)) {
        i <- which(uid == u)
        out[map$pos[i]] <- as.numeric(payload[[u]])[map$key_id[i]]
    }
    out
}

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "zscoreScaleParam"),
    function(x, param, ...) {
        margin <- param$MARGIN %null% 2
        if (!isTRUE(margin %in% c(1, 2))) {
            stop("[processData(parquetExprBase, zscoreScaleParam)] `MARGIN` ",
                 "must be 1 (rows) or 2 (columns).", call. = FALSE)
        }
        center <- isTRUE(param$center %null% TRUE)
        scale  <- isTRUE(param$scale %null% TRUE)
        if (!center && !scale) return(x)
        # MARGIN names a matrix axis; which store axis that is follows the
        # store's orientation (adr/0016)
        axis <- .pe_orient(x, "feat", "cell")[[margin]]

        m <- .pe_axis_moments(x, axis)
        n <- m$n
        mean <- m$sum / n
        # As base::scale(): the sd after centring, the root mean square
        # without it, both over n - 1. A constant gives NaN, as in memory.
        sd <- if (center) sqrt((m$sumsq - n * mean^2) / (n - 1))
            else sqrt(m$sumsq / (n - 1))
        if (center) sd[!is.na(sd) & sd < 0] <- 0   # rounding below 0

        if (scale) {
            x <- .pe_push_op_last(x, list(type = "multiply", axis = axis,
                factors = .pe_axis_factor_payload(x, 1 / sd, axis, fill = NA_real_)))
        }
        if (center) {
            terms <- if (scale) -mean / sd else -mean
            x <- .pe_push_op(x, list(type = "add", axis = axis,
                terms = .pe_axis_factor_payload(x, terms, axis, fill = NA_real_)),
                phase = "post")
        }
        x
    }
)

# As Giotto's defaultScaleParam: z-score features (MARGIN 1) and cells
# (MARGIN 2), in `scale_order`
#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "defaultScaleParam"),
    function(x, param, ...) {
        feats <- Giotto::scaleParam("zscore", center = TRUE, scale = TRUE, MARGIN = 1)
        cells <- Giotto::scaleParam("zscore", center = TRUE, scale = TRUE, MARGIN = 2)
        steps <- list()
        if (isTRUE(param$scale_feats) && isTRUE(param$scale_cells)) {
            order <- match.arg(param$scale_order, c("first_feats", "first_cells"))
            steps <- if (order == "first_feats") list(feats, cells) else list(cells, feats)
        } else if (isTRUE(param$scale_feats)) {
            steps <- list(feats)
        } else if (isTRUE(param$scale_cells)) {
            steps <- list(cells)
        }
        for (p in steps) x <- processData(x, p)
        x
    }
)
