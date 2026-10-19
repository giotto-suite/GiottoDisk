#' @include pkg_imports.R class-parquetExprStore.R
NULL

# Binarization on an expression store: thresholdParam("binarize") from Giotto.
# A numeric threshold (one, or one per feature) is the `compare` op as it
# stands. A rule ("kmeans", "rank") is solved per feature with Giotto's
# histThreshold() on the streamed value histogram, so the store is read once
# for the histogram and the result stays lazy.
#
# The rule's threshold is a bin edge, so it is applied in the histogram's own
# units: the same per-feature `multiply` that binned the values, then a
# `compare` against the edge. Converting the edge back to value units instead
# would round differently from the binning and move entries lying on an edge.
# The binned split equals the exact in-memory one whenever each feature's
# distinct positive values fall in distinct bins -- raw counts up to `n_bins`,
# for instance.

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "binarizeThreshParam"),
    function(x, param, ...) {
        # thresholds are per feature, whichever way the store is presented
        if (isTRUE(x@transposed)) return(t(processData(t(x), param, ...)))
        thr <- param$threshold %null% 0
        if (!is.character(thr)) return(x > thr)
        r <- .pe_rule_thresholds(x, thr,
            percentage_rank = param$percentage_rank %null% 30,
            n_bins = param$n_bins %null% 1024L)
        scaled <- .pe_push_op_last(x, list(type = "multiply", axis = "feat",
            factors = .pe_axis_factor_payload(x, r$scale, "feat")))
        scaled > r$edge
    }
)

# Per feature, the scale onto histogram bins (`n_bins / max`, 0 for a feature
# storing nothing positive) and the threshold as a bin edge in those units.
# Bins are represented by their midpoints; the cells storing nothing (or <= 0)
# are a bin at 0.
.pe_rule_thresholds <- function(x, rule, percentage_rank = 30, n_bins = 1024L) {
    h <- .pe_value_hist(x, n_bins = n_bins)
    edge <- rep(Inf, length(h$max))
    by_pos <- split(h$counts, h$counts$pos)
    for (key in names(by_pos)) {
        d <- by_pos[[key]]
        i <- as.integer(key)
        edge[i] <- histThreshold(
            value = c(0, d$bin + 0.5),
            count = c(h$n_zero[i], d$n),
            upper = c(0, d$bin + 1),
            method = rule, percentage_rank = percentage_rank)
    }
    list(scale = ifelse(h$max > 0, n_bins / h$max, 0), edge = edge)
}
