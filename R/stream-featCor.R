#' @include pkg_imports.R class-parquetExprStore.R
NULL

# analyzeData(parquetExprBase, featCorParam): the feature x feature Pearson
# correlation from x x^T and the row sums (GiottoClass's
# feat_cor_from_moments). Plain values take one Gram pass (.pe_gram_pass,
# shared with stream PCA) over the store materialized with its chain baked.
# With `smooth =` the sums are taken window by window over the smoothed
# values (.pe_smooth_moments), which are never written. The result is dense
# features x features, so pass a feature subset (`x[feats, ]`).

#' @rdname analyzeData
#' @param network for a `featCorParam` with `smooth`: the spatial network as
#'   a `parquetEdgeStore`.
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "featCorParam"),
    function(x, param, ..., network = NULL) {
        if (!identical(param$method, "pearson")) {
            stop("[analyzeData(parquetExprBase, featCorParam)] only ",
                "method = \"pearson\" streams; ranks need every value of a ",
                "feature at once.", call. = FALSE)
        }
        n <- as.numeric(x@n_cells)
        if (!is.null(param$smooth)) {
            prep <- .pe_smooth_prepare(x, param$smooth, network)
            on.exit(unlink(prep$tmp, recursive = TRUE), add = TRUE)
            m <- .pe_smooth_moments(prep, x@feat_ids)
        } else {
            tmp <- .dump_tempfile()
            on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
            pe <- storeWrite(parquetExprStore(path = tmp), x)
            m <- .pe_gram_pass(pe, length(pe@feat_ids))
        }
        feat_cor_from_moments(m$G, m$s, n, x@feat_ids)
    }
)
