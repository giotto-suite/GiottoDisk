#' @include pkg_imports.R class-parquetExprStore.R
NULL

# analyzeData(parquetExprBase, featCorParam): the feature x feature Pearson
# correlation from one Gram pass (.pe_gram_pass, shared with stream PCA)
# over the store materialized with its chain baked. The result is dense
# features x features, so pass a feature subset (`x[feats, ]`).

#' @rdname analyzeData
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "featCorParam"),
    function(x, param, ...) {
        if (!identical(param$method, "pearson")) {
            stop("[analyzeData(parquetExprBase, featCorParam)] only ",
                "method = \"pearson\" streams; ranks need every value of a ",
                "feature at once.", call. = FALSE)
        }
        tmp <- .dump_tempfile()
        on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
        pe <- storeWrite(parquetExprStore(path = tmp), x)
        P <- length(pe@feat_ids)
        gp <- .pe_gram_pass(pe, P)
        n <- as.numeric(pe@n_cells)
        cov <- (gp$G - tcrossprod(gp$s) / n) / (n - 1)
        v <- diag(cov)
        sdv <- sqrt(ifelse(v > 0, v, NA_real_))
        r <- pmin(pmax(cov / outer(sdv, sdv), -1), 1)
        diag(r)[!is.na(sdv)] <- 1
        dimnames(r) <- list(pe@feat_ids, pe@feat_ids)
        r
    }
)
