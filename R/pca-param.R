#' @include pkg_imports.R
NULL

# gramEigenPcaParam — streaming-only PCA flavor. Class + constructor
# live here (not in Giotto) because the technique's advantage (2 disk
# passes for exact top-k) is meaningful only out-of-core. Constructor
# is exported directly rather than routed through Giotto::pcaParam()
# because it carries a knob the factory doesn't expose (fallback_relerr).

#' @rdname gramEigenPcaParam
#' @exportClass gramEigenPcaParam
setClass("gramEigenPcaParam", contains = "pcaParam")


#' @name gramEigenPcaParam
#' @title Gram-eigen streaming PCA parameter
#' @description
#' Streaming PCA that accumulates the P × P Gram matrix `AᵀA` in two
#' disk passes and eigendecomposes it in memory. Forming `AᵀA` squares
#' the condition number; when the predicted rel-err on d_k exceeds
#' `fallback_relerr`, delegates to the Halko path.
#'
#' @param ncp number of components. Default `50`.
#' @param center logical. Center columns. Default `TRUE`.
#' @param scale logical. Per-gene z-score (implicit — σ rescaling is
#'   absorbed into the projection matrix; sparsity preserved, no extra
#'   passes). Default `TRUE`.
#' @param feats_to_use character vector of feature IDs. `NULL` = use all.
#' @param set_seed,seed_number Halko-fallback seed control.
#' @param n_oversamples,n_power_iter Halko-fallback knobs. Defaults `10`, `2`.
#' @param fallback_relerr predicted-relerr threshold above which the
#'   method delegates to Halko. Default `0.01`.
#' @param ... reserved.
#' @return A `gramEigenPcaParam` object.
#' @seealso [pestore-chunking] for how the passes are windowed and which
#'   options bound their memory.
#' @examples
#' # p <- gramEigenPcaParam(ncp = 30, feats_to_use = hvg_ids)
#' # res <- reduceData(pes, p)
#' @export
gramEigenPcaParam <- function(ncp             = 50L,
                              center          = TRUE,
                              scale           = TRUE,
                              feats_to_use    = NULL,
                              set_seed        = TRUE,
                              seed_number     = 1234L,
                              n_oversamples   = 10L,
                              n_power_iter    = 2L,
                              fallback_relerr = 0.01,
                              ...) {
    p <- new("gramEigenPcaParam", param = list(...))
    p$method          <- "gram-eigen"
    p$ncp             <- as.integer(ncp)
    p$center          <- isTRUE(center)
    p$scale           <- isTRUE(scale)
    p$feats_to_use    <- feats_to_use
    p$n_oversamples   <- as.integer(n_oversamples)
    p$n_power_iter    <- as.integer(n_power_iter)
    p$set_seed        <- isTRUE(set_seed)
    p$seed_number     <- as.integer(seed_number)
    p$fallback_relerr <- as.numeric(fallback_relerr)
    p
}
