# stream-network ####
# Approximate nearest-neighbor search for network construction.
#
# GiottoClass builds kNN / sNN networks through `.net_dt_knn()` /
# `.net_dt_snn()`, both of which call `dbscan::kNN()`. That is exact and
# single-threaded, and it degrades toward brute force as dimensionality rises
# -- at the ~15-50 PCs a network is normally built on it is the dominant cost
# of `createNearestNetwork()`.
#
# `hnswKNN()` is the drop-in alternative: an HNSW index (hnswlib, via
# RcppHNSW) built once and queried in one shot, parallel over `n_threads`.
# It returns the same shape `dbscan::kNN()` does -- a `c("kNN", "NN")` object
# carrying `id` / `dist` / `k` / `sort` / `metric` -- so `dbscan::sNN()`
# consumes it directly and nothing downstream of the search changes.
#
# The trade is exactness: HNSW is approximate, so recall is high but below
# 1.0 and results can shift with thread count. It is therefore opt-in
# (`engine = "hnsw"` on the network params), never the default.

#' @title Approximate k-nearest neighbors via HNSW
#' @name hnswKNN
#' @description
#' Find the `k` nearest neighbors of every row of `x` using an HNSW index
#' (Hierarchical Navigable Small World), returning the same structure as
#' [dbscan::kNN()] so the two are interchangeable as the search step of
#' network construction.
#'
#' HNSW is *approximate*. Recall is high but not guaranteed to be 1.0, and
#' because the index build and search are multithreaded, results may differ
#' slightly between runs with different `n_threads`. Use [dbscan::kNN()] when
#' exactness matters.
#'
#' @param x numeric matrix. Rows are observations (cells), columns are
#'   dimensions (typically PCA coordinates).
#' @param k integer. Number of neighbors to return per observation, excluding
#'   the observation itself.
#' @param distance character. Metric, one of `"euclidean"` (default),
#'   `"cosine"`, `"l2"` (squared euclidean) or `"ip"` (inner product).
#' @param M integer. HNSW graph degree (default 16). Higher improves recall
#'   at the cost of memory and build time.
#' @param ef_construction integer. Beam width during index construction
#'   (default 200). Higher improves recall at the cost of build time.
#' @param ef integer. Beam width during search (default 50). Higher improves
#'   recall at the cost of query time. Raised to at least `k + 1`.
#' @param n_threads integer. Threads for build and search. Defaults to
#'   [GiottoUtils::determine_cores()]. Both phases scale near-linearly, and
#'   leaving this at 1 makes the search the bottleneck.
#' @param ... unused, for signature compatibility with [dbscan::kNN()].
#' @returns object of class `c("kNN", "NN")` with elements `id` (integer
#'   matrix, `nrow(x)` x `k`), `dist` (numeric matrix, same shape), `k`,
#'   `sort` and `metric`.
#' @examples
#' \dontrun{
#' m <- matrix(rnorm(1000 * 20), nrow = 1000)
#' nn <- hnswKNN(m, k = 30)
#' str(nn$id)
#' }
#' @export
hnswKNN <- function(x,
    k,
    distance = c("euclidean", "cosine", "l2", "ip"),
    M = 16L,
    ef_construction = 200L,
    ef = 50L,
    n_threads = NULL,
    ...
) {
    GiottoUtils::package_check("RcppHNSW")
    distance <- match.arg(distance)

    if (!is.matrix(x)) x <- as.matrix(x)
    checkmate::assert_matrix(x, mode = "numeric")
    k <- as.integer(k)
    n <- nrow(x)
    if (k >= n) {
        stop("[hnswKNN] k (", k, ") must be less than nrow(x) (", n, ").",
             call. = FALSE)
    }

    n_threads <- as.integer(n_threads %null% GiottoUtils::determine_cores())
    # A self-query returns the point itself, so ask for one extra and drop it
    # below. ef must cover the widened request or recall degrades at the tail.
    k_query <- k + 1L
    ef <- max(as.integer(ef), k_query)

    ann <- RcppHNSW::hnsw_build(x,
        distance = distance,
        M = as.integer(M),
        ef = as.integer(ef_construction),
        n_threads = n_threads
    )
    res <- RcppHNSW::hnsw_search(x,
        ann = ann,
        k = k_query,
        ef = ef,
        n_threads = n_threads
    )

    keep <- .hnsw_drop_self(res$idx)
    list_out <- list(
        dist = .hnsw_compact(res$dist, keep, k),
        id = .hnsw_compact(res$idx, keep, k, as_int = TRUE),
        k = k,
        sort = TRUE,
        metric = distance
    )
    structure(list_out, class = c("kNN", "NN"))
}


# Index of the entries to keep after removing each row's self-hit.
#
# The self-hit is normally column 1, but with duplicate coordinates it can
# land anywhere in the row, and with enough duplicates it can be missing
# entirely. So this locates it per row rather than assuming, and falls back to
# dropping the last (furthest) entry when it is absent -- which keeps every
# row exactly k wide either way.
#
# Returns a logical matrix in column-major order, suitable for indexing
# `idx` / `dist` directly.
.hnsw_drop_self <- function(idx) {
    n <- nrow(idx)
    keep <- matrix(TRUE, nrow = n, ncol = ncol(idx))
    self_col <- max.col(idx == seq_len(n), ties.method = "first")
    has_self <- idx[cbind(seq_len(n), self_col)] == seq_len(n)
    # absent self -> drop the furthest neighbor instead
    self_col[!has_self] <- ncol(idx)
    keep[cbind(seq_len(n), self_col)] <- FALSE
    keep
}


# Apply a per-row keep mask, dropping one entry per row.
#
# Done on the transpose: `keep` holds exactly k TRUEs per ROW, so `t(keep)`
# holds exactly k per COLUMN, and a column-major extract from `t(m)` yields
# each row's kept values contiguously. Extracting from `m` directly would
# read column-major across rows whose dropped position differs, so columns
# contribute unequal counts and the reshape silently misaligns rows -- which
# is invisible whenever the self-hit happens to be column 1 in every row.
.hnsw_compact <- function(m, keep, k, as_int = FALSE) {
    out <- t(matrix(t(m)[t(keep)], nrow = k, ncol = nrow(m)))
    if (as_int) storage.mode(out) <- "integer"
    out
}
