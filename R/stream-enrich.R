#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-enrich ####
#
# analyzeData(parquetExprBase, pageEnrichParam) -- streaming PAGE enrichment.
#
# What makes PAGE streamable is that its statistic is a function of moments,
# not of the matrix. In memory it builds
#
#     geneFold = expr - mean_gene_expr        # G x N, DENSE
#
# and then only ever asks three things of it: the per-cell mean, the per-cell
# sd, and the per-cell mean over each cell type's marker rows. All three are
# recoverable from per-cell sums of the STORED (sparse) values plus the
# per-gene vector m, so `geneFold` never has to exist:
#
#     S1_c = sum_g (x_gc - m_g)   = colsum_c - sum(m)
#     S2_c = sum_g (x_gc - m_g)^2 = colsumsq_c - 2 * sum_g x_gc m_g + sum_g m_g^2
#     V_ct = mean_{g in M_ct} (x_gc - m_g)
#          = (sum_{g in M_ct} x_gc - sum_{g in M_ct} m_g) / |M_ct|
#
# The only term that is not a plain marginal is `sum_g x_gc m_g`, which is a
# per-cell sum after a per-feature `multiply` -- sparsity-preserving, so it
# lowers into the same accumulator as the others.
#
# Passes over the store: one on the feature axis for `m`, one on the cell axis
# for sum and sumsq, one on the cell axis for the m-weighted sum, and one per
# cell type over that type's marker rows only. Cell-type passes read |M_ct|
# rows, not G, so they are cheap next to the three full ones.
#
# NOT streamed here:
#
#   p_value = TRUE. The permutation branch computes a marker-set mean for
#   n_times x cell_types random gene sets -- 20,000 additional per-cell sums
#   at the defaults. That is a different algorithm, not a different backend,
#   and it errors rather than quietly costing a thousand passes.
#
#   rank and hypergeometric. Both are refused with a clear message; see the
#   methods at the bottom of this file for why each one is not a moment.
#
# Exactness: the accumulators are bit-identical to their in-memory
# counterparts, and `.op_transform_expm1` was written as `base^value - 1`
# precisely so that `m` is too. The sd is not: computing sum((d - dbar)^2)
# from raw moments sums in a different order than `stats::sd`, so it agrees to
# floating point rather than to the last bit. Everything downstream divides by
# it, so the parity tests assert a relative tolerance and say what it is.


# The topic owns the `@name` on its own block. A `@name` on a setMethod block
# suppresses roxygen's \alias{analyzeData,...-method}, and R CMD check then
# reports that method as undocumented -- which is why the two refusal methods
# below, which use only `@rdname`, were fine and this one was not.

#' @name analyzeData-pageEnrichParam
#' @title Streaming PAGE enrichment
NULL

#' @rdname analyzeData-pageEnrichParam
#' @title Streaming PAGE enrichment
#' @description
#' [GiottoClass::analyzeData()] method for a [Giotto::pageEnrichParam-class] on
#' a disk-backed expression store. Nothing is densified: PAGE's statistic is a
#' function of per-gene and per-cell moments, which stream.
#'
#' `p_value = TRUE` is not supported here -- it is a permutation over
#' `n_times` x cell-types random gene sets, which is thousands of extra passes
#' rather than a heavier one. Run it on an in-memory object.
#' @param x a `parquetExprBase` store.
#' @param param a [Giotto::pageEnrichParam-class].
#' @param sign_matrix binary sign matrix, genes x cell types.
#' @param ... additional arguments (none used).
#' @returns a `data.table` of `cell_ID` and one column per cell type, the same
#'   contract the in-memory method returns.
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "pageEnrichParam"),
    function(x, param, ..., sign_matrix) {
        .stream_page(pe = x, sign_matrix = sign_matrix, param = param)
    }
)


# Build a `multiply` payload for a per-feature vector given in VIEW order.
#
# The op registry keys factors by on-disk id per substore uid (see the
# `multiply` note in utils-pestore-ops.R), and `.pe_axis_pos_map()` is the
# view-position -> on-disk-key map that every other consumer uses. A single
# store's map carries no `source_id`, so the payload is keyed by its own uid.
.pe_feat_factor_payload <- function(pe, w) {
    map <- .pe_axis_pos_map(pe, "feat")
    subs <- .exprbase_substores(pe)
    if (!"source_id" %in% names(map)) {
        uid <- as.character(subs[[1L]]$store@uid)
        map <- data.table::copy(map)[, "source_id" := uid]
    }
    split_map <- split(map, map$source_id)
    lapply(split_map, function(d) {
        v <- numeric(max(d$key_id))
        v[d$key_id] <- w[d$pos]
        v
    })
}


.stream_page <- function(pe, sign_matrix, param) {
    if (isTRUE(param$p_value)) {
        stop("[analyzeData(parquetExprBase, pageEnrichParam)] ",
            "`p_value = TRUE` is not supported on a disk-backed store. The ",
            "permutation branch draws n_times x cell_types random gene sets ",
            "and takes a marker-set mean for each, which is thousands of ",
            "additional passes rather than one heavier one. Run it on an ",
            "in-memory object, or use the scores with p_value = FALSE.",
            call. = FALSE)
    }

    verbose <- isTRUE(param$verbose)
    all_genes <- pe@feat_ids
    cell_ids <- pe@cell_ids
    n_cells <- length(cell_ids)
    n_genes <- length(all_genes)

    sign_matrix <- as.matrix(sign_matrix)

    # ---- cell-type selection, identical to the in-memory method ------------
    # detected = column sums of the sign matrix over genes the store holds
    in_store <- rownames(sign_matrix) %in% all_genes
    detected <- colSums(sign_matrix[in_store, , drop = FALSE])
    lost <- names(detected)[detected <= param$min_overlap_genes]
    for (ct in lost) {
        if (verbose) {
            print(paste0("Warning, ", ct, " only has ", detected[[ct]],
                " overlapping genes. Will be removed."))
        }
    }
    available_ct <- names(detected)[detected > param$min_overlap_genes]
    if (length(available_ct) == 1L) stop("Only one cell type available.")

    interGene <- intersect(rownames(sign_matrix), all_genes)
    filterSig <- sign_matrix[interGene, available_ct, drop = FALSE]

    # ---- pass 1: the per-gene reference level ------------------------------
    if (isTRUE(param$reverse_log_scale)) {
        base <- param$logbase
        src <- .pe_push_op(pe, list(type = "expm1", base = base))
        m <- log(.stream_expr_accum(src, axis = "feat", stats = "sum")$sum /
                     n_cells + 1)
    } else {
        m <- .stream_expr_accum(pe, axis = "feat", stats = "sum")$sum / n_cells
    }
    names(m) <- all_genes
    sum_m <- sum(m)
    sum_m2 <- sum(m * m)

    # ---- pass 2: per-cell sum and sum of squares ---------------------------
    acc <- .stream_expr_accum(pe, axis = "cell", stats = c("sum", "sumsq"))

    # ---- pass 3: per-cell sum weighted by the reference level --------------
    weighted <- .pe_push_op(pe, list(
        type = "multiply", axis = "feat",
        factors = .pe_feat_factor_payload(pe, m)
    ))
    xm <- .stream_expr_accum(weighted, axis = "cell", stats = "sum")$sum

    # geneFold's per-cell mean and sd, from the moments above. `n_genes - 1`
    # matches stats::sd, which is what apply(geneFold, 2, sd) calls.
    s1 <- acc$sum - sum_m
    s2 <- acc$sumsq - 2 * xm + sum_m2
    colmean <- s1 / n_genes
    colsd <- sqrt(pmax(s2 - s1 * s1 / n_genes, 0) / (n_genes - 1))

    # ---- pass 4..T+3: one per cell type, over its marker rows only ---------
    zs <- vapply(available_ct, function(ct) {
        markers <- interGene[filterSig[, ct] == 1]
        k <- length(markers)
        sub <- pe[markers, ]
        got <- .stream_expr_accum(sub, axis = "cell", stats = "sum")$sum
        v1 <- (got - sum(m[markers])) / k
        (v1 - colmean) * sqrt(k) / colsd
    }, FUN.VALUE = numeric(n_cells))
    zs <- matrix(zs, nrow = n_cells, ncol = length(available_ct),
                 dimnames = list(cell_ids, available_ct))

    if (identical(param$output_enrichment, "zscore")) {
        zs <- apply(zs, 2, function(v) as.numeric(scale(v)))
        dimnames(zs) <- list(cell_ids, available_ct)
    }

    out <- data.table::data.table(cell_ID = cell_ids)
    out <- cbind(out, data.table::as.data.table(zs))
    # Ordered by cell_ID, because the in-memory method's final dcast() is: the
    # two backends must hand back the same table, not the same table up to a
    # permutation.
    data.table::setorderv(out, "cell_ID")
    out[]
}


# ---- rank and hypergeometric: refused, with the reason ---------------------

#' @rdname analyzeData-pageEnrichParam
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "rankEnrichParam"),
    function(x, param, ..., sign_matrix) {
        stop("[analyzeData(parquetExprBase, rankEnrichParam)] ",
            "rank enrichment is not supported on a disk-backed store. It ",
            "ranks each gene across every cell (sparseMatrixStats::rowRanks), ",
            "so no cell chunk can be scored without the others -- unlike ",
            "PAGE, whose statistic is a function of moments. Use ",
            "enrich_method = \"PAGE\", or materialize the store first.",
            call. = FALSE)
    }
)

#' @rdname analyzeData-pageEnrichParam
#' @export
setMethod("analyzeData",
    signature(x = "parquetExprBase", param = "hyperEnrichParam"),
    function(x, param, ..., sign_matrix) {
        stop("[analyzeData(parquetExprBase, hyperEnrichParam)] ",
            "hypergeometric enrichment is not supported on a disk-backed ",
            "store yet. Its per-cell top-percentage cutoff is a quantile ",
            "over all genes, which is chunkable by cell but is not a moment ",
            "and so shares nothing with the accumulators PAGE uses. Use ",
            "enrich_method = \"PAGE\", or materialize the store first.",
            call. = FALSE)
    }
)
