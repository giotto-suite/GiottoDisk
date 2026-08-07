#' @include class-parquetExprStore.R
#' @include utils-pestore-ops.R
NULL

# stream-normalize ####
# Streaming JIT normalization for parquetExprBase. Plugs into Giotto's
# existing processData(x, param) dispatch via two setMethod calls on the
# shared `parquetExprBase` virtual:
#
#   processData(parquetExprBase, libraryNormParam)
#       -> appends a `multiply` op (axis = "cell") to x@post_ops
#   processData(parquetExprBase, logNormParam)
#       -> appends an independent `log` op to x@post_ops
#
# The two are separate op records with no ordering requirement between them:
# @post_ops is a chain applied in order, so the chain itself supplies the
# sequencing. Either can be used alone -- log-only on raw counts is a
# legitimate request.
#
# Single (`parquetExprStore`) and union (`unionParquetExprStore`) collapse to
# one algorithm: `.stream_expr_accum(axis = "cell")` returns per-cell sums
# over the whole store in one pass -- a union included, since that verb covers
# every substore in a single Acero plan. The result is sliced back per
# substore by `cell_offset` to build the composite
# `(source_id, orig_row_id, scalef)` lookup table the op record carries.
#
# Neither method rewrites the Parquet file. The recipe lives as a pure-data
# record on @post_ops and is applied R-side after collect (see
# .pe_apply_post_op_multiply_df). The recipe survives saveRDS / load
# cycles without special handling — no closures.
#
# zscoreScaleParam is intentionally NOT implemented for parquetExprBase:
# per-cell / per-gene centering+scaling densifies the sparse matrix and
# breaks the O(N*k) streaming guarantee. normalizeGiotto already errors
# upstream when scale_cells / scale_feats = TRUE on a streaming backend.

# ---- libraryNormParam: compute & store JIT scale factors -------------------

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "libraryNormParam"),
    function(x, param, ...) {
        scalefactor <- param$scalefactor %null% 6e3

        # Append, never edit. The record does its multiplication at the
        # position it occupies; a later call must not reach back and rewrite
        # an earlier one, because arbitrary steps may sit in between. So the
        # factors come from the values as the whole current chain leaves them
        # -- exactly what a new terminal op will multiply.
        #
        # Re-running is therefore self-correcting rather than special-cased:
        # normalizing an already-normalized store yields factors of ~1, and
        # normalizing to a new scalefactor yields the ratio.
        libsizes <- .stream_expr_accum(x, axis = "cell", stats = "sum")$sum
        libsizes[libsizes == 0] <- 1   # guard div-by-zero
        scalef <- as.numeric(scalefactor) / libsizes

        # Payload: one full-length vector per substore, indexed by ON-DISK
        # row_id. Invariant under `[`, so nothing needs slicing later.
        factors <- list()
        for (sub in .exprbase_substores(x)) {
            store <- sub$store
            off   <- as.integer(sub$cell_offset)
            n_sub <- as.integer(store@n_cells)
            ci    <- if (length(store@cell_idx) > 0L) store@cell_idx
                     else seq_len(n_sub)
            v <- rep(NA_real_, max(ci))
            v[ci] <- scalef[off + seq_len(n_sub)]
            factors[[as.character(store@uid)]] <- v
        }

        .pe_push_op(x, 
            list(type = "multiply", axis = "cell", factors = factors),
            phase = "lazy"
        )
    }
)

# ---- logNormParam: set JIT log flag ----------------------------------------

#' @rdname processData
#' @export
setMethod("processData",
    signature(x = "parquetExprBase", param = "logNormParam"),
    function(x, param, ...) {
        base   <- param$base   %null% 2
        offset <- param$offset %null% 1

        if (!isTRUE(offset == 1)) {
            stop("[processData(parquetExprBase, logNormParam)] ",
                 "offset != 1 is not supported for streaming because it ",
                 "would densify the sparse representation. Use offset = 1 ",
                 "(log1p) to preserve sparsity.", call. = FALSE)
        }

        # Stands alone: log-only on raw counts is a legitimate request, and
        # the chain supplies the ordering -- whatever ops precede this one
        # have already been applied by the time the log runs.
        .pe_push_op(x, 
            list(type = "log", base = as.numeric(base)),
            phase = "lazy"
        )
    }
)
