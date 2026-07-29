#' @include class-parquetExprStore.R
NULL

# parquetExprStore op chain — two-phase design.
#
# @ops       arrow-lazy phase. Ops that translate to arrow-dplyr steps
#            on the lazy query (filter, distinct, head, tail, sample,
#            future arrow-native transforms). Composed at storeRead time
#            via .pe_apply_ops before collect.
#
# @post_ops  R-side post-collect phase. Ops that operate on the materialized
#            data.table, for every output mode and for streaming consumers
#            alike -- they all reach it through storeRead. Composed at
#            output-materialization time via .pe_apply_post_ops_df.
#
# Each op is a pure-data record `list(type = <character>, ...params)`.
# No closures, no phase field on the record — the phase is determined by
# which slot the op lives in.
#
# Monotonic phase rule: once @post_ops has an entry, subsequent pushes go
# to @post_ops regardless of the op's natural phase. Enforced by
# .pe_push_op. Users needing a lazy op after a post op must materialize
# first (storeWrite bakes @post_ops into on-disk values; new pe starts
# with empty chains).
#
# Op types currently supported:
#
#   norm_libsize_log  (phase: post)
#     Fused library-size scaling + optional log transform. Mutates
#     `value` in place, R-side.
#     Params:
#       scalef   data.table(source_id, orig_row_id, scalef). One row per
#                cell present in the store's current view at processData
#                time. Composite (source_id, orig_row_id) key so the same
#                op record carries scalefs for all substores of a union.
#       log      logical. Apply log1p / log(base) after scaling.
#       base     numeric. log base (default 2).
#
# Extension protocol:
#   - Lazy op: add a branch to .pe_apply_op (arrow-side switch), and a
#     .pe_op_table_keys entry if the op carries axis-keyed state that
#     `[`-subset should narrow.
#   - Post op: add branches to .pe_apply_post_op_df AND
#     .pe_apply_post_op_mat (both shapes). Add .pe_op_table_keys entry
#     if axis-keyed. Add a natural-phase-"post" case in the verb that
#     produces it.


# ---- @ops arrow-side executor ---------------------------------------------

# Apply a single arrow-lazy op record to a lazy arrow query.
# Returns the augmented query.
.pe_apply_op <- function(atab, op) {
    switch(op$type,
        # Future arrow-native op types dispatch here. norm_libsize_log
        # moved to @post_ops; not an arrow-side op.
        stop("[.pe_apply_op] unknown arrow-side op type: ", op$type,
            call. = FALSE)
    )
}

# Fold the arrow-side chain — composes all @ops into one lazy query.
.pe_apply_ops <- function(atab, ops) {
    if (length(ops) == 0L) return(atab)
    for (op in ops) atab <- .pe_apply_op(atab, op)
    atab
}


# ---- @post_ops R-side executor (data.table shape) --------------------------
#
# Used by materializing output paths (tibble / data.table / dgcmatrix) and
# by any consumer that collects a triplet chunk into a data.table. Ops
# mutate `df$value` in place.

.pe_apply_post_op_df <- function(df, op) {
    switch(op$type,
        "norm_libsize_log" = .pe_apply_post_op_norm_libsize_log_df(df, op),
        stop("[.pe_apply_post_op_df] unknown post op type: ", op$type,
            call. = FALSE)
    )
}

.pe_apply_post_ops_df <- function(df, post_ops) {
    if (length(post_ops) == 0L) return(df)
    for (op in post_ops) df <- .pe_apply_post_op_df(df, op)
    df
}

.pe_apply_post_op_norm_libsize_log_df <- function(df, op) {
    # NSE bindings
    row_id <- source_id <- value <- scalef <- NULL
    # For a single-source scalef payload, replace the (source_id,
    # orig_row_id) update-join with a plain positional vector index:
    # build a numeric vector where `scalef_by_row_id[orig_row_id]` is
    # the scalef for that on-disk cell, then look up per df row via
    # `scalef_by_row_id[df$row_id]`.  No composite hash, no per-row
    # join dispatch — a native vector index.
    # Union collected chunks (multi-source scalef) still use the
    # update-join because source_id disambiguates overlapping row_ids
    # across substores.
    sids <- unique(op$scalef$source_id)
    if (length(sids) == 1L) {
        n_slots <- max(op$scalef$orig_row_id)
        scalef_by_row_id <- numeric(n_slots)
        scalef_by_row_id[op$scalef$orig_row_id] <- op$scalef$scalef
        df[, value := value * scalef_by_row_id[row_id]]
    } else {
        df[op$scalef,
            on = c("source_id", "row_id" = "orig_row_id"),
            scalef := i.scalef]
        df[, value := value * scalef]
        df[, scalef := NULL]
    }
    if (isTRUE(op$log)) {
        df[, value := log1p(value) / log(op$base %null% 2)]
    }
    df
}


# ---- Push helper: monotonic phase enforcement ------------------------------

# Route an op record to the appropriate slot on the store. Enforces the
# monotonic-phase rule: once @post_ops has an entry, subsequent lazy
# pushes are ALSO routed to @post_ops (bumped) or rejected — we currently
# choose reject (strict), which is easier to debug and steers users to
# materialize when they hit the boundary.
#
# For our current op inventory (norm on @post_ops, no lazy pestore ops
# yet), the "phase = lazy" path is unused. When we add filter or similar
# lazy pestore ops, this helper is the enforcement point.
.pe_push_op <- function(pe, op, phase = c("lazy", "post")) {
    phase <- match.arg(phase)
    if (phase == "lazy" && length(pe@post_ops) > 0L) {
        stop("[.pe_push_op] cannot queue a lazy op after a post op is ",
             "already on @post_ops. Materialize first (via storeWrite) ",
             "to reset the chain, then push the lazy op.", call. = FALSE)
    }
    if (phase == "post") {
        pe@post_ops <- c(pe@post_ops, list(op))
    } else {
        pe@ops <- c(pe@ops, list(op))
    }
    pe
}


# ---- parquetExprBase substore iteration protocol ---------------------------
#
# Stream-pipeline methods (filterData, processData, analyzeData, ...) that
# work uniformly over both `parquetExprStore` and `unionParquetExprStore`
# dispatch on the shared `parquetExprBase` virtual and iterate via this
# protocol.
#
# Returns a list of substore-entry records. Each entry is a list:
#
#   $store        : the substore (always a parquetExprStore)
#   $cell_offset  : 0-based offset of this substore's cells in the
#                   union's `@cell_ids` axis (always 0 for a single
#                   parquetExprStore; cumulative substore offset for a
#                   unionParquetExprStore)

# Project a parent (union) @ops chain onto a single substore so its
# `storeRead()` carries the same arrow-lazy recipe restricted to this
# substore's rows. For arrow-native ops with source-keyed payload the
# per-substore filter is applied here. norm_libsize_log now lives on
# @post_ops (not @ops), so this projection currently no-ops for it —
# @post_ops are consumed by streaming consumers via .pe_scalef_vec_for_sub
# and don't need to travel with the substore's own @ops.
# Project a union parent's phase chains onto one of its substores, so the
# substore alone is enough to read from.  Union substores carry no ops by
# constraint (see the `[` method for unionParquetExprStore) -- only the union
# does -- so anything reading a substore directly has to transplant them.
#
# `parent_post_ops` is optional and only injected when the substore has none
# of its own: for a single (non-union) store `.exprbase_substores()` yields
# the store itself, which already carries its @post_ops, and re-adding them
# would double-apply.
#
# Injecting @post_ops here rather than at chunk-read time is what lets `[`
# slice them: narrowing the substore's cell axis runs the surviving keys
# (source_id = this substore's uid, orig_row_id = kept cells) through
# `.pe_slice_ops`, which inner-joins the parent's global scalef table down to
# just this chunk's rows.  That leaves a single source_id, so the executor
# takes the positional-lookup branch instead of the keyed-join branch.
.exprbase_inject_parent_ops <- function(sub, parent_ops,
                                         parent_post_ops = list()) {
    if (length(parent_ops) > 0L) {
        sub@ops <- c(sub@ops, parent_ops)
    }
    if (length(parent_post_ops) > 0L && length(sub@post_ops) == 0L) {
        sub@post_ops <- c(sub@post_ops, parent_post_ops)
    }
    sub
}

.exprbase_substores <- function(x) {
    if (inherits(x, "unionParquetExprStore")) {
        offsets <- c(0L, cumsum(vapply(x@stores,
            function(s) as.integer(s@n_cells), integer(1L))))
        return(lapply(seq_along(x@stores), function(i) {
            list(store = x@stores[[i]],
                 cell_offset = offsets[i])
        }))
    }
    if (inherits(x, "parquetExprStore")) {
        return(list(list(store = x, cell_offset = 0L)))
    }
    stop("[.exprbase_substores] expected a parquetExprBase, got ",
        toString(class(x)), call. = FALSE)
}


# ---- shared chunk reader (used by storeWrite baking) ------------------------
#
# Reads cells [sub_cs, sub_ce] from `info$sub` and returns the normalized
# chunk, or NULL when the chunk holds no nonzeros.
#
# Goes through the framework verbs -- `[` to narrow the cell axis, `storeRead`
# to filter, materialize and apply @post_ops -- rather than hand-rolling a
# query.  Three things fall out of that:
#   * the gene predicate comes from `sub@gene_idx`, which `storeRead` already
#     applies; the old explicit `col_id %in% hvg_orig` filter duplicated it.
#   * a contiguous cell chunk is the gapless case in `.pe_axis_pred()`, so the
#     cell predicate becomes a pure `row_id >= lo, row_id <= hi` range that
#     prunes parquet row groups, instead of an `is_in` over the chunk's ids.
#   * @post_ops application lives in one place instead of two.
#
# `info$sub` must already carry both phase chains --
# `.exprbase_inject_parent_ops()` transplants a union parent's @ops and
# @post_ops onto the substore -- so `[` slices @post_ops to this chunk's cells
# and `storeRead` applies them.
#
# Returns genes x cells (Bioconductor convention); callers index accordingly
# rather than materializing a `t()`.
#
# `info` is a list whose only field this reader touches is `$sub`, the
# parquetExprStore substore with both op chains already injected. Callers also
# carry `$hvg_orig` and `$scalef_vecs`, and pass `post_ops` / `P_hvg`, none of
# which are read here now that `storeRead` owns the gene filter and the
# @post_ops apply; they are kept so `info` and the call signature stay one
# shape across readers.
.pe_read_chunk_sub <- function(info, sub_cs, sub_ce, post_ops, P_hvg) {
    M <- storeRead(info$sub[, sub_cs:sub_ce], output = "dgcmatrix",
                   max_rows = Inf, max_cols = Inf)
    if (length(M@x) == 0L) return(NULL)
    M
}


# ---- introspection helpers -------------------------------------------------

# Find the index of the first op of a given type in a chain. Returns
# NA_integer_ when not present. Chain-agnostic — pass @ops or @post_ops.
.pe_find_op_type <- function(ops, type) {
    if (length(ops) == 0L) return(NA_integer_)
    idx <- which(vapply(ops, function(op) identical(op$type, type),
        logical(1L)))
    if (length(idx) == 0L) NA_integer_ else as.integer(idx[1L])
}

# Predicate: does the store carry a normalization recipe? Norm now lives
# on @post_ops.
.pe_has_norm_op <- function(pe) {
    !is.na(.pe_find_op_type(pe@post_ops, "norm_libsize_log"))
}


# ---- subset slice ----------------------------------------------------------
#
# Ops (both phases) are frozen snapshots of population stats captured at
# trigger time. Subsetting cells / genes never invalidates them — every
# op with axis-keyed lookup tables gets those tables filtered by the
# surviving axis keys. Survivors retain their captured-population stats;
# the user re-runs processData(...) if they want stats over the new
# population.
#
# Slice dispatch walks the registry and filters rows.

# Registry: per op type, list of (table_field_name -> list of axes -> key
# columns). The key columns are the field names in the lookup table that
# match the parquet's hive-partition / row-id schema.
#   cell axis key: c("source_id", "orig_row_id") composite — row_id
#                  restarts per substore so source_id disambiguates
#   gene axis key: "feat_id" — names, since col_id positional layout can
#                  differ across substores (alignment not guaranteed)
.pe_op_table_keys <- list(
    norm_libsize_log = list(
        scalef = list(cell = c("source_id", "orig_row_id"))
    )
    # Future kinds add their entries here.
)


# Slice one op record along one axis. Returns the updated op (with axis-
# keyed tables filtered to surviving keys). Ops that don't have a table
# on the given axis are returned unchanged.
#
# `surviving_keys` is a data.table:
#   axis == "cell": data.table(source_id, orig_row_id) of surviving cells
#   axis == "gene": data.table(feat_id) of surviving features
.pe_slice_op <- function(op, axis, surviving_keys) {
    spec <- .pe_op_table_keys[[op$type]]
    if (is.null(spec)) return(op)  # unknown kind: leave untouched
    for (tname in names(spec)) {
        axis_spec <- spec[[tname]][[axis]]
        if (is.null(axis_spec)) next  # this table isn't axis-keyed
        tbl <- op[[tname]]
        if (is.null(tbl) || nrow(tbl) == 0L) next
        join_cols <- axis_spec
        op[[tname]] <- tbl[surviving_keys, on = join_cols,
                            nomatch = NULL, allow.cartesian = FALSE]
    }
    op
}


# Slice every op in a chain along one axis. Chain-agnostic; caller passes
# @ops or @post_ops as `ops`.
.pe_slice_ops <- function(ops, axis, surviving_keys) {
    if (length(ops) == 0L) return(ops)
    lapply(ops, .pe_slice_op, axis = axis, surviving_keys = surviving_keys)
}
