#' @include class-parquetExprStore.R
NULL

# parquetExprStore op chain — one sequence, split at materialization.
#
# The chain is ONE ordered sequence. The two slots record WHERE that sequence
# materializes, not which ops happen to be lowerable:
#
# @ops       the prefix that runs BEFORE materialization. Folded into the lazy
#            arrow query at storeRead time via .pe_apply_ops and executed by
#            Acero as a single plan. Necessarily composed only of ops that
#            lower to arrow -- but that is a CONSEQUENCE of the position, not
#            the definition of the slot.
#
# @post_ops  everything from the first step that cannot run in Acero onward.
#            Applied R-side to the collected data.table via
#            .pe_apply_post_ops_df, for every output mode and for streaming
#            consumers alike -- they all reach it through storeRead.
#
# So a perfectly lowerable op can legitimately sit on @post_ops: if it comes
# after a step that forced materialization, it has nowhere else to go. Reading
# @post_ops as "the ops that cannot be lowered" is the mistake -- it is "the
# suffix that runs after we left Acero."
#
# Each op is a pure-data record `list(type = <character>, ...params)`.
# No closures, no phase field on the record — the phase is determined by
# which slot the op lives in.
#
# Monotonic phase rule: once @post_ops has an entry, a later op cannot go on
# @ops -- it would execute before the post op rather than after it, since the
# fold applies all of @ops, collects, then all of @post_ops. Enforced by
# .pe_push_op. Materializing resets the split (storeWrite bakes the chain into
# on-disk values; the new store starts with both slots empty).
#
# A consumer that wants a lowerable op run R-side edits the chain rather than
# routing around it: .pe_demote_ops moves the op and everything after it to
# @post_ops (same monotonic rule, expressed as a rewrite).
#
# The roles in this file (producer / record / payload / carrier / executor /
# fold / chain editor / consumer) and the invariants connecting them are
# recorded in adr/0004-op-machinery-roles.md. The two that bite most often:
# executors are keyed by (record type, carrier) and NOT by phase, and a
# consumer must never infer what an op is capable of from which slot it is in.
#
# Why the slot means position rather than capability: adr/0005. Why payloads
# are keyed by on-disk id: adr/0003. Why an op whose meaning depends on the
# current window must freeze that statistic when it is pushed, rather than
# consulting the window at read time: adr/0006.
#
# Op types currently supported:
#
#   multiply          (phase: either)
#     Multiply `value` by a per-axis factor. Sparsity-preserving, so it
#     lowers to Acero and applies over a collected triplet frame alike.
#     Params:
#       axis     "cell" | "feat" | "all"
#       factors  scalar (axis "all"), or a named list mapping a substore uid
#                to a numeric vector INDEXED BY ON-DISK ID. Invariant under
#                `[` -- on-disk ids do not move when a view narrows.
#
#   add               (STUB -- recorded and refused, not implemented)
#     Add a per-axis offset. Params mirror `multiply`, with `terms` in place
#     of `factors`. Both executors refuse it; no verb emits one. Intended
#     shape, scope and why it is deferred are in
#     vignettes/articles/roadmap.Rmd, "An `add` op for centred display values".
#
#   log               (phase: post)
#     log1p / log(base). Carries no axis-keyed state.
#     Params:
#       base     numeric. log base (default 2).
#
# Records are positional and self-contained: each does its work at the
# position it occupies, and a verb appends rather than revisiting anything it
# wrote earlier. Nothing needs to be applied in a particular order or to be
# present at all -- the chain supplies the sequencing.
#
# Extension protocol:
#   - Add a branch to .pe_apply_op (arrow) and/or .pe_apply_post_op_df
#     (triplets), depending on which engines can express it.
#   - Key any payload by ON-DISK id, not by view position. That is what makes
#     it invariant under `[` -- see the subset-slice note at the bottom.
#   - Have the producing verb append the record; never edit an existing one.


# ---- @ops arrow-side executor ---------------------------------------------

# Apply a single prefix op record to the lazy arrow query.
# Returns the augmented query.
.pe_apply_op <- function(atab, op) {
    switch(op$type,
        "log"      = .op_transform_log(atab, op),
        "multiply" = .op_multiply(atab, op),
        "add"      = .op_add_refuse(op),
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

.op_transform_log <- function(x, op) {
    value <- NULL # NSE
    base <- op$base %||% 2
    if (data.table::is.data.table(x)) {
        return(x[, value := log1p(value) / log(base)])
    }
    dplyr::mutate(x, value = log1p(value) / log(!!base))
}

# ---- multiply / add ---------------------------------------------------------
#
# Two primitives, matching the vocabulary the rest of the suite already uses
# for the same job (`BPCells::multiply_rows` / `add_rows`, and ScaledMatrix's
# `scale` / `center`). The axis lives on the record, so no axis suffix here.
#
#   list(type = "multiply", axis = "cell"|"feat"|"all", factors = <payload>)
#   list(type = "add",      axis = "cell"|"feat"|"all", terms   = <payload>)
#
# `<payload>` is either a scalar (axis "all") or a named list mapping a
# substore's uid to a numeric vector INDEXED BY ON-DISK ID -- `factors[[uid]][id]`
# is the multiplier for that row_id / col_id. Same shape as the `@stats`
# marginals, and invariant under `[` for the same reason: on-disk ids do not
# move when a view narrows, so nothing has to be sliced or re-derived.
#
# NOT named "scale": at the workflow tier that word means standardize, centring
# included (`scaleParam`, the "scaled" expression slot, `ScaledMatrix` itself),
# while at the operation tier it means multiply only. `multiply` is unambiguous
# and leaves `add` free for its counterpart.
#
# The two are NOT interchangeable in where they can run:
#
#   multiply  preserves sparsity -- an implicit zero stays zero -- so it lowers
#             to Acero, applies over a collected triplet frame, and survives
#             any output mode.
#   add       destroys it -- every implicit zero becomes the offset -- so it
#             CANNOT be expressed over triplets at all. It is honored only when
#             a bounded chunk is materialized, by wrapping rather than by
#             mutating values (see .pe_add_wrap and the ScaledMatrix note on
#             .pe_check_dgc_dims).

.op_multiply <- function(atab, op) {
    value <- w <- NULL   # NSE
    axis <- op$axis %||% "cell"
    if (identical(axis, "all")) {
        k <- as.numeric(op$factors)
        return(dplyr::mutate(atab, value = value * !!k))
    }
    key <- if (identical(axis, "feat")) "col_id" else "row_id"
    tbl <- .pe_axis_payload_table(op$factors, key)
    tbl_a <- arrow::as_arrow_table(data.frame(
        source_id = as.character(tbl$source_id),
        key_id    = as.integer(tbl$key_id),
        w         = as.numeric(tbl$w),
        stringsAsFactors = FALSE
    ))
    by <- c("source_id" = "source_id"); by[key] <- "key_id"
    atab |>
        dplyr::left_join(tbl_a, by = by) |>
        dplyr::mutate(value = value * w) |>
        dplyr::select(-w)
}

# `add` is recorded but not yet executable. It cannot be lowered to Acero
# (arrow has no way to synthesize the implicit zeros), and the R-side executor
# would need to densify the slice into triplet form first -- see the op
# inventory above for the intended shape.
.op_add_refuse <- function(op) {
    stop("[.pe_apply_op] `add` ops are recorded but not yet executable. ",
         "Adding a per-", op$axis %||% "cell", " offset requires densifying ",
         "the slice (every implicit zero becomes the offset), which no ",
         "executor does yet.", call. = FALSE)
}


# ---- @post_ops R-side executor (data.table shape) --------------------------
#
# Used by materializing output paths (tibble / data.table / dgcmatrix) and
# by any consumer that collects a triplet chunk into a data.table. Ops
# mutate `df$value` in place.

.pe_apply_post_op_df <- function(df, op) {
    switch(op$type,
        "log"      = .op_transform_log(df, op),
        "multiply" = .pe_apply_post_op_multiply_df(df, op),
        "add"      = .op_add_refuse(op),
        stop("[.pe_apply_post_op_df] unknown post op type: ", op$type,
            call. = FALSE)
    )
}

.pe_apply_post_ops_df <- function(df, post_ops) {
    if (length(post_ops) == 0L) return(df)
    for (op in post_ops) df <- .pe_apply_post_op_df(df, op)
    df
}

.pe_apply_post_op_multiply_df <- function(df, op) {
    value <- source_id <- NULL   # NSE
    axis <- op$axis %||% "cell"
    if (identical(axis, "all")) {
        k <- as.numeric(op$factors)
        df[, value := value * k]
        return(df)
    }
    key <- if (identical(axis, "feat")) "col_id" else "row_id"

    # Positional index rather than a join: the payload is already a vector
    # keyed by on-disk id, so `w[id]` is the whole lookup. A union carries one
    # vector per substore, so split on source_id and index within each.
    #
    # `uniqueN` rather than `length(unique(...))`: the single-source branch
    # only needs the COUNT, and materializing every distinct string first
    # costs 20x more than counting them (0.061 s vs 0.003 s over 9.6M rows) --
    # which was two thirds of this executor's runtime.
    if (data.table::uniqueN(df$source_id) == 1L) {
        w <- .pe_axis_payload_vec(op$factors, df$source_id[1L], 0L)
        df[, value := value * w[get(key)]]
    } else {
        for (u in unique(df$source_id)) {
            w <- .pe_axis_payload_vec(op$factors, u, 0L)
            if (is.null(w)) next
            df[source_id == u, value := value * w[get(key)]]
        }
    }
    df
}

# ---- Chain edit helpers ----------------------------------------------------

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


# Read the store as it stands with NOTHING queued -- the values on disk,
# ignoring every op in either phase.
#
# For a consumer whose statistic is defined on the stored values rather than
# on whatever the chain produces. `.stream_filter_masks()` is the case: Giotto
# filter thresholds are count thresholds, and nothing stops a caller from
# filtering after normalizing, so the chain is suppressed explicitly rather
# than relied on to be absent.
#
# Note there is deliberately no "prefix up to op i" variant. Ops are
# positional and self-contained: a verb writes its record at the end of the
# chain and never revisits it, so no producer needs a basis cut partway
# through. An earlier version grew one to serve replace-in-place, which is
# exactly the reach-back that discipline forbids.
.pe_chain_none <- function(pe) {
    pe@ops <- list()
    pe@post_ops <- list()
    pe
}


# Demote an @ops entry, and every op after it, to the front of @post_ops.
#
# For a consumer that wants a lowerable op run R-side instead of in Acero:
# rather than routing around the chain (read `output = "query"` and reimplement
# the steps), it edits the chain and lets storeRead execute what it is given.
#
# The cascade is not optional. Materialization is one-way, so once op i runs
# R-side every op after it must too — this is the monotonic rule of
# .pe_push_op() expressed as a rewrite instead of a prohibition. Order within
# the chain is preserved end to end: the moved block keeps its internal order
# and lands ahead of whatever @post_ops already held, which by construction
# came after all of @ops.
#
# `from` selects the first op to demote, as either an @ops index or an op type
# (first match). A type with no match is an error — a consumer meaning "demote
# it if present" should guard with .pe_find_op_type().
.pe_demote_ops <- function(pe, from) {
    n <- length(pe@ops)
    if (is.character(from)) {
        idx <- .pe_find_op_type(pe@ops, from)
        if (is.na(idx)) {
            stop("[.pe_demote_ops] no op of type '", from, "' on @ops.",
                 call. = FALSE)
        }
        from <- idx
    }
    from <- as.integer(from)
    if (length(from) != 1L || is.na(from) || from < 1L || from > n) {
        stop("[.pe_demote_ops] `from` must select one of the ", n,
             " ops on @ops.", call. = FALSE)
    }

    moved <- pe@ops[from:n]
    pe@ops <- if (from == 1L) list() else pe@ops[seq_len(from - 1L)]
    pe@post_ops <- c(moved, pe@post_ops)
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
# `storeRead()` carries the same pre-materialization recipe restricted to this
# substore's rows. For arrow-native ops with source-keyed payload the
# per-substore filter is applied here. The norm ops now live on
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
# Injecting @post_ops here rather than at chunk-read time keeps the substore
# self-sufficient: `[` narrows its cell axis and `storeRead` applies the chain,
# with the payload carried through untouched (it is keyed by on-disk id).
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
