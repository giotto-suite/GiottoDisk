#' @include class-parquetExprStore.R
NULL

# The op chain model, the op registry and the extension protocol are the
# header of R/kernel-ops.R, beside the executors they describe. This file
# holds the chain editors and the window walk, which need the store.

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


# ---- cell windowing --------------------------------------------------------
#
# THE seam for streaming an expression store a cell window at a time. Every
# bounded pass over expression values goes through here: the statistic
# accumulators, the PCA forward/backward/gram/coords passes, and the
# `storeWrite` bake. Do not hand-roll a `while (cs <= n)` walk -- that loop
# existed in eight places before these two helpers, and the copies had already
# drifted (one folded its partials eagerly, one retained one per window).
#
# Why cell ranges and not row counts, or a cell SET: a contiguous `row_id`
# range is the gapless case in `.pe_axis_pred()`, so it lowers to a pure
# `row_id >= lo & row_id <= hi` that prunes parquet row groups -- the store is
# written cell-major (`setorder(row_id, col_id)`), so this is the only axis
# where narrowing prunes. A scattered set lowers to `is_in` and reads
# everything. It is also why windowing the FEATURE axis is a trap: it prunes
# nothing, so each batch rescans the store in full.
#
# Two shapes, because the call sites genuinely differ:
#
#   .pe_chunk_ranges()  the primitive -- chunk boundaries within one range.
#                       For a caller that already has its substore and a
#                       sub-range of it, as the parallel PCA band workers do.
#   .pe_windows()       substores x their full cell range, as descriptors.
#                       For a caller that means "the whole view".
#
# Neither owns reduction, and that is deliberate: the four reductions in the
# package are irreconcilable (eager fold, scatter into a preallocated matrix,
# matrix accumulation, write a part-file). An iterator that tried to own them
# would grow a mode argument per caller.

# Chunk boundaries covering `[from, to]`. Returns a list of `c(cs, ce)`.
.pe_chunk_ranges <- function(from, to, chunk_size) {
    from <- as.integer(from)
    to   <- as.integer(to)
    if (is.na(from) || is.na(to) || to < from) return(list())
    chunk_size <- max(1L, as.integer(chunk_size))
    starts <- seq.int(from, to, by = chunk_size)
    lapply(starts, function(cs) c(cs, min(cs + chunk_size - 1L, to)))
}

# Cell-window descriptors over a whole store or union. Each is
#
#   list(sub = <parquetExprStore>, cs = , ce = , offset = , index = )
#
# `sub` carries both op chains already (a union parent's are transplanted by
# `.exprbase_inject_parent_ops`, which is what makes the substore
# self-sufficient). `cs`/`ce` are positions in THAT substore's cell axis, not
# the view's -- `offset` converts, and is what the write path and the PCA passes
# use to place a window's rows globally. `index` is the substore's ordinal, for
# a caller carrying its own per-substore struct to look up (PCA's `sub_infos`).
#
# Unions iterate substores rather than taking a global cell range because
# `row_id` restarts per substore, so a global range is not a contiguous range
# on either side of the boundary and would prune nothing.
#
# `inject_ops = FALSE` skips the parent-op transplant, for a caller that has
# already prepared its substores.
.pe_windows <- function(pe, chunk_size, inject_ops = TRUE) {
    is_union <- inherits(pe, "unionParquetExprStore")
    out <- list()
    subs <- .exprbase_substores(pe)
    for (i in seq_along(subs)) {
        sub <- subs[[i]]$store
        # Only a union needs the transplant: for a single store
        # `.exprbase_substores()` yields the store itself, which already
        # carries its chains, and re-adding them would double-apply.
        if (is_union && isTRUE(inject_ops)) {
            sub <- .exprbase_inject_parent_ops(sub, pe@ops, pe@post_ops)
        }
        n_sub <- as.integer(sub@n_cells)
        for (rng in .pe_chunk_ranges(1L, n_sub, chunk_size)) {
            out[[length(out) + 1L]] <- list(
                sub    = sub,
                cs     = rng[[1L]],
                ce     = rng[[2L]],
                offset = as.integer(subs[[i]]$cell_offset),
                index  = i
            )
        }
    }
    out
}

# The window as a store to read from. Skips `[` when the window covers the
# whole substore: the narrowing would add an exact-range predicate that admits
# every row anyway, and a store with no `@cell_idx` is the cheaper plan.
#
# Index with `cs:ce`, never `seq.int(cs, ce)` -- an ALTREP compact seq becomes
# one hyperslab where a materialized vector becomes a point selection.
.pe_window_store <- function(d) {
    if (d$cs == 1L && d$ce == as.integer(d$sub@n_cells)) return(d$sub)
    d$sub[, d$cs:d$ce]
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


# Orientation ####

# Swap a (row, col) pair iff the store is transposed. Every site that turns a
# logical position into a semantic axis goes through this, so the flip is
# spelled once and the op machinery never has to know about it: ops name
# semantic axes bound to on-disk columns, so orientation is invisible to them.
#' @keywords internal
#' @noRd
.pe_orient <- function(x, row_side, col_side) {
    if (x@transposed) list(col_side, row_side) else list(row_side, col_side)
}
