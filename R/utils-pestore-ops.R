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
# Why the slot means position rather than capability: adr/0005 (which
# supersedes adr/0002, the original phase split and its measurements). Why payloads
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
#     Note for whoever implements it: every sufficient-statistic verb here
#     (`.pe_accum_raw` and its callers -- QC stats, grouped feature stats, HVF,
#     marker moments) is correct only because an absent entry means 0 in the
#     value space being aggregated. `multiply` has f(0) = 0 and `log` is log1p
#     with `offset != 1` refused precisely to keep it. `add` is the first op
#     that would break it, and it would break it SILENTLY -- wrong means, no
#     error. The roadmap's densification into triplet form is what preserves
#     the invariant, by making the zero block explicit before f(0) != 0 can
#     matter. It is not an optimisation detail; ship it with the op.
#
#   log               (phase: lazy or post)
#     log(value + 1) / log(base). Carries no axis-keyed state.
#     Params:
#       base     numeric. log base (default 2).
#
# Records are positional and self-contained: each does its work at the
# position it occupies, and a verb appends rather than revisiting anything it
# wrote earlier. Nothing needs to be applied in a particular order or to be
# present at all -- the chain supplies the sequencing.
#
# Extension protocol:
#   - Add a branch to .pe_apply_op (lazy) and/or .pe_apply_post_op_df
#     (triplets), depending on which engines can express it.
#   - .pe_apply_op serves BOTH lazy carriers -- Acero and a DuckDB tbl_dbi --
#     so write the branch in plain dplyr and it lowers to both. Reach for a
#     carrier test only where an engine cannot accept the other's data, as
#     .pe_payload_carrier does; a branch per engine is how the two outputs
#     drift apart.
#   - Key any payload by ON-DISK id, not by view position. That is what makes
#     it invariant under `[` -- see the subset-slice note at the bottom.
#   - Have the producing verb append the record; never edit an existing one.


# ---- @ops lazy-side executor ----------------------------------------------
#
# Carrier-agnostic: the same fold runs over an Acero query and over a DuckDB
# `tbl_dbi`, because every record here has a dplyr form and dplyr targets both.
# That is what keeps `output = "query"` and `output = "duckdb"` returning the
# same values -- the equivalence is structural, not something the tests police.
# Only `.op_multiply`'s payload has to know which engine it landed on.

# Apply a single prefix op record to the lazy query.
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

# `log(value + 1)` rather than `log1p(value)`: DuckDB has no log1p and dbplyr
# has no translation for it, so log1p would reach the engine verbatim and fail
# at collect. log() is native to Acero, data.table and dbplyr (-> LN) alike,
# which keeps this one expression across all three carriers rather than a
# branch per engine. Cost on Acero is nil -- measured slightly faster
# single-threaded and indistinguishable at the default thread count, where the
# transform is memory-bandwidth bound. The accuracy difference is confined to
# value << 1 and is ~1e-16 absolute at value = 1e-9, orders below anything a
# library-normalized count reaches.
.op_transform_log <- function(x, op) {
    value <- NULL # NSE
    base <- op$base %||% 2
    if (data.table::is.data.table(x)) {
        return(x[, value := log(value + 1) / log(base)])
    }
    dplyr::mutate(x, value = log(value + 1) / log(!!base))
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
    tbl_a <- .pe_payload_carrier(atab, op$factors, key)
    by <- c("source_id" = "source_id"); by[key] <- "key_id"
    atab |>
        dplyr::left_join(tbl_a, by = by) |>
        dplyr::mutate(value = value * w) |>
        dplyr::select(-w)
}

# Put the payload on the same engine as the query it joins into. Acero cannot
# read a tbl_dbi and DuckDB cannot read an arrow Table, so this is the one
# place in the chain where the carrier matters -- the join / multiply / drop
# above stays a single expression for both.
#
# The carrier is read off `x` rather than passed in. That matches
# `.op_transform_log`'s existing shape and keeps `.pe_apply_ops(x, ops)`
# callable unchanged from every site that already calls it.
#
# duckdb_register_arrow, not dbplyr::copy_inline: copy_inline writes the
# payload into the query TEXT as a literal VALUES list (one row per cell or
# per feature) and casts `w` to NUMERIC, which would move `value` off DOUBLE.
# Registration is zero-copy and leaves the types alone -- `key_id` stays int32
# against the int32 row_id/col_id from parquet, so the join needs no cast.
#
# No COALESCE on `w`: an unmatched key yields NA here, from arrow's
# `value * NA` and from DuckDB's NULL alike, matching the out-of-range index
# in .pe_apply_post_op_multiply_df. Defaulting the factor to 1 would instead
# return that entry's RAW value dressed as a normalized one.
.pe_payload_carrier <- function(x, factors, key) {
    tbl <- .pe_axis_payload_table(factors, key)
    tab <- arrow::as_arrow_table(data.frame(
        source_id = as.character(tbl$source_id),
        key_id    = as.integer(tbl$key_id),
        w         = as.numeric(tbl$w),
        stringsAsFactors = FALSE
    ))
    if (!inherits(x, "tbl_dbi")) return(tab)
    name <- tolower(paste0("gd_pew_", .make_uid()))
    duckdb::duckdb_register_arrow(dbplyr::remote_con(x), name, tab)
    dplyr::tbl(dbplyr::remote_con(x), name)
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
