# Kernel: the op chain interpreter -- executors for @ops and @post_ops records.
#
# R/kernel-*.R is the part of GiottoDisk that runs where GiottoDisk is not
# loaded. `.kernel_bundle()` (R/utils-isolate.R) copies these functions out of
# the namespace for a worker that has only arrow, dplyr and data.table, so each
# one takes plain data (op records, axis plans, lookup tables, paths) plus a
# lazy query or a data.table, and:
#   * calls other packages only as `pkg::fn`;
#   * calls no GiottoDisk function outside R/kernel-*.R;
#   * reaches no S4 generic or class -- including a base name GiottoDisk masks
#     with a generic, such as `unique`, which is why the kernel writes
#     `base::unique`;
#   * reads options only with a default (the worker gets a giottodisk.*
#     snapshot, not the parent's session).
# test-kernel-isolation.R enforces the list. See AGENTS.md, "Kernel".

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
#   compare           (phase: lazy or post)
#     An indicator: keep the stored entries where `value <op> e2` holds and
#     set them to 1. Emitted by the `Compare` group methods, so `x >= 1` is
#     this record and `rowSums(x >= 1)` counts through the accumulator.
#     Unlike every op above it CHANGES WHICH ENTRIES ARE STORED -- it is the
#     first record that drops rows -- but it keeps the absent-means-0
#     invariant, because the producer refuses any comparison that is TRUE at 0
#     (that result would be dense). 1 rather than TRUE because `value` is a
#     double column on every carrier.
#     Params:
#       op       one of ">", ">=", "<", "<=", "==", "!=", store on the left
#       axis     "all" only, today. Shaped like `multiply` so a per-axis
#                threshold can be added without changing the record: base
#                R's `m >= v` with length(v) == nrow(m) is a per-row
#                threshold, which here would be `axis = "feat"` (or
#                "cell" when transposed) with `e2` a payload keyed by
#                on-disk id. Comparing against another store is a join,
#                not this record.
#       e2       numeric scalar (axis "all")
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
        "expm1"    = .op_transform_expm1(atab, op),
        "multiply" = .op_multiply(atab, op),
        "compare"  = .op_compare(atab, op),
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

# The inverse of `log`: value -> base^value - 1.
#
# Sparsity-preserving, which is the only reason it can live here at all: an
# absent entry is 0 and base^0 - 1 is 0, so implicit zeros stay implicit and
# the op lowers to Acero like `multiply` rather than being refused like `add`.
#
# Written as `base^value - 1` rather than `expm1(value * log(base))` because
# the in-memory backends compute the former and the two are not bit-identical
# in general. Arrow's `^` agrees with R's elementwise, which is what makes the
# streamed and in-memory PAGE means comparable rather than merely close.
.op_transform_expm1 <- function(x, op) {
    value <- NULL # NSE
    base <- op$base %||% 2
    if (data.table::is.data.table(x)) {
        return(x[, value := base^value - 1])
    }
    dplyr::mutate(x, value = (!!base)^value - 1)
}

# Indicator of a comparison against a scalar. One predicate call serves every
# carrier: spliced into `filter()` it lowers to Acero and dbplyr, and evaluated
# over the column it is the data.table row selector. `which()` drops NA the way
# `filter()` does, so the two phases agree on a stored NA as well.
#
# Returns a new table rather than mutating in place, since it drops rows --
# every caller already reassigns `df <- .pe_apply_post_ops_df(df, ...)`.
.op_compare <- function(x, op) {
    value <- NULL # NSE
    if (!identical(op$axis %||% "all", "all")) {
        stop("[.op_compare] only `axis = \"all\"` (a scalar threshold) is ",
             "implemented; got axis \"", op$axis, "\".", call. = FALSE)
    }
    pred <- call(op$op, quote(value), op$e2)
    if (data.table::is.data.table(x)) {
        x <- x[which(eval(pred, x))]
        return(x[, value := 1])
    }
    x |>
        dplyr::filter(!!pred) |>
        dplyr::mutate(value = 1)
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
        "expm1"    = .op_transform_expm1(df, op),
        "multiply" = .pe_apply_post_op_multiply_df(df, op),
        "compare"  = .op_compare(df, op),
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
        for (u in base::unique(df$source_id)) {
            w <- .pe_axis_payload_vec(op$factors, u, 0L)
            if (is.null(w)) next
            df[source_id == u, value := value * w[get(key)]]
        }
    }
    df
}


# ---- op payloads -----------------------------------------------------------
#
# Same indexing convention as the `@stats` marginals above: a numeric vector
# whose POSITION is the on-disk row_id / col_id. Kept together because that
# convention is the thing to preserve -- it is what makes both invariant under
# `[`, and what lets the R-side executor index directly where arrow has to
# join.

# Resolve a payload to a full-length numeric vector for one substore.
.pe_axis_payload_vec <- function(payload, uid, n) {
    if (is.null(payload)) return(NULL)
    if (!is.list(payload)) return(rep_len(as.numeric(payload), n))
    v <- payload[[as.character(uid)]]
    if (is.null(v)) return(NULL)
    as.numeric(v)
}

# Build the joinable (source_id, <key>, w) table an Acero plan needs. Arrow
# cannot index an R vector from inside a query, so per-axis state has to arrive
# as a table -- rebuilt per call, which is what the previous executor did too.
.pe_axis_payload_table <- function(payload, key) {
    src <- names(payload)
    data.table::rbindlist(lapply(src, function(u) {
        v <- as.numeric(payload[[u]])
        data.table::data.table(
            source_id = rep_len(as.character(u), length(v)),
            key_id    = seq_along(v),
            w         = v
        )
    }))[!is.na(w)]
}
