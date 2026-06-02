#' @include class-parquetExprStore.R
NULL

# parquetExprStore @ops — lazy arrow step recipes.
#
# Each op is a pure-data record `list(type = <character>, ...params)`.
# No closures, no fn field. At storeRead time, the executor translates the
# whole chain into composed arrow steps applied to the lazy query, then
# executes once via collect / output dispatch — matching the parquetGeomBase
# @ops contract.
#
# Op types currently supported:
#
#   norm_libsize_log
#     Fused library-size scaling + optional log transform. Adds `v_norm`
#     column to the lazy query.
#     Params:
#       scalef   data.table(source_id, orig_row_id, scalef). One row per
#                cell present in the store's current view at processData
#                time. The (source_id, orig_row_id) composite is required
#                so the same op record carries scalefs for all substores
#                of a unionParquetExprStore in a single table — row_id
#                restarts per substore so source_id is the disambiguator.
#                Single stores collapse to a constant source_id (= the
#                store's @uid).
#       log      logical. Apply log1p / log(base) after scaling.
#       base     numeric. log base (default 2).
#     Arrow translation:
#       atab |>
#         left_join(scalef_tab,
#                   by = c("source_id" = "source_id",
#                          "row_id"    = "orig_row_id")) |>
#         mutate(v_norm = value * scalef [|> log1p(.) * inv_log_base]) |>
#         select(-scalef)
#
# Extension protocol: add a new branch to .pe_do_op and a `slice_*`
# branch to .pe_slice_op_cells / .pe_slice_op_genes (commit 2). No other
# code paths change — the executor is the only place that knows about
# specific op types.


# ---- executor: translate one op record into arrow steps --------------------

# Apply a single op record to a lazy arrow query (or arrow_dplyr_query).
# Returns the augmented query.
.pe_do_op <- function(atab, op) {
    switch(op$type,
        "norm_libsize_log" = .pe_do_op_norm_libsize_log(atab, op),
        stop("[.pe_do_op] unknown op type: ", op$type, call. = FALSE)
    )
}

# Fold the chain — composes all ops into one lazy query.
.pe_apply_ops <- function(atab, ops) {
    if (length(ops) == 0L) return(atab)
    for (op in ops) atab <- .pe_do_op(atab, op)
    atab
}


# ---- per-type translators --------------------------------------------------

.pe_do_op_norm_libsize_log <- function(atab, op) {
    # NSE bindings
    row_id <- source_id <- value <- v_norm <- scalef <- NULL

    scalef_tab <- arrow::as_arrow_table(op$scalef)
    inv_log_base <- 1 / log(op$base)

    atab <- atab |>
        dplyr::left_join(scalef_tab,
                          by = c("source_id" = "source_id",
                                 "row_id"    = "orig_row_id")) |>
        dplyr::mutate(v_norm = value * scalef)
    if (isTRUE(op$log)) {
        atab <- dplyr::mutate(atab, v_norm = log1p(v_norm) * inv_log_base)
    }
    # Drop the joined scalef column so downstream consumers see only the
    # schema the recipe extended. v_norm sticks around as the recipe's
    # contribution.
    dplyr::select(atab, -scalef)
}


# ---- introspection helpers -------------------------------------------------

# Find the index of the first op of a given type in @ops. Returns NA_integer_
# when not present. Used by processData(logNormParam) to detect-and-fuse with
# an upstream libsize op rather than appending a redundant entry.
.pe_find_op_type <- function(ops, type) {
    if (length(ops) == 0L) return(NA_integer_)
    idx <- which(vapply(ops, function(op) identical(op$type, type),
        logical(1L)))
    if (length(idx) == 0L) NA_integer_ else as.integer(idx[1L])
}

# Predicate: does @ops carry a normalization recipe?
.pe_has_norm_op <- function(ops) {
    !is.na(.pe_find_op_type(ops, "norm_libsize_log"))
}


# ---- subset slice ----------------------------------------------------------
#
# Ops are frozen snapshots of population stats captured at trigger time.
# Subsetting cells / genes never invalidates them — every op carries axis-
# keyed lookup tables, and subset filters those tables by the surviving
# axis keys. Survivors retain their captured-population stats; the user
# re-runs processData(...) if they want stats over the new population.
#
# Each op type declares its lookup tables and axis keys via the registry
# below. Subset dispatch walks the registry and filters rows.

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
        # Semi-join the table on the axis key columns against surviving_keys.
        # data.table's `tbl[surviving_keys, on = .(...), nomatch = NULL]`
        # is an inner-join that yields all tbl columns; equivalent semi-join
        # by selecting tbl-side fields after the join.
        join_cols <- axis_spec  # named c(parquet_col = tbl_col) or simple char vec
        # For composite cell axis, surviving_keys columns are
        # (source_id, orig_row_id); axis_spec is c("source_id",
        # "orig_row_id") — same names on both sides.
        op[[tname]] <- tbl[surviving_keys, on = join_cols,
                            nomatch = NULL, allow.cartesian = FALSE]
        # Drop any join-introduced columns from surviving_keys side if they
        # differ from tbl columns (none expected for the current schemas).
    }
    op
}


# Slice every op in a chain along one axis.
.pe_slice_ops <- function(ops, axis, surviving_keys) {
    if (length(ops) == 0L) return(ops)
    lapply(ops, .pe_slice_op, axis = axis, surviving_keys = surviving_keys)
}
