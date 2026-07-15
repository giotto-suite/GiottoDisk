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
#       scalef   data.table(orig_row_id, scalef). One row per cell present
#                in the store's current view at processData time.
#                orig_row_id is the parquet's on-disk row_id (cell index
#                in the original parquet, surviving subset-slice cleanly
#                via row filtering — see .pe_slice_op_cells).
#       log      logical. Apply log1p / log(base) after scaling.
#       base     numeric. log base (default 2).
#     Arrow translation:
#       atab |>
#         left_join(scalef_tab, by = c("row_id" = "orig_row_id")) |>
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
    row_id <- value <- v_norm <- scalef <- NULL

    scalef_tab <- arrow::as_arrow_table(op$scalef)
    inv_log_base <- 1 / log(op$base)

    atab <- atab |>
        dplyr::left_join(scalef_tab, by = c("row_id" = "orig_row_id")) |>
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
