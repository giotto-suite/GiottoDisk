# Kernel: scan construction -- axis predicates, the composed scan, the id remap.
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

# ---- shared axis predicates -------------------------------------------------
#
# `cell_idx` / `gene_idx` narrowing becomes an arrow predicate. Three shapes,
# all exact -- the predicate admits precisely the in-view entries, so nothing
# downstream re-filters:
#
#   1. gapless               -> range alone
#   2. gaps, few dropped     -> range AND `!(x %in% dropped)`
#                               (the range is required for CORRECTNESS here)
#   3. gaps, few kept        -> `x %in% kept` alone, NO range
#
# Case 3 looks like an omission and is not: a `col_id` range prunes no row
# groups (the file is cell-major) while still costing a comparison per row.
# Adding it "for symmetry" is a measured regression. adr/0008 has the numbers
# and the sort-order argument.
#
# Bounds come from min/max, never first/last -- `idx` is not guaranteed sorted
# (`feats_to_use` may be HVG-rank ordered). Gap detection runs on unique values
# so duplicates cannot make `n == span` accidentally true, and `dropped` is
# materialized only when case 2 wins.

#' @keywords internal
#' @noRd
.pe_axis_pred <- function(idx) {
    if (length(idx) == 0L) return(NULL)
    u  <- base::unique(as.integer(idx))
    lo <- min(u)
    hi <- max(u)
    span <- hi - lo + 1L
    gapless  <- (length(u) == span)
    use_anti <- !gapless && (span - length(u)) < length(u)
    list(lo = lo, hi = hi,
         gapless   = gapless,
         use_anti  = use_anti,
         use_range = gapless || use_anti,
         dropped   = if (use_anti) setdiff(lo:hi, u) else integer(0),
         kept      = u)
}

# Emit a plan as quoted predicate expressions over column `col`, with values
# inlined as literals. One emitter for both consumers: `storeRead` chains them
# onto a lazy query, `.union_substore_filter_expr()` ANDs them into a
# per-substore clause. Returns list() when there is nothing to filter.

#' @keywords internal
#' @noRd
.pe_axis_pred_exprs <- function(plan, col) {
    if (is.null(plan)) return(list())
    sym <- as.name(col)
    out <- list()
    if (plan$use_range) {
        out[[length(out) + 1L]] <-
            bquote(.(sym) >= .(plan$lo) & .(sym) <= .(plan$hi))
    }
    if (!plan$gapless) {
        out[[length(out) + 1L]] <- if (plan$use_anti) {
            bquote(!(.(sym) %in% .(plan$dropped)))
        } else {
            bquote(.(sym) %in% .(plan$kept))
        }
    }
    out
}

# Apply an axis plan to a lazy query, on either carrier. One site for both, so
# a subset can never be applied to one output and skipped on the other.
#
# The shapes come from `.pe_axis_pred_exprs()` in every case; this only decides
# HOW a membership set reaches the engine. Acero takes it as a hash set, which
# is what adr/0008 measured and there is nothing to fix. dbplyr inlines it into
# the query TEXT, so on a tbl_dbi a large set becomes a large SQL string:
# measured over a 500k-row scan, 1k ids cost the same as a registered
# semi-join, 20k cost 8x, and 100k cost 38x on top of 778 KB of SQL. Above the
# threshold the ids are registered instead and the test becomes a semi/anti
# join, which dbplyr renders as EXISTS / NOT EXISTS -- the same shape the
# tabular path's `id_filter` already uses, and NULL-safe where NOT IN is not.
#
# The range half is kept as literals either way: shape 2 needs it for
# correctness (adr/0008), and it is also the only half DuckDB can turn into a
# row-group prune.

#' @keywords internal
#' @noRd
.pe_apply_axis_pred <- function(x, plan, col) {
    if (is.null(plan)) return(x)
    ids <- if (plan$use_anti) plan$dropped else if (!plan$gapless) plan$kept
           else integer(0L)
    thresh <- getOption("giottodisk.duckdb_in_subquery_threshold", 1000L)
    if (!inherits(x, "tbl_dbi") || length(ids) <= thresh) {
        for (p in .pe_axis_pred_exprs(plan, col)) x <- dplyr::filter(x, !!p)
        return(x)
    }
    if (plan$use_range) {
        sym <- as.name(col)
        x <- dplyr::filter(x,
            !!bquote(.(sym) >= .(plan$lo) & .(sym) <= .(plan$hi)))
    }
    ids_tbl <- .pe_register_ids(dbplyr::remote_con(x), ids, col)
    if (plan$use_anti) {
        dplyr::anti_join(x, ids_tbl, by = col)
    } else {
        dplyr::semi_join(x, ids_tbl, by = col)
    }
}

#' @keywords internal
#' @noRd
.pe_register_ids <- function(conn, ids, col) {
    name <- tolower(paste0("gd_peid_", .make_uid()))
    tab <- do.call(arrow::arrow_table,
        stats::setNames(list(as.integer(ids)), col))
    duckdb::duckdb_register_arrow(conn, name, tab)
    dplyr::tbl(conn, name)
}

# A store's lazy scan from its base scan: both axis plans, then the @ops
# prefix. storeRead() folds this into @read_fun; a fanned-out window runs it in
# a worker from the same plans, lowered by `.pe_lower_read()`.
.pe_compose_scan <- function(ds, ci_plan, gi_plan, ops) {
    ds <- .pe_apply_axis_pred(ds, ci_plan, "row_id")
    ds <- .pe_apply_axis_pred(ds, gi_plan, "col_id")
    if (length(ops) > 0L) ds <- .pe_apply_ops(ds, ops)
    ds
}

# The default parquetExprStore @read_fun. Named so `.pe_lower_read()` can tell
# a store that reads its path as a plain dataset from one given its own reader.
.pe_read_dataset <- function(x, ...) arrow::open_dataset(sources = x, ...)

# Join the lookup tables from `.pestore_remap_luts()` onto a scan, renumbering
# row_id / col_id to positions in the view, and sort cell-major.
.pestore_apply_remap <- function(q, luts) {
    row_id <- col_id <- value <- row_id_new <- col_id_new <- NULL  # NSE
    cell_remap <- arrow::as_arrow_table(luts$cell)
    # `by` must be fully named -- arrow's dplyr join handler trips on the
    # mixed-named form `c("source_id", "row_id" = "row_id_orig")` because
    # the unnamed element parses with an empty name on the right side.
    if (isTRUE(luts$by_source)) {
        q <- dplyr::left_join(q, cell_remap,
            by = c("source_id" = "source_id",
                   "row_id" = "row_id_orig"))
    } else {
        q <- dplyr::left_join(q, cell_remap,
            by = c("row_id" = "row_id_orig"))
    }
    q <- dplyr::mutate(q, row_id = row_id_new)
    q <- dplyr::select(q, -dplyr::any_of(c("row_id_new", "source_id")))
    if (!is.null(luts$gene)) {
        q <- dplyr::left_join(q, arrow::as_arrow_table(luts$gene),
            by = c("col_id" = "col_id_orig"))
        q <- dplyr::mutate(q, col_id = col_id_new)
        q <- dplyr::select(q, -col_id_new)
    }
    q <- dplyr::arrange(q, row_id, col_id)
    dplyr::select(q, row_id, col_id, value)
}

# Unique ids for artifacts, and for the tables a duckdb carrier registers.
.make_uid <- function(n = 8L, include_node = NULL, include_pid = NULL) {
    include_node <- include_node %||% 
        getOption("giottodisk.uid_include_node", FALSE)
    include_pid <- include_pid %||% 
        getOption("giottodisk.uid_include_pid", TRUE)
    count <- getOption("giottodisk.uid_count", 1L)
    on.exit({
        options("giottodisk.uid_count" = count + 1L)
    })
    GiottoUtils::gwith_seed(seed = Sys.time() + count, {
        sampleset <- sample(c(letters, LETTERS, as.character(0:9)), 
            size = n,
            replace = TRUE
        )
        rand <- paste(sampleset, collapse = "")
    })

    parts <- rand
    if (include_pid) {
        pid <- Sys.getpid()
        parts <- c(pid, parts)
    }
    if (include_node) {
        host <- substr(Sys.info()[["nodename"]], 1, 8)
        parts <- c(host, parts)
    }
    paste0(parts, collapse = "_")
}
