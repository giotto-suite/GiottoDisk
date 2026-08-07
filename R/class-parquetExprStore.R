#' @include class-dataStore.R
NULL

# docs ####

#' @name parquetExprStore-class
#' @title Parquet Expression Matrix Store (streaming)
#' @description
#' S4 class for **disk-backed expression matrices** stored as long-format
#' Apache Parquet, designed for streaming read access. Unlike
#' [parquetStore-class] (which stores arbitrary tabular data with a
#' `row_index` column), `parquetExprStore` represents a sparse cell x gene
#' expression matrix using three integer / float columns:
#'
#' | Column   | Type    | Meaning              |
#' |----------|---------|----------------------|
#' | `row_id` | int32   | 1-based cell index   |
#' | `col_id` | int32   | 1-based gene index   |
#' | `value`  | float64 | expression count     |
#'
#' Parquet files are sorted by `row_id` so that Arrow predicate-pushdown
#' row-group skipping makes chunked reads fast on large datasets.
#'
#' Cell barcodes (`cell_IDs`) and gene names (`feat_IDs`) are stored as
#' character slots on the S4 object -- the Parquet payload itself stays
#' minimal. The slot vectors act as a lookup table from integer index to
#' character ID.
#'
#' @section Use case:
#' This store is the streaming expression backend for datasets too large to
#' hold in memory as a sparse matrix. It is slotted into `exprObj@exprMat`
#' like any other disk-backed store; downstream Giotto methods that recognize
#' the class dispatch to streaming implementations, which read the triplet
#' payload in cell chunks rather than materializing the whole matrix.
#'
#' @slot path character. Local file path (single Parquet) or directory
#'   (one Parquet per chunk for very large datasets). Arrow's `open_dataset`
#'   handles both transparently.
#' @slot uid character. Auto-generated unique ID for artifact tracking.
#' @slot read_fun function. Preset to `arrow::open_dataset()`.
#' @slot params list. Reserved for downstream pipeline metadata
#'   (e.g. HVG indices after `sc_hvg`). Not used for normalization recipes
#'   anymore — those live on `@ops`.
#' @slot ops list. The part of the op chain that runs **before**
#'   materialization. `@ops` and `@post_ops` are one ordered sequence split at
#'   the point where execution leaves Acero — this is the prefix, not a
#'   collection of whichever steps happen to be lowerable.
#'
#'   Each entry is a pure-data `list(type, ...params)` record (no closures),
#'   so the chain survives `saveRDS` cleanly. At `storeRead()` time
#'   `.pe_apply_op()` translates each into arrow-dplyr steps on the lazy
#'   query, and the whole prefix compiles into one plan executed once at
#'   collect. Every record here must therefore lower to arrow — a consequence
#'   of the position rather than the slot's definition.
#'
#'   Independent of `@cell_idx` / `@gene_idx`: `[` narrows the window and never
#'   touches the chain, and the chain never consults the window at read time. An
#'   op whose meaning depends on the window (library normalization, whose
#'   factors come from column sums over the features in view) freezes that
#'   statistic into its payload when the producing verb runs. Re-running the
#'   verb is how you ask for a statistic over a new population; subsetting is
#'   not. See `adr/0006`.
#'
#'   Empty by default; populated by `processData()` methods. See
#'   `R/utils-pestore-ops.R` for the op type registry.
#' @slot n_cells numeric. Number of cells in the dataset
#'   (length of `cell_ids`).
#' @slot n_genes numeric. Number of genes / features
#'   (length of `feat_ids`).
#' @slot cell_ids character. Cell barcodes; index `i` corresponds to
#'   `row_id == i` in the Parquet file.
#' @slot feat_ids character. Gene / feature IDs; index `j` corresponds to
#'   `col_id == j` in the Parquet file.
#' @slot stats list. Cached marginal counts for the Parquet payload, filled
#'   by `storeWrite()`. Two integer vectors:
#'
#'   \describe{
#'     \item{`col_nnz`}{stored-entry count per feature, indexed by on-disk
#'       `col_id`.}
#'     \item{`row_nnz`}{stored-entry count per cell, indexed by on-disk
#'       `row_id`.}
#'   }
#'
#'   Keyed by **on-disk id**, not by identifier name and not by view position,
#'   which is what makes them invariant under `[`: subsetting only narrows
#'   `@cell_idx` / `@gene_idx` against the same file, so
#'   `sum(stats$col_nnz[gene_idx])` is the exact nonzero count of a
#'   gene-narrowed view. Names would break under feature renaming; view
#'   positions would be invalidated by every subset.
#'
#'   Their lengths are the *file's* dimensions, which after a subset is the
#'   only place that survives — `@n_genes` / `@n_cells` describe the view.
#'
#'   A subset on one axis is exact; a subset on both scales the exact axis by
#'   the other's kept fraction, which assumes nonzeros spread uniformly across
#'   the cell axis. Empty for a store whose Parquet was not written through
#'   `storeWrite()`; consumers fall back to counting.
#'
#'   Note `storeWrite()` **renumbers** ids when the input is subset, so a
#'   written store's marginals are computed against the new file rather than
#'   inherited from its parent.
#' @family store types
#' @seealso [parquetExprStore()], [mtxInput()]
NULL

# definitions ####

# Virtual base for any cell x gene expression store backed by Parquet.
# Used as a dispatch root so streaming pipeline methods (filter, normalize,
# HVF, QC) can be written once for both the single-file `parquetExprStore`
# and the multi-substore `unionParquetExprStore`. Per-substore iteration
# is exposed via the `.exprbase_substores()` protocol (see
# `R/utils-pestore-ops.R`); methods loop over substores, aggregate, and
# resolve `(source, row_id)` or `col_id` back to union cell / feat
# positions.
#
# Concrete subclasses (`parquetExprStore`, `unionParquetExprStore`) keep
# their own slot declarations; this is a tag-only virtual to avoid any
# slot-relocation churn or RDS deserialization risk.
#' @rdname parquetExprStore-class
#' @exportClass parquetExprBase
setClass("parquetExprBase",
    contains = "VIRTUAL"
)

#' @rdname parquetExprStore-class
#' @section Subset semantics:
#' `parquetExprStore` supports lightweight subsetting via the `[` operator.
#' `pe[i, j]` returns a new store whose `feat_ids` / `cell_ids` are narrowed
#' to the kept rows / columns and whose `gene_idx` / `cell_idx` slots record
#' the *original* parquet positions of the kept entries. The Parquet file
#' on disk is **not** rewritten -- `storeRead()` filters lazily via Arrow
#' using the recorded indices, so chained subsets stay cheap.
#' @slot cell_idx integer. Active cell positions in the original Parquet
#'   (length 0 = no subset, all cells are active). When non-empty:
#'   `length(cell_idx) == length(cell_ids) == n_cells`. This is *view* state and
#'   is independent of the op chain — narrowing it never invalidates a queued
#'   op, because window-dependent ops froze their statistic at push time
#'   (`adr/0006`).
#' @slot gene_idx integer. Active gene positions in the original Parquet
#'   (length 0 = no subset).
setClass("parquetExprStore",
    contains = c("queryableStore", "parquetExprBase"),
    slots = list(
        n_cells    = "numeric",
        n_genes    = "numeric",
        cell_ids   = "character",
        feat_ids   = "character",
        cell_idx   = "integer",
        gene_idx   = "integer",
        stats      = "list",   # cached on-disk marginals (see @stats)
        ops        = "list",   # chain prefix, run in Acero
        post_ops   = "list"    # chain suffix, run R-side after collect
    ),
    prototype = list(
        n_cells    = 0,
        n_genes    = 0,
        cell_ids   = character(0L),
        feat_ids   = character(0L),
        cell_idx   = integer(0L),
        gene_idx   = integer(0L),
        stats      = list(),
        ops        = list(),
        post_ops   = list()
    )
)

# constructor ####

#' @name parquetExprStore
#' @title Create a Parquet Expression Matrix Store
#' @description
#' Construct a [parquetExprStore-class] handle around a long-format Parquet
#' file or directory.  The Parquet must contain `row_id`, `col_id`, `value`
#' columns and be sorted by `row_id`.
#'
#' To populate the Parquet from a 10x / Xenium MatrixMarket triple, build
#' an [mtxInput()] and pass it to [storeWrite()].
#'
#' @param path character. Path to a Parquet file or directory of Parquet
#'   chunks.
#' @param cell_ids character. Cell barcodes. Length must equal `n_cells`.
#' @param feat_ids character. Gene / feature IDs. Length must equal
#'   `n_genes`.
#' @param n_cells numeric. Total number of cells. Defaults to
#'   `length(cell_ids)`.
#' @param n_genes numeric. Total number of genes. Defaults to
#'   `length(feat_ids)`.
#' @param scan_stats logical. Scan the Parquet at `path` to cache its marginal
#'   nonzero counts on `@stats` (default `FALSE`). Only meaningful when `path`
#'   already holds data: the usual pattern is to construct an empty handle and
#'   populate it with [storeWrite()], which caches the marginals itself. Set
#'   `TRUE` when attaching to a Parquet written elsewhere and you would rather
#'   pay the scan now than have the first consumer pay it. Leaving it `FALSE`
#'   costs correctness nothing — consumers that need the counts fall back to
#'   counting on demand — so this is purely about when the scan happens.
#'   Reachable through [storeCreate()], which forwards `...` here.
#' @param ... additional slots passed through to `new()`.
#' @return A `parquetExprStore` S4 object.
#' @family store constructors
#' @export
parquetExprStore <- function(
    path       = .dump_tempfile(),
    cell_ids   = character(0L),
    feat_ids   = character(0L),
    n_cells    = length(cell_ids),
    n_genes    = length(feat_ids),
    scan_stats = FALSE,
    ...
) {
    if (length(cell_ids) > 0L && length(cell_ids) != n_cells) {
        stop("[parquetExprStore] length(cell_ids) (", length(cell_ids),
             ") != n_cells (", n_cells, ").", call. = FALSE)
    }
    if (length(feat_ids) > 0L && length(feat_ids) != n_genes) {
        stop("[parquetExprStore] length(feat_ids) (", length(feat_ids),
             ") != n_genes (", n_genes, ").", call. = FALSE)
    }
    store <- new("parquetExprStore",
        path       = path,
        n_cells    = as.numeric(n_cells),
        n_genes    = as.numeric(n_genes),
        cell_ids   = as.character(cell_ids),
        feat_ids   = as.character(feat_ids),
        ...
    )
    if (isTRUE(scan_stats)) store <- .pestore_finalize_stats(store)
    store
}

# initialize ####

setMethod("initialize", signature("parquetExprStore"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    # default reader is arrow::open_dataset (handles file or directory)
    if (.is_empty_fun(.Object@read_fun)) {
        .Object@read_fun <- function(x, ...) arrow::open_dataset(sources = x, ...)
    }
    .Object
})


# unionParquetExprStore ####

#' @name unionParquetExprStore-class
#' @title Virtual Union of Expression Stores
#' @description
#' Lazy column-wise (cell-wise) concatenation of N
#' [parquetExprStore-class] objects. Substores must share an identical
#' `feat_ids` vector (same panel, same ordering); `cell_ids` accumulate
#' across substores and must be globally unique (caller pre-prefixes if
#' needed).
#'
#' Construction is O(1) — the union is purely virtual via Arrow's
#' `UnionDataset` over the substores' on-disk hive-partitioned datasets.
#' No data is rewritten; substores remain independently usable.
#' @slot stores list of [parquetExprStore-class] objects
#' @slot cell_ids character. Concatenated cell barcodes.
#' @slot feat_ids character. Shared feature IDs.
#' @slot n_cells numeric. Sum of substore n_cells.
#' @slot n_genes numeric. Shared feature count.
#' @slot params list. Reserved for downstream pipeline metadata.
#' @slot ops list. The pre-materialization prefix of the chain. Mirrors
#'   `parquetExprStore@ops` — same record schema, same `.pe_apply_op`
#'   executor. Axis-keyed payloads carry one entry per substore keyed by
#'   `uid`, so a single union-level record covers every substore. Substores
#'   must have empty `@ops` at union construction time (see constructor); the
#'   union's own `@ops` carries any subsequent recipes.
#' @slot post_ops list. The suffix of the chain, running from the first step
#'   that cannot execute in Acero onward — applied R-side to the materialized
#'   `data.table` after the `@ops` plan has run. A lowerable record can sit
#'   here legitimately: once something has forced materialization, everything
#'   after it must follow. Same pure-data record schema as `@ops`; the
#'   executor is `.pe_apply_post_ops_df()`. Axis-keyed payloads (e.g. a
#'   `multiply` op's factor vectors) carry one entry per substore keyed by
#'   `uid`. Substores must have empty `@post_ops` at construction
#'   (see constructor); at read time the union transplants its own `@ops`
#'   and `@post_ops` onto each substore via
#'   `.exprbase_inject_parent_ops()` so per-substore chunk reads stay
#'   self-sufficient. `storeWrite()` bakes the chain into on-disk values,
#'   leaving the output store with empty chains. Once `@post_ops` is
#'   non-empty, subsequent op pushes are routed here regardless of their
#'   natural phase (monotonic phase rule). See `R/utils-pestore-ops.R`.
#' @family store types
NULL

setClass("unionParquetExprStore",
    contains = c("dataStore", "parquetExprBase"),
    slots = list(
        stores   = "list",
        cell_ids = "character",
        feat_ids = "character",
        n_cells  = "numeric",
        n_genes  = "numeric",
        params   = "list",
        ops      = "list",     # chain prefix, run in Acero
        post_ops = "list"      # chain suffix, run R-side after collect
    ),
    prototype = list(
        stores   = list(),
        cell_ids = character(0L),
        feat_ids = character(0L),
        n_cells  = 0,
        n_genes  = 0,
        params   = list(),
        ops      = list(),
        post_ops = list()
    )
)

#' @name unionParquetExprStore
#' @title Construct a unionParquetExprStore
#' @description
#' Validates that all substores share an identical `feat_ids` vector,
#' concatenates `cell_ids`, and returns the union handle. Substore files
#' on disk are untouched.
#'
#' Equivalent to `cbind2()` over `parquetExprStore` objects — see also
#' the `cbind2` methods.
#' @param stores list of `parquetExprStore` objects
#' @return A [unionParquetExprStore-class] object.
#' @family store constructors
#' @export
unionParquetExprStore <- function(stores) {
    if (!is.list(stores)) {
        stop("[unionParquetExprStore] `stores` must be a list of ",
             "parquetExprStore objects", call. = FALSE)
    }
    if (length(stores) == 0L) {
        stop("[unionParquetExprStore] at least one substore is required",
             call. = FALSE)
    }
    if (!all(vapply(stores, inherits, logical(1L), "parquetExprStore"))) {
        stop("[unionParquetExprStore] all substores must be ",
             "parquetExprStore objects", call. = FALSE)
    }
    # Substores must be ops-clean. Otherwise per-substore op chains
    # would have to compose against a union view, which complicates
    # arrow translation and breaks the "ops are frozen population
    # snapshots" contract (per-substore norm would be tuned to per-
    # substore populations, not the union's). The canonical workflow
    # is: cbind raw substores → run processData on the union.
    if (any(vapply(stores, function(s)
        length(s@ops) > 0L || length(s@post_ops) > 0L,
        logical(1L)))) {
        stop("[unionParquetExprStore] one or more substores has queued ",
             "@ops or @post_ops. Materialize via ",
             "storeWrite(parquetExprStore(), s) first to bake the chain ",
             "into a fresh raw store, or run processData() on the union ",
             "after cbind. Per-substore ops are not composed across ",
             "unions.", call. = FALSE)
    }
    # feat_ids must be identical and in identical order across substores.
    f0 <- stores[[1L]]@feat_ids
    for (i in seq_along(stores)[-1L]) {
        if (!identical(stores[[i]]@feat_ids, f0)) {
            stop("[unionParquetExprStore] substore ", i,
                 "'s feat_ids differ from substore 1 (must be identical ",
                 "and in identical order).", call. = FALSE)
        }
    }
    cell_ids <- unlist(lapply(stores, function(s) s@cell_ids),
                       use.names = FALSE)
    if (anyDuplicated(cell_ids) > 0L) {
        stop("[unionParquetExprStore] duplicate cell_ids across ",
             "substores — caller must pre-prefix to ensure uniqueness.",
             call. = FALSE)
    }
    new("unionParquetExprStore",
        stores   = stores,
        cell_ids = as.character(cell_ids),
        feat_ids = f0,
        n_cells  = as.numeric(sum(vapply(stores,
                       function(s) s@n_cells, numeric(1L)))),
        n_genes  = as.numeric(length(f0))
    )
}
