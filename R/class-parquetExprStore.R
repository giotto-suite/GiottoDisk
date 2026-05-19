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
#' This store is the streaming-friendly expression backend designed for the
#' scstream pipeline. It is slotted into `exprObj@exprMat` like any other
#' disk-backed store; downstream Giotto methods that recognize the class
#' will dispatch to streaming implementations.
#'
#' @slot path character. Local file path (single Parquet) or directory
#'   (one Parquet per chunk for very large datasets). Arrow's `open_dataset`
#'   handles both transparently.
#' @slot uid character. Auto-generated unique ID for artifact tracking.
#' @slot read_fun function. Preset to `arrow::open_dataset()`.
#' @slot params list. Reserved for downstream pipeline metadata
#'   (e.g. JIT scale factors after `sc_normalize`, HVG indices after
#'   `sc_hvg`).
#' @slot n_cells numeric. Number of cells in the dataset
#'   (length of `cell_ids`).
#' @slot n_genes numeric. Number of genes / features
#'   (length of `feat_ids`).
#' @slot cell_ids character. Cell barcodes; index `i` corresponds to
#'   `row_id == i` in the Parquet file.
#' @slot feat_ids character. Gene / feature IDs; index `j` corresponds to
#'   `col_id == j` in the Parquet file.
#' @slot chunk_size numeric. Default Arrow read-chunk size in cells
#'   (default 250,000).
#' @family store types
#' @seealso [parquetExprStore()], [mtx_to_parquetExprStore()]
NULL

# definitions ####

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
#'   `length(cell_idx) == length(cell_ids) == n_cells`.
#' @slot gene_idx integer. Active gene positions in the original Parquet
#'   (length 0 = no subset).
setClass("parquetExprStore",
    contains = "queryableStore",
    slots = list(
        n_cells    = "numeric",
        n_genes    = "numeric",
        cell_ids   = "character",
        feat_ids   = "character",
        cell_idx   = "integer",
        gene_idx   = "integer",
        chunk_size = "numeric"
    ),
    prototype = list(
        n_cells    = 0,
        n_genes    = 0,
        cell_ids   = character(0L),
        feat_ids   = character(0L),
        cell_idx   = integer(0L),
        gene_idx   = integer(0L),
        chunk_size = 250000
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
#' To create the Parquet from a 10x / Xenium MatrixMarket triple, use
#' [mtx_to_parquetExprStore()].
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
#' @param chunk_size numeric. Default Arrow read-chunk size in cells
#'   (default 250,000).
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
    chunk_size = 250000,
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
    new("parquetExprStore",
        path       = path,
        n_cells    = as.numeric(n_cells),
        n_genes    = as.numeric(n_genes),
        cell_ids   = as.character(cell_ids),
        feat_ids   = as.character(feat_ids),
        chunk_size = as.numeric(chunk_size),
        ...
    )
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
#' @family store types
NULL

setClass("unionParquetExprStore",
    contains = "dataStore",
    slots = list(
        stores   = "list",
        cell_ids = "character",
        feat_ids = "character",
        n_cells  = "numeric",
        n_genes  = "numeric",
        params   = "list"
    ),
    prototype = list(
        stores   = list(),
        cell_ids = character(0L),
        feat_ids = character(0L),
        n_cells  = 0,
        n_genes  = 0,
        params   = list()
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
