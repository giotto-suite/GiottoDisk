#' @include class-dataStore.R
NULL

# Marker classes for raw expression-matrix inputs. Each one wraps an
# on-disk source (10x mtx triple, 10x h5, GEF, wide CSV) and exposes a
# uniform `storeRead()` that yields a stateful batch iterator over
# `(row_id, col_id, value)` triplets. `storeWrite(parquetExprStore, .)`
# consumes the iterator and produces a sorted long-format parquet at the
# pre-allocated store path.
#
# Inheriting from `fileStore` gives us @path, @uid, @params, @read_fun
# for free, plus automatic compatibility with gDirSource adoption / hashing.

# exprInput (virtual base) ####

#' @name exprInput-class
#' @title Expression Matrix Input (virtual)
#' @description
#' Virtual base for raw expression-matrix sources. Subclasses describe a
#' format-specific on-disk layout and expose a batch iterator via
#' [storeRead()]; [storeWrite()] on a [parquetExprStore-class] consumes
#' that iterator to materialise a sorted long-format parquet.
#'
#' Metadata slots (`cell_ids`, `feat_ids`, `n_cells`, `n_genes`) are
#' populated eagerly at construction when the format allows (e.g. mtx
#' reads its sidecars). For formats where cell identity is only known
#' during the stream (e.g. wide CSV), the iterator returned by
#' `storeRead()` is responsible for mutating these as it advances.
#' @slot cell_ids character. Cell barcodes (length `n_cells`).
#' @slot feat_ids character. Gene / feature IDs (length `n_genes`).
#' @slot n_cells integer. Total cells.
#' @slot n_genes integer. Total features.
#' @family store types
NULL

setClass("exprInput",
    contains = c("fileStore", "VIRTUAL"),
    slots = list(
        cell_ids = "character",
        feat_ids = "character",
        n_cells  = "integer",
        n_genes  = "integer"
    ),
    prototype = list(
        cell_ids = character(0L),
        feat_ids = character(0L),
        n_cells  = 0L,
        n_genes  = 0L
    )
)


# mtxInput ####

#' @name mtxInput-class
#' @title 10x / Xenium MatrixMarket Triple Input
#' @description
#' Wraps a `matrix.mtx[.gz]` + `barcodes.tsv[.gz]` + `features.tsv[.gz]`
#' triple. Sidecars are read at construction to populate `cell_ids` /
#' `feat_ids`; the mtx itself is opened lazily by `storeRead()` and
#' streamed in batches of `batch_lines` triplets via a single
#' long-lived gzip / file connection.
#'
#' @slot batch_lines integer. Triplets per batch. Default 5,000,000
#'   (~120 MB peak per batch).
#' @slot feature_id_col integer. Column of `features.tsv` to use as the
#'   gene identifier. `1L` for Ensembl ID, `2L` for gene symbol
#'   (default).
#' @section Input layout:
#' `@path` is the matrix.mtx[.gz] file. The two sidecar paths live in
#' `@params$barcodes_path` and `@params$features_path` (set by the
#' constructor; can be overridden manually).
#' @family store types
NULL

setClass("mtxInput",
    contains = "exprInput",
    slots = list(
        batch_lines    = "integer",
        feature_id_col = "integer"
    ),
    prototype = list(
        batch_lines    = 5000000L,
        feature_id_col = 2L
    )
)

#' @name mtxInput
#' @title Create a 10x / Xenium MatrixMarket triple input
#' @description
#' Eagerly reads `barcodes.tsv` and `features.tsv` so `cell_ids`,
#' `feat_ids`, `n_cells`, `n_genes` are known up-front. The mtx file
#' itself is not opened until [storeRead()].
#'
#' Three input modes are supported:
#'
#' 1. **Directory** (10x / Xenium `cell_feature_matrix/` layout):
#'    `mtxInput(dir)` — auto-resolves `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`,
#'    `features.tsv[.gz]` inside the directory.
#' 2. **Matrix file with sidecars resolvable in its parent**:
#'    `mtxInput(mtx_path)` — used when the three files live as siblings.
#' 3. **Bare `.mtx` with explicit IDs**:
#'    `mtxInput(mtx_path, cell_ids = ..., feat_ids = ...)` — for non-10x
#'    MatrixMarket files where the barcodes/features sidecars aren't
#'    present. `cell_ids` and `feat_ids` must both be supplied.
#'
#' Explicit `cell_ids` / `feat_ids` always override sidecar resolution.
#'
#' @param mtx_path character. Path to a directory (10x layout) or a
#'   `matrix.mtx[.gz]` file.
#' @param barcodes_path character. Path to `barcodes.tsv[.gz]`. Default:
#'   auto-resolved relative to `mtx_path`. Ignored if `cell_ids` is
#'   supplied.
#' @param features_path character. Path to `features.tsv[.gz]`. Default:
#'   auto-resolved relative to `mtx_path`. Ignored if `feat_ids` is
#'   supplied.
#' @param cell_ids character. Explicit cell barcodes. If supplied,
#'   `barcodes_path` is not consulted.
#' @param feat_ids character. Explicit feature IDs. If supplied,
#'   `features_path` and `feature_id_col` are not consulted.
#' @param feature_id_col integer. `1L` = Ensembl ID, `2L` = gene symbol
#'   (default). Used only when reading from `features_path`.
#' @param batch_lines integer. Triplets per batch. Default 5,000,000.
#' @return An `mtxInput` object.
#' @family store constructors
#' @export
mtxInput <- function(
    mtx_path,
    barcodes_path  = NULL,
    features_path  = NULL,
    cell_ids       = NULL,
    feat_ids       = NULL,
    feature_id_col = 2L,
    batch_lines    = 5000000L
) {
    stopifnot(file.exists(mtx_path))
    # If `mtx_path` is a directory (10x / Xenium cell_feature_matrix layout),
    # auto-resolve matrix.mtx[.gz] inside it; otherwise treat as the mtx file
    # and resolve sidecars from its parent dir.
    if (dir.exists(mtx_path)) {
        mtx_dir  <- mtx_path
        mtx_path <- .resolve_sidecar(mtx_dir, "matrix.mtx")
    } else {
        mtx_dir  <- dirname(mtx_path)
    }
    # Partial-explicit IDs is almost always a mistake.
    if (xor(is.null(cell_ids), is.null(feat_ids))) {
        stop("[mtxInput] supply both `cell_ids` and `feat_ids`, or neither.",
             call. = FALSE)
    }

    if (is.null(cell_ids)) {
        barcodes_path <- barcodes_path %||% .resolve_sidecar(mtx_dir, "barcodes.tsv")
        stopifnot(file.exists(barcodes_path))
        cell_ids <- readLines(barcodes_path)
        barcodes_path <- normalizePath(barcodes_path)
    } else {
        cell_ids <- as.character(cell_ids)
        barcodes_path <- NA_character_
    }

    if (is.null(feat_ids)) {
        features_path <- features_path %||% .resolve_sidecar(mtx_dir, "features.tsv")
        stopifnot(file.exists(features_path))
        features <- data.table::fread(features_path, header = FALSE,
                                      sep = "\t", quote = "")
        if (feature_id_col > ncol(features)) {
            stop("[mtxInput] feature_id_col (", feature_id_col,
                 ") exceeds number of columns in ", features_path,
                 " (", ncol(features), ").", call. = FALSE)
        }
        feat_ids <- .disambiguate_feat_ids(
            as.character(features[[feature_id_col]])
        )
        features_path <- normalizePath(features_path)
    } else {
        feat_ids <- .disambiguate_feat_ids(as.character(feat_ids))
        features_path <- NA_character_
    }

    new("mtxInput",
        path           = normalizePath(mtx_path),
        params         = list(
            barcodes_path = barcodes_path,
            features_path = features_path
        ),
        cell_ids       = cell_ids,
        feat_ids       = feat_ids,
        n_cells        = length(cell_ids),
        n_genes        = length(feat_ids),
        batch_lines    = as.integer(batch_lines),
        feature_id_col = as.integer(feature_id_col)
    )
}

# helpers ####

# Accept either .gz or plain sibling.
.resolve_sidecar <- function(dir, base) {
    cands <- file.path(dir, c(paste0(base, ".gz"), base))
    hit   <- cands[file.exists(cands)]
    if (length(hit) == 0L) {
        stop("[exprInput] could not find ", toString(cands), call. = FALSE)
    }
    hit[[1L]]
}

# True gzip vs .gz-extension-but-actually-plain.
.is_real_gz <- function(path) {
    if (!grepl("\\.gz$", path, ignore.case = TRUE)) return(FALSE)
    magic <- readBin(path, what = "raw", n = 2L)
    length(magic) == 2L && magic[1L] == as.raw(0x1f) && magic[2L] == as.raw(0x8b)
}