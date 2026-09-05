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


# tenxH5Input ####

#' @name tenxH5Input-class
#' @title 10x cell_feature_matrix.h5 Input
#' @description
#' Wraps a 10x HDF5 sparse-matrix file (CSC layout with `data`,
#' `indices`, `indptr`, `barcodes`, `features/{id,name}`, `shape`).
#' Metadata is read eagerly at construction; the actual sparse data is
#' streamed in cell-chunks by [storeRead()] via hyperslab reads on a
#' long-lived `hdf5r::H5File` handle.
#'
#' @slot batch_cells integer. Cells per batch. Default 250,000.
#' @slot feature_id_col integer. `1L` for Ensembl ID, `2L` for gene
#'   symbol (default).
#' @family store types
NULL

setClass("tenxH5Input",
    contains = "exprInput",
    slots = list(
        batch_cells    = "integer",
        feature_id_col = "integer"
    ),
    prototype = list(
        batch_cells    = 250000L,
        feature_id_col = 2L
    )
)

#' @name tenxH5Input
#' @title Create a 10x cell_feature_matrix.h5 input
#' @description
#' Opens the h5 briefly to read `barcodes`, `features/{id,name}`, and
#' `shape`; closes immediately. The handle is reopened by `storeRead()`
#' for streaming.
#'
#' @param h5_path character. Path to `cell_feature_matrix.h5`.
#' @param feature_id_col integer. `1L` = Ensembl ID, `2L` = gene symbol
#'   (default).
#' @param batch_cells integer. Cells per batch. Default 250,000.
#' @return A `tenxH5Input` object.
#' @family store constructors
#' @export
tenxH5Input <- function(
    h5_path,
    feature_id_col = 2L,
    batch_cells    = 250000L
) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("[tenxH5Input] hdf5r is required to read 10x .h5 files. ",
             "Install with: install.packages(\"hdf5r\").", call. = FALSE)
    }
    stopifnot(file.exists(h5_path))
    h5 <- hdf5r::H5File$new(h5_path, mode = "r")
    on.exit(h5$close_all())

    root      <- names(h5)[1L]
    cell_ids  <- as.character(h5[[paste0(root, "/barcodes")]][])
    feat_id   <- as.character(h5[[paste0(root, "/features/id")]][])
    feat_name <- as.character(h5[[paste0(root, "/features/name")]][])
    feat_ids  <- .disambiguate_feat_ids(
        if (feature_id_col == 1L) feat_id else feat_name
    )
    shape     <- as.integer(h5[[paste0(root, "/shape")]][])

    new("tenxH5Input",
        path           = normalizePath(h5_path),
        params         = list(root = root),
        cell_ids       = cell_ids,
        feat_ids       = feat_ids,
        n_cells        = shape[2L],
        n_genes        = shape[1L],
        batch_cells    = as.integer(batch_cells),
        feature_id_col = as.integer(feature_id_col)
    )
}


# cellbinGefInput ####

#' @name cellbinGefInput-class
#' @title Stereo-seq cellbin GEF Input
#' @description
#' Wraps a Stereo-seq cellbin `.gef` file (HDF5 compound datasets under
#' `cellBin/`). The `cell` and `gene` tables are read in full at
#' construction (small); the compound `geneExp` is streamed gene-chunk-
#' wise via [storeRead()] using rhdf5 hyperslab reads, respecting the
#' safe-boundary chunk plan that keeps duplicate-named genes together.
#'
#' @slot batch_genes integer. Approximate raw-gene rows per batch
#'   (default 500). Actual boundaries may be expanded to keep duplicate-
#'   named gene groups intact.
#' @slot gene_column character. `"geneName"` (default) or `"geneID"`.
#' @family store types
NULL

setClass("cellbinGefInput",
    contains = "exprInput",
    slots = list(
        batch_genes = "integer",
        gene_column = "character"
    ),
    prototype = list(
        batch_genes = 500L,
        gene_column = "geneName"
    )
)

#' @name cellbinGefInput
#' @title Create a Stereo-seq cellbin GEF input
#' @param gef_path character. Path to a Stereo-seq cellbin `.gef` file.
#' @param gene_column character. `"geneName"` (default) or `"geneID"`.
#' @param batch_genes integer. Approximate raw-gene rows per batch.
#'   Default 500.
#' @return A `cellbinGefInput` object.
#' @family store constructors
#' @export
cellbinGefInput <- function(
    gef_path,
    gene_column = c("geneName", "geneID"),
    batch_genes = 500L
) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("[cellbinGefInput] rhdf5 is required for Stereo-seq .gef. ",
             "Install with: BiocManager::install(\"rhdf5\").",
             call. = FALSE)
    }
    stopifnot(file.exists(gef_path))
    gene_column <- match.arg(gene_column, c("geneName", "geneID"))

    cellDT <- data.table::setDT(rhdf5::h5read(gef_path, "cellBin/cell"))
    geneDT <- data.table::setDT(rhdf5::h5read(gef_path, "cellBin/gene"))

    n_cells  <- nrow(cellDT)
    cell_ids <- paste0("cell_", cellDT$id)

    cnt             <- as.integer(geneDT$cellCount)
    all_names       <- as.character(geneDT[[gene_column]])
    expressed_names <- all_names[cnt > 0]
    feat_ids        <- sort(unique(expressed_names))

    new("cellbinGefInput",
        path        = normalizePath(gef_path),
        # Stash the parsed tables + name_to_row + cumulative offsets so
        # storeRead doesn't re-read them. Cheap (tens of KB each).
        params      = list(
            cell_id_map = data.table::setattr(
                seq_len(n_cells), "names", as.character(cellDT$id)
            ),
            gene_cnt    = cnt,
            name_to_row = match(all_names, feat_ids),
            cum_offsets = c(0L, as.integer(cumsum(as.numeric(cnt))))
        ),
        cell_ids    = cell_ids,
        feat_ids    = feat_ids,
        n_cells     = as.integer(n_cells),
        n_genes     = length(feat_ids),
        batch_genes = as.integer(batch_genes),
        gene_column = gene_column
    )
}


# binGefInput ####

#' @name binGefInput-class
#' @title Stereo-seq bin GEF Input
#' @description
#' Wraps a Stereo-seq bin `.gef` file (`geneExp/<bin_size>/expression`
#' compound dataset). Unlike cellbin, the cell identity universe (one per
#' unique `(x, y)` coord) is not known up-front — `(x, y) -> bin_ID` is
#' assigned as new coords are encountered during the gene-chunk stream.
#' [storeRead()]'s iterator publishes the accumulated `cell_ids` /
#' `n_cells` via its accessors after iteration completes.
#'
#' @slot bin_size character. Bin size key under `geneExp/` (e.g. `"50"`).
#' @slot batch_genes integer. Approximate raw-gene rows per batch.
#' @slot gene_column character.
#' @family store types
NULL

setClass("binGefInput",
    contains = "exprInput",
    slots = list(
        bin_size    = "character",
        batch_genes = "integer",
        gene_column = "character"
    ),
    prototype = list(
        bin_size    = "50",
        batch_genes = 500L,
        gene_column = "geneName"
    )
)

#' @name binGefInput
#' @title Create a Stereo-seq bin GEF input
#' @param gef_path character.
#' @param bin_size character or integer. Bin size key under `geneExp/`
#'   (e.g. `50`, `"50"`, `"100"`).
#' @param gene_column character. `"geneName"` (default) or `"geneID"`.
#' @param batch_genes integer. Default 500.
#' @return A `binGefInput` object. `cell_ids` / `n_cells` are empty until
#'   [storeRead()] has been driven to completion.
#' @family store constructors
#' @export
binGefInput <- function(
    gef_path,
    bin_size,
    gene_column = c("geneName", "geneID"),
    batch_genes = 500L
) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
        stop("[binGefInput] rhdf5 is required for Stereo-seq .gef. ",
             "Install with: BiocManager::install(\"rhdf5\").",
             call. = FALSE)
    }
    stopifnot(file.exists(gef_path))
    gene_column <- match.arg(gene_column, c("geneName", "geneID"))
    bin_size <- as.character(bin_size)

    geneDT <- data.table::setDT(rhdf5::h5read(
        gef_path, paste0("geneExp/", bin_size, "/gene")
    ))
    cnt             <- as.integer(geneDT$count)
    all_names       <- as.character(geneDT[[gene_column]])
    expressed_names <- all_names[cnt > 0]
    feat_ids        <- sort(unique(expressed_names))

    new("binGefInput",
        path        = normalizePath(gef_path),
        params      = list(
            gene_cnt    = cnt,
            name_to_row = match(all_names, feat_ids),
            cum_offsets = c(0L, as.integer(cumsum(as.numeric(cnt)))),
            # Reference cell for the (x, y) -> bin_ID map the iterator
            # accumulates. An environment because the object is copied on
            # the way into storeWrite(), so a plain slot could not carry a
            # value back out. Bin coordinates only exist inside the
            # expression records, so this is the one chance to capture them
            # without a second full read of the gef.
            coord_env   = new.env(parent = emptyenv())
        ),
        cell_ids    = character(0L),
        feat_ids    = feat_ids,
        n_cells     = 0L,                 # accumulated during iteration
        n_genes     = length(feat_ids),
        bin_size    = bin_size,
        batch_genes = as.integer(batch_genes),
        gene_column = gene_column
    )
}


# csvWideInput ####

#' @name csvWideInput-class
#' @title Wide-format CSV Input
#' @description
#' Wraps a wide-format CSV (one row per cell, one column per feature,
#' plus a cell-ID column). The header is read at construction to
#' establish `feat_ids` and `n_genes`; cells are streamed row-chunk-wise
#' via [storeRead()] and `cell_ids` accumulate on the iterator as the
#' file is consumed.
#'
#' @slot cell_id_col character. CSV column carrying the cell barcode.
#' @slot skip_cols character. Additional non-feature columns to ignore.
#' @slot row_filter_fun function. Optional per-chunk row filter (e.g.
#'   CosMx `cell_ID != 0`).
#' @slot batch_rows integer. Cells per batch. Default 5,000.
#' @family store types
NULL

setClassUnion("FunctionOrNull", c("function", "NULL"))

setClass("csvWideInput",
    contains = "exprInput",
    slots = list(
        cell_id_col    = "character",
        skip_cols      = "character",
        row_filter_fun = "FunctionOrNull",
        batch_rows     = "integer"
    ),
    prototype = list(
        cell_id_col    = "cell_ID",
        skip_cols      = character(0L),
        row_filter_fun = NULL,
        batch_rows     = 5000L
    )
)

#' @name csvWideInput
#' @title Create a wide-format CSV input
#' @param csv_path character. Path to `.csv` or `.csv.gz`.
#' @param cell_id_col character. CSV column carrying the cell barcode
#'   (default `"cell_ID"`).
#' @param skip_cols character. Additional non-feature columns to ignore.
#' @param row_filter_fun function. Optional per-chunk row filter.
#' @param batch_rows integer. Cells per batch. Default 5,000.
#' @return A `csvWideInput` object. `cell_ids` / `n_cells` are empty
#'   until [storeRead()] has been driven to completion.
#' @family store constructors
#' @export
csvWideInput <- function(
    csv_path,
    cell_id_col    = "cell_ID",
    skip_cols      = character(0L),
    row_filter_fun = NULL,
    batch_rows     = 5000L
) {
    stopifnot(file.exists(csv_path))

    # Read header to determine feat_ids without scanning the body.
    con <- if (.is_real_gz(csv_path)) gzfile(csv_path, "r") else file(csv_path, "r")
    header_line <- readLines(con, n = 1L, warn = FALSE)
    close(con)
    if (length(header_line) == 0L) {
        stop("[csvWideInput] empty CSV at ", csv_path, call. = FALSE)
    }
    all_cols <- strsplit(header_line, ",", fixed = TRUE)[[1L]]
    all_cols <- gsub("^\"|\"$", "", all_cols)
    if (!cell_id_col %in% all_cols) {
        stop("[csvWideInput] cell_id_col '", cell_id_col,
             "' not in CSV header. Header: ",
             toString(head(all_cols, 5)), "...", call. = FALSE)
    }
    drop_cols <- unique(c(cell_id_col, skip_cols))
    feat_cols <- setdiff(all_cols, drop_cols)
    feat_ids  <- .disambiguate_feat_ids(feat_cols)

    new("csvWideInput",
        path           = normalizePath(csv_path),
        params         = list(all_cols = all_cols, feat_cols = feat_cols),
        cell_ids       = character(0L),
        feat_ids       = feat_ids,
        n_cells        = 0L,
        n_genes        = length(feat_ids),
        cell_id_col    = cell_id_col,
        skip_cols      = as.character(skip_cols),
        row_filter_fun = row_filter_fun,
        batch_rows     = as.integer(batch_rows)
    )
}

# tenxZarrInput ####

#' @name tenxZarrInput-class
#' @title 10x Xenium/Atera Zarr Expression Input
#' @description
#' Wraps the `cell_feature_matrix.zarr.zip` archive (or unzipped `.zarr`
#' tree) shipped by 10x Xenium and Atera. Construction is metadata-only:
#' the feature catalog (`cell_features/.zattrs`), the encoded cell
#' barcodes, and the CSC-by-FEATURE `indptr` are read eagerly; the
#' `indices` / `data` arrays are streamed by [storeRead()].
#'
#' Because the on-disk matrix is feature-major and the iterator contract
#' is cell-ordered, [storeRead()] reorders on the fly: a two-pass
#' cell-major placement when the triplet buffers fit the RAM budget
#' (`mode = "full"`), or bounded per-cell-block rescans when they do not
#' (`mode = "cellblock"`; forced when nnz exceeds int32).
#'
#' @slot feat_types character. Feature classes as 10x display strings
#'   ("Gene Expression", "Negative Control Probe", ...), aligned 1:1
#'   with `feat_ids`.
#' @slot mode character. "auto" (default), "full", or "cellblock".
#' @slot cells_per_block integer. Cells per window in cellblock mode
#'   (0 = derive from the RAM budget at read time).
#' @slot nnz numeric. Stored values across kept features (double — may
#'   exceed int32 at Atera scale).
#' @family store types
NULL

setClass("tenxZarrInput",
    contains = "exprInput",
    slots = list(
        feat_types      = "character",
        mode            = "character",
        cells_per_block = "integer",
        nnz             = "numeric"
    ),
    prototype = list(
        feat_types      = character(0L),
        mode            = "auto",
        cells_per_block = 0L,
        nnz             = 0
    )
)

# zarr feature_types are snake_case; translate to the 10x display strings
# the Xenium/Atera readers key their feat_type split on.
.zarr_feat_class_display <- function(x) {
    known <- c(
        gene = "Gene Expression",
        protein = "Protein Expression",
        negative_control_probe = "Negative Control Probe",
        negative_control_codeword = "Negative Control Codeword",
        unassigned_codeword = "Unassigned Codeword",
        deprecated_codeword = "Deprecated Codeword",
        genomic_control = "Genomic Control",
        aggregate_gene = "Aggregate Gene"
    )
    out <- unname(known[x])
    miss <- is.na(out)
    if (any(miss)) {
        # unseen class: Title Case the snake_case name
        out[miss] <- vapply(strsplit(x[miss], "_", fixed = TRUE),
            function(w) {
                paste(toupper(substring(w, 1L, 1L)),
                    substring(w, 2L), sep = "", collapse = " ")
            }, character(1L))
    }
    out
}

#' @name tenxZarrInput
#' @title Create a 10x Xenium/Atera zarr expression input
#' @description
#' Eagerly reads the feature catalog, cell barcodes and `indptr` from a
#' `cell_feature_matrix.zarr.zip` archive so `cell_ids` / `feat_ids` /
#' `n_cells` / `n_genes` are known up-front; `indices` / `data` are
#' streamed by [storeRead()].
#' @param zarr_path character. Path to `cell_feature_matrix.zarr.zip` or
#'   an unzipped `cell_feature_matrix.zarr` directory.
#' @param feature_id_col integer. `1L` = Ensembl ID (zarr `feature_ids`),
#'   `2L` = gene symbol (zarr `feature_keys`; default). Mirrors
#'   [mtxInput()]'s features.tsv column semantics.
#' @param mode character. "auto" (default; pick from the RAM budget),
#'   "full" (two-pass cell-major placement over full triplet buffers), or
#'   "cellblock" (bounded per-cell-block rescans).
#' @param cells_per_block integer. Cells per window in cellblock mode.
#'   `NULL` (default) derives it from the RAM budget at read time.
#' @param drop_aggregate logical. Drop `aggregate_gene` summary features
#'   (the "Total transcripts" row) from the matrix (default `TRUE`).
#' @return A `tenxZarrInput` object.
#' @family store constructors
#' @export
tenxZarrInput <- function(
    zarr_path,
    feature_id_col  = 2L,
    mode            = c("auto", "full", "cellblock"),
    cells_per_block = NULL,
    drop_aggregate  = TRUE
) {
    mode <- match.arg(mode)
    src <- .zarr_open(zarr_path)
    on.exit(.zarr_close(src), add = TRUE)

    cat_info <- .load_cfm_feature_catalog(src)
    cid <- .zarr_array(src, "cell_features/cell_id")
    cell_ids <- .encode_xenium_id(cid[, 1L], cid[, 2L])
    # indptr may be uint32 today, uint64 at Atera scale — always doubles
    indptr <- as.double(.zarr_array(src, "cell_features/indptr"))
    n_feat_orig <- length(indptr) - 1L

    feat_ids_all <- switch(as.integer(feature_id_col),
        cat_info$feature_ids,
        cat_info$feature_keys,
        stop("[tenxZarrInput] feature_id_col must be 1 (Ensembl) or 2 ",
            "(symbols)", call. = FALSE)
    )
    if (!length(feat_ids_all)) feat_ids_all <- cat_info$feature_ids
    if (length(feat_ids_all) != n_feat_orig ||
        length(cat_info$feature_types) != n_feat_orig) {
        stop("[tenxZarrInput] feature catalog length (",
            length(feat_ids_all), ") disagrees with indptr (",
            n_feat_orig, " features)", call. = FALSE)
    }

    keep <- rep(TRUE, n_feat_orig)
    if (isTRUE(drop_aggregate)) {
        keep <- cat_info$feature_types != "aggregate_gene"
    }
    col_remap <- integer(n_feat_orig) # 0 marks dropped features
    col_remap[keep] <- seq_len(sum(keep))
    # nnz over KEPT features only (what the iterator will emit)
    nnz <- sum((indptr[-1L] - indptr[-length(indptr)])[keep])

    new("tenxZarrInput",
        path            = normalizePath(zarr_path),
        params          = list(
            indptr = indptr,
            keep = keep,
            col_remap = col_remap
        ),
        cell_ids        = cell_ids,
        feat_ids        = .disambiguate_feat_ids(feat_ids_all[keep]),
        feat_types      = .zarr_feat_class_display(
            cat_info$feature_types[keep]),
        n_cells         = length(cell_ids),
        n_genes         = sum(keep),
        nnz             = nnz,
        mode            = mode,
        cells_per_block = as.integer(cells_per_block %||% 0L)
    )
}
