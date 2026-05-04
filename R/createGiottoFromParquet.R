#' @include class-parquetExprStore.R utils-parquetExprStore.R
NULL

# createGiottoFromParquet ####

#' @name createGiottoFromParquet
#' @title Create a Giotto object backed by a parquet expression matrix
#' @description
#' Convenience wrapper that constructs a `giotto` object whose expression
#' backend is a [parquetExprStore-class]. Accepts either an existing sorted
#' long-format Parquet file (or directory of chunks), or a 10x / Xenium
#' MatrixMarket triple — in the latter case the streaming converter
#' [mtx_to_parquetExprStore()] is called and **no `dgCMatrix` is ever
#' materialized**.
#'
#' Internally:
#' 1. Builds the [parquetExprStore-class] (via direct construction or via
#'    [mtx_to_parquetExprStore()]).
#' 2. Wraps it in an `exprObj` via `new("exprObj", ...)` to bypass
#'    GiottoClass's `.evaluate_expr_matrix` validator (which doesn't
#'    recognize `parquetExprStore`).
#' 3. Creates a `giotto` skeleton with a 1×1 stub matrix shaped like the
#'    parquet store, then `setExpression(g, eo)` swaps the parquet store
#'    in.
#' 4. Optionally calls [sc_recommend_chunk()] to print a chunk-size
#'    recommendation based on currently free RAM.
#'
#' @section Two ways to call:
#' **From an existing parquet** — pass `parquet_path`, `cell_ids`,
#' `feat_ids`:
#' ```r
#' g <- createGiottoFromParquet(
#'     parquet_path = "counts.parquet",
#'     cell_ids     = barcode_vec,
#'     feat_ids     = gene_vec,
#'     spatial_locs = locs_df
#' )
#' ```
#'
#' **From a 10x MatrixMarket triple** — pass `mtx_dir`:
#' ```r
#' g <- createGiottoFromParquet(
#'     mtx_dir       = "Xenium_output/cell_feature_matrix/",
#'     output_path   = "Xenium_output/expression.parquet",
#'     spatial_locs  = locs_df
#' )
#' ```
#' Streams `matrix.mtx.gz` to a sorted long-format Parquet, reads
#' barcodes/features from sidecars, then proceeds as above.
#'
#' @param parquet_path character. Path to an existing sorted long-format
#'   parquet file or directory. Mutually exclusive with `mtx_dir`.
#' @param mtx_dir character. Directory containing `matrix.mtx.gz`,
#'   `barcodes.tsv.gz`, `features.tsv.gz` (10x / Xenium layout).
#'   Mutually exclusive with `parquet_path`.
#' @param output_path character. Required when `mtx_dir` is supplied:
#'   destination path for the parquet file or directory produced by
#'   [mtx_to_parquetExprStore()].
#' @param cell_ids character. Cell barcodes, length = n_cells. Required
#'   when `parquet_path` is supplied; ignored (read from `barcodes.tsv.gz`)
#'   when `mtx_dir` is supplied.
#' @param feat_ids character. Gene / feature IDs, length = n_genes.
#'   Required when `parquet_path` is supplied; ignored (read from
#'   `features.tsv.gz`) when `mtx_dir` is supplied.
#' @param spatial_locs spatial locations passed through to
#'   [Giotto::createGiottoObject()].
#' @param chunk_size numeric. Streaming chunk size (cells per Arrow
#'   read). `NULL` (default) means call [sc_recommend_chunk()] using
#'   currently free RAM and print a guidance line. Pass an integer to
#'   skip the recommendation.
#' @param feature_id_col integer. Which column of `features.tsv` to use
#'   as the gene identifier when `mtx_dir` is supplied. Default `2L`
#'   (gene symbol).
#' @param verbose logical. Print progress messages. Default `TRUE`.
#' @param ... additional arguments passed through to
#'   [Giotto::createGiottoObject()] (e.g. `spat_unit`, `feat_type`,
#'   `instructions`).
#' @return A `giotto` object whose `@exprMat` slot in the "raw"
#'   `exprObj` is a [parquetExprStore-class].
#' @seealso [parquetExprStore()], [mtx_to_parquetExprStore()],
#'   [sc_recommend_chunk()].
#' @export
createGiottoFromParquet <- function(
    parquet_path   = NULL,
    mtx_dir        = NULL,
    output_path    = NULL,
    cell_ids       = NULL,
    feat_ids       = NULL,
    spatial_locs   = NULL,
    chunk_size     = NULL,
    feature_id_col = 2L,
    verbose        = TRUE,
    ...
) {
    if (is.null(parquet_path) && is.null(mtx_dir)) {
        stop("[createGiottoFromParquet] one of `parquet_path` or `mtx_dir` ",
             "must be supplied.", call. = FALSE)
    }
    if (!is.null(parquet_path) && !is.null(mtx_dir)) {
        stop("[createGiottoFromParquet] supply either `parquet_path` or ",
             "`mtx_dir`, not both.", call. = FALSE)
    }

    # ---- Branch 1: stream from a 10x / Xenium MatrixMarket triple --------
    if (!is.null(mtx_dir)) {
        if (is.null(output_path)) {
            stop("[createGiottoFromParquet] `output_path` is required when ",
                 "`mtx_dir` is supplied.", call. = FALSE)
        }
        # Accept either .gz or plain text variants for each sidecar
        .pick <- function(dir, base) {
            cands <- file.path(dir, c(paste0(base, ".gz"), base))
            hit   <- cands[file.exists(cands)]
            if (length(hit) == 0L) {
                stop("[createGiottoFromParquet] could not find ",
                     toString(cands), call. = FALSE)
            }
            hit[[1L]]
        }
        if (isTRUE(verbose)) {
            message("[createGiottoFromParquet] streaming MatrixMarket -> parquet ...")
        }
        pe <- mtx_to_parquetExprStore(
            mtx_path       = .pick(mtx_dir, "matrix.mtx"),
            barcodes_path  = .pick(mtx_dir, "barcodes.tsv"),
            features_path  = .pick(mtx_dir, "features.tsv"),
            output_path    = output_path,
            feature_id_col = feature_id_col,
            overwrite      = TRUE
        )
    } else {
        # ---- Branch 2: existing parquet ----------------------------------
        if (is.null(cell_ids) || is.null(feat_ids)) {
            stop("[createGiottoFromParquet] `cell_ids` and `feat_ids` are ",
                 "required when `parquet_path` is supplied.", call. = FALSE)
        }
        pe <- parquetExprStore(
            path     = normalizePath(parquet_path, mustWork = TRUE),
            cell_ids = as.character(cell_ids),
            feat_ids = as.character(feat_ids)
        )
    }

    # ---- Chunk-size recommendation --------------------------------------
    if (is.null(chunk_size)) {
        rec <- tryCatch(
            sc_recommend_chunk(
                n_cells = pe@n_cells,
                n_genes = pe@n_genes,
                density = .estimate_density(pe),
                verbose = verbose
            ),
            error = function(e) NULL
        )
        if (!is.null(rec)) pe@chunk_size <- as.numeric(rec)
    } else {
        pe@chunk_size <- as.numeric(chunk_size)
    }

    # ---- Build giotto with a stub matrix, then swap in pe ---------------
    if (isTRUE(verbose)) {
        message("[createGiottoFromParquet] assembling giotto skeleton ...")
    }
    stub <- Matrix::sparseMatrix(
        i = 1L, j = 1L, x = 1,
        dims = c(as.integer(pe@n_genes), as.integer(pe@n_cells))
    )
    rownames(stub) <- pe@feat_ids
    colnames(stub) <- pe@cell_ids

    g <- Giotto::createGiottoObject(
        expression   = stub,
        spatial_locs = spatial_locs,
        verbose      = verbose,
        ...
    )

    eo <- methods::new("exprObj",
        name      = "raw",
        exprMat   = pe,
        spat_unit = "cell",
        feat_type = "rna"
    )
    g <- GiottoClass::setExpression(g, x = eo, name = "raw", verbose = FALSE)

    if (isTRUE(verbose)) {
        message(sprintf(
            "[createGiottoFromParquet] done — %s cells x %s genes, backend = parquetExprStore",
            format(pe@n_cells, big.mark = ","),
            format(pe@n_genes, big.mark = ",")
        ))
    }
    g
}


# Helper: rough density estimate for sc_recommend_chunk. Counts entries
# in the parquet via Arrow (cheap — single scan, no value reads).
.estimate_density <- function(pe) {
    ds <- storeRead(pe)
    nnz <- tryCatch({
        ds |> dplyr::count() |> dplyr::collect() |> as.numeric()
    }, error = function(e) NA_real_)
    if (is.na(nnz) || nnz == 0) return(0.1)   # fallback
    as.numeric(nnz) / (as.numeric(pe@n_cells) * as.numeric(pe@n_genes))
}
