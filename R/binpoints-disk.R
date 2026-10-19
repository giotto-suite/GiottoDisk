#' @include class-parquetStore.R class-parquetExprStore.R
NULL

# Store-backed giottoBinPoints ####
#
# On a backed object a giottoBinPoints holds its bins and counts as stores:
# `@spatial` a parquetGeomTileStore with one point per bin with counts, whose
# `feat_ID` is the bin's on-disk row_id in `@counts`, a parquetExprStore of
# features x bins. `@bid` and `@pmap` stay empty: the counts store's cell IDs
# are the bin IDs. GiottoClass's giottoBinPoints methods hand store-backed
# objects to the methods here.

.gbp_store_backed <- function(x) {
    !is.null(x@spatial) && !inherits(x@spatial, "SpatVector")
}

# A store-backed giottoBinPoints from a bins' expression store and their
# positions. `xy` has `sdimx`, `sdimy` and either `cell_ID` or `row_id` (the
# bin's on-disk row in `expr`); bins without counts, or not in `expr`, are
# left out. `bins` reuses the points store of another feature type of the
# same matrix instead of writing one.
.binpoints_disk <- function(expr, xy, gsource, feat_type = "rna",
                            bins = NULL, verbose = NULL) {
    checkmate::assert_class(expr, "parquetExprStore")
    checkmate::assert_string(feat_type)
    if (is.null(bins)) bins <- .binpoints_store(expr, xy, gsource, verbose)
    checkmate::assert_class(bins, "parquetGeomTileStore")
    new("giottoBinPoints",
        spatial = bins,
        counts = expr,
        fid = expr@feat_ids,
        bid = character(0L),
        pmap = integer(0L),
        feat_type = feat_type,
        compact = TRUE
    )
}

# One point per bin with counts, keyed by the bin's on-disk row in `expr`
.binpoints_store <- function(expr, xy, gsource, verbose = NULL) {
    if (inherits(xy, "spatLocsObj")) xy <- xy[]
    checkmate::assert_names(names(xy), must.include = c("sdimx", "sdimy"))

    if ("row_id" %in% names(xy)) {
        row_id <- as.integer(xy$row_id)
        keep <- which(!is.na(row_id))
        row_id <- row_id[keep]
    } else {
        checkmate::assert_names(names(xy), must.include = "cell_ID")
        pos <- if (identical(xy$cell_ID, expr@cell_ids)) {
            seq_along(expr@cell_ids)
        } else {
            match(xy$cell_ID, expr@cell_ids)
        }
        keep <- which(!is.na(pos))
        row_id <- if (length(expr@cell_idx)) {
            expr@cell_idx[pos[keep]]
        } else {
            pos[keep]
        }
        rm(pos)
    }
    # bins without counts add nothing to an overlap
    has <- .pe_rows_with_counts(expr, row_id)
    keep <- keep[has]
    row_id <- row_id[has]

    f <- tempfile("bins_", fileext = ".parquet")
    on.exit(unlink(f), add = TRUE)
    arrow::write_parquet(arrow::arrow_table(
        feat_ID = row_id,
        x = xy$sdimx[keep],
        y = xy$sdimy[keep]
    ), f)
    rm(keep, row_id)
    # a package-level read_fun: a closure here would carry `xy` to the tile
    # workers
    qs <- as(fileStore(path = f, read_fun = .overlap_read_fun),
        "queryableStore")
    sourceWrite(gsource, qs,
        store_type = "parquetGeomTile",
        type = "points",
        id_col = "feat_ID",
        sdimx = "x",
        sdimy = "y",
        verbose = verbose
    )
}

# For on-disk row ids `rows`: does the row hold any counts. From the cached
# marginals when the store has them, otherwise one scan of its row ids.
.pe_rows_with_counts <- function(pe, rows) {
    if (!length(rows)) return(logical(0L))
    nnz <- pe@stats$row_nnz
    if (length(nnz) >= max(rows)) return(nnz[rows] > 0L)
    row_id <- NULL  # NSE
    ids <- storeRead(pe, output = "query") |>
        dplyr::distinct(row_id) |>
        dplyr::collect()
    ids <- as.integer(ids$row_id)
    present <- logical(max(c(rows, ids)))
    present[ids] <- TRUE
    present[rows]
}

#' @title Table of store-backed bin points
#' @name binPointsTable
#' @description
#' The counts of a store-backed `giottoBinPoints`, one row per feature and
#' bin, as `as.data.table()` gives them for in-memory bin points (and calls
#' this for store-backed ones). A feature selection or crop pending on the
#' stores applies.
#' @param x store-backed `giottoBinPoints`
#' @param geom `"XY"` adds the bins' `x` and `y`; `""` (default) leaves them
#'   out
#' @returns a `data.table` with columns `feat_ID`, `count` and, for
#'   `geom = "XY"`, `x` and `y`
#' @export
binPointsTable <- function(x, geom = "") {
    checkmate::assert_class(x, "giottoBinPoints")
    if (!.gbp_store_backed(x)) {
        stop("[binPointsTable] `x` is not store-backed", call. = FALSE)
    }
    conn <- .duckdb_connect()
    on.exit(duckdb::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    .duckdb_load_spatial(conn) # a queued crop on the bins is spatial SQL
    scan <- function(store) dbplyr::sql_render(storeRead(store,
        output = "duckdb", duckdb_params = list(conn = conn)))
    pe <- x@counts
    duckdb::duckdb_register(conn, "gd_bp_feat", data.frame(
        col_id = .pe_orig_col(seq_along(pe@feat_ids), pe),
        feat_ID = pe@feat_ids))
    xy <- if (identical(geom, "XY")) ', b.x_index AS "x", b.y_index AS "y"'
    dt <- DBI::dbGetQuery(conn, sprintf(paste(
        'SELECT f.feat_ID AS "feat_ID", e.value AS "count"%s',
        "FROM (%s) AS b",
        "JOIN (%s) AS e ON e.row_id = CAST(b.feat_ID AS INTEGER)",
        "JOIN gd_bp_feat AS f ON e.col_id = f.col_id"),
        xy %||% "", scan(x@spatial), scan(pe)))
    data.table::setDT(dt)[]
}

#' @rdname sourceWrite-gDirSource
#' @details A `giottoBinPoints` is written as its two stores: the counts as a
#'   `parquetExprStore` of features x bins, and one point per bin with counts.
#'   The result is the store-backed `giottoBinPoints`.
#' @export
setMethod("sourceWrite", signature("gDirSource", "giottoBinPoints"),
    function(src, data, meta = NULL, ...) {
        if (.gbp_store_backed(data)) return(data)
        d <- data@counts
        m <- Matrix::sparseMatrix(i = d$i, j = d$j, x = d$x,
            dims = c(length(data@fid), length(data@bid)),
            dimnames = list(data@fid, data@bid))
        pe <- sourceWrite(src, m, meta = meta, ...)
        xy <- terra::crds(data@spatial)
        .binpoints_disk(pe,
            xy = data.frame(cell_ID = data@bid[data@pmap],
                            sdimx = xy[, 1L], sdimy = xy[, 2L]),
            gsource = src,
            feat_type = data@feat_type)
    }
)
