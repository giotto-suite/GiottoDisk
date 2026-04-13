#' @name calculateOverlap
#' @title Calculate Overlap
#' @description
#' Calculate features `y` that are overlapped by polygons `x`. GiottoDisk
#' provides methods for [GiottoClass::calculateOverlap()] for operating on
#' disk backed stores.
#'
#' Both dispatch paths produce a flat `parquetStore` with a unified schema:
#' `poly_ID`, `feat_ID`, any `keep_cols`, auto-included `count` (if present),
#' `pt_tile_index`, `pt_row_index`, `row_index`. The result is self-contained
#' for [overlapToMatrix()] and retains integer join keys for optional joins back
#' to the point store.
#'
#' When `y` is a `parquetGeomTileStore`, iteration is driven by the point tiles.
#' Each point belongs to exactly one tile (no deduplication needed).
#'
#' When `y` is a flat `parquetGeomStore`, iteration is driven by adaptive polygon
#' tiles built via `quadtreePlan`.
#' @param x `parquetGeomStore`-inheriting object containing polygons
#' @param y `parquetGeomStore`-inheriting object containing features (points)
#' @param method (`engine = "terra"`) `character`. One of `"vector"` or `"raster"`.
#'   Method of overlap calculation. See [GiottoClass::calculateOverlap()] for details.
#' @param threshold (optional) `numeric` maximum number of polygons per tile
#'   when `y` is a flat store. `NULL` auto-computes via `.auto_threshold()`.
#' @param tiles (optional) seed `tilePlan` for quadtree planning when `y` is a
#'   flat store. `NULL` builds a default grid from the data extent and aspect
#'   ratio. A `freeTilePlan` (e.g. from `dry_run`) is used directly.
#' @param pad_y (optional) `numeric`. Fixed spatial padding (in data units)
#'   around each tile extent when fetching points.
#'   * `parquetGeomStore` dispatch: padding around each polygon tile. Default `500`.
#'   * `parquetGeomTileStore` dispatch: padding around each point tile. When
#'     `NULL` (default), derived from `x@params$max_poly_radius * (1 + poly_buf_factor)`.
#'     Supply a value directly when `max_poly_radius` was not recorded at write time.
#' @param poly_buf_factor (optional) `numeric`. Fractional buffer beyond
#'   `x@params$max_poly_radius` (auto-computed at write time) applied around
#'   each point tile extent when fetching points. Default `0.15`. Ignored when
#'   `pad_y` is supplied directly.
#' @param tile_idx (optional) `integerlike`. When `y` is a
#'   `parquetGeomTileStore`, restrict processing to these tile indices only.
#'   `NULL` (default) processes all tiles. Ignored for `engine = "duckdb"`.
#' @param engine `character` one of `"terra"` or `"duckdb"` (default = `"terra"`).
#'   `"terra"` iterates over adaptive polygon tiles (or point tiles for
#'   `parquetGeomTileStore` `y`). `"duckdb"` performs a single full-dataset
#'   spatial join via DuckDB's spatial extension — tiling params
#'   (`threshold`, `tiles`, `pad_y`, `poly_buf_factor`) are ignored.
#'   Requires the DuckDB spatial extension (`INSTALL spatial`).
#' @param path `character` filepath for the result store
#' @param poly_id_col `character` column in `x` holding polygon IDs. Written as
#'   `poly_ID` in the result (default = `"poly_ID"`).
#' @param feat_id_col `character` column in `y` holding feature IDs. Written as
#'   `feat_ID` in the result (default = `"feat_ID"`).
#' @param keep_cols `character` (optional) additional columns from `y` to carry
#'   into the result. `"count"` is auto-included if present in `y`.
#' @param ... additional params to pass
NULL

#' @name overlapToMatrix
#' @title Aggregate Overlap Results to Sparse Matrix
#' @description
#' Aggregates the output of [calculateOverlap()] into a feature × cell sparse
#' count matrix, written as a Matrix Market (.mtx) directory. The directory
#' layout is 10x-compatible (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) and
#' can be loaded directly by [BPCells::import_matrix_market()],
#' [Matrix::readMM()], or Python's `scipy.io.mmread`.
#'
#' To obtain a BPCells on-disk matrix pass the returned path to
#' [BPCells::import_matrix_market()] followed by [BPCells::write_matrix_dir()].
#' @param overlap_store `parquetStore` output from [calculateOverlap()]
#' @param path `character` output directory for the Matrix Market files
#' @param feat_id_col `character` feature ID column name (default `"feat_ID"`)
#' @param poly_id_col `character` polygon ID column name (default `"poly_ID"`)
#' @param count_col `character` (optional) column to sum instead of counting
#'   rows. Useful when feature detections carry a `count` field.
#' @param ... additional params to pass
#' @returns `character` path to the output directory (invisibly)
NULL

# calculateOverlap #####

#' @rdname calculateOverlap
#' @export
setMethod("calculateOverlap", signature("parquetGeomStore", "parquetGeomStore"),
    function(x, y,
        method = c("vector", "raster"),
        threshold = NULL,
        tiles = NULL,
        pad_y = 500,
        engine = c("terra", "duckdb"),
        path = .dump_tempfile(),
        poly_id_col = "poly_ID",
        feat_id_col = "feat_ID",
        keep_cols = NULL,
        ...
    ) {
    method <- match.arg(method, c("vector", "raster"))
    engine <- match.arg(engine, c("terra", "duckdb"))
    if (inherits(x, "unionParquetStore") || inherits(y, "unionParquetStore")) {
        stop(
            "[calculateOverlap] union stores are not supported; ",
            "calculate per substore",
            call. = FALSE
        )
    }
    if (engine == "duckdb") {
        return(.calculate_overlap_duckdb(x, y,
            dir = path,
            poly_id_col = poly_id_col,
            feat_id_col = feat_id_col,
            keep_cols = keep_cols
        ))
    }
    if (method == "raster") {
        stop("[calculateOverlap] raster method not yet implemented", call. = FALSE)
    }
    .calculate_overlap_terra(x, y,
        dir = path,
        threshold = threshold,
        tiles = tiles,
        pad_y = pad_y,
        poly_id_col = poly_id_col,
        feat_id_col = feat_id_col,
        keep_cols = keep_cols
    )
})

#' @rdname calculateOverlap
#' @export
setMethod("calculateOverlap", signature("parquetGeomStore", "parquetGeomTileStore"),
    function(x, y,
        method = c("vector", "raster"),
        poly_buf_factor = 0.15,
        pad_y = NULL,
        tile_idx = NULL,
        engine = c("terra", "duckdb"),
        path = .dump_tempfile(),
        poly_id_col = "poly_ID",
        feat_id_col = "feat_ID",
        keep_cols = NULL,
        ...
    ) {
    method <- match.arg(method, c("vector", "raster"))
    engine <- match.arg(engine, c("terra", "duckdb"))
    if (engine == "duckdb") {
        return(.calculate_overlap_duckdb(x, y,
            dir = path,
            poly_id_col = poly_id_col,
            feat_id_col = feat_id_col,
            keep_cols = keep_cols
        ))
    }
    if (method == "raster") {
        stop("[calculateOverlap] raster method not yet implemented", call. = FALSE)
    }
    .calculate_overlap_terra_tiled(x, y,
        dir = path,
        poly_buf_factor = poly_buf_factor,
        pad_y = pad_y,
        tile_idx = tile_idx,
        poly_id_col = poly_id_col,
        feat_id_col = feat_id_col,
        keep_cols = keep_cols
    )
})

## internals ####

.calculate_overlap_terra_tiled <- function(x, y,
    dir = file.path(tempdir(), .make_uid()),
    poly_buf_factor = 0.15,
    pad_y = NULL,
    tile_idx = NULL,
    poly_id_col = "poly_ID",
    feat_id_col = "feat_ID",
    keep_cols = NULL
) {
    result_store <- parquetStore(path = dir)
    write_dir <- file.path(dir, paste0("source_id=", result_store@uid))

    tile_sel <- y@tiles
    if (length(tile_sel) == 0L) {
        warning("[calculateOverlap] no point tiles found", call. = FALSE)
        return(result_store)
    }
    if (!is.null(tile_idx)) {
        tile_sel <- tile_sel[i = as.integer(tile_idx), drop = FALSE]
    }

    # Expand outermost tile bounds to cover polygon centroid extent.
    # Catches polygons whose centroids are outside the point tile plan but
    # whose geometry still overlaps points in the outermost tiles.
    poly_atab <- storeRead(x, output = "query")
    poly_data_ext <- .ext_to_num_vec(.dplyr_ext(poly_atab, sdimx = "x_index", sdimy = "y_index"))
    b <- tile_sel$bounds  # n x 4: xmin, xmax, ymin, ymax
    plan_ext <- c(min(b[, 1L]), max(b[, 2L]), min(b[, 3L]), max(b[, 4L]))
    expand <- c(
        poly_data_ext[[1L]] < plan_ext[[1L]],  # left
        poly_data_ext[[2L]] > plan_ext[[2L]],  # right
        poly_data_ext[[3L]] < plan_ext[[3L]],  # bottom
        poly_data_ext[[4L]] > plan_ext[[4L]]   # top
    )
    if (any(expand)) {
        if (expand[[1L]]) b[b[, 1L] == plan_ext[[1L]], 1L] <- poly_data_ext[[1L]]
        if (expand[[2L]]) b[b[, 2L] == plan_ext[[2L]], 2L] <- poly_data_ext[[2L]]
        if (expand[[3L]]) b[b[, 3L] == plan_ext[[3L]], 3L] <- poly_data_ext[[3L]]
        if (expand[[4L]]) b[b[, 4L] == plan_ext[[4L]], 4L] <- poly_data_ext[[4L]]
        tile_sel$bounds <- b
    }

    # polygon buffer: explicit pad_y takes priority; otherwise derive from
    # max_poly_radius. Error loudly if neither is available.
    poly_buffer <- if (!is.null(pad_y)) {
        pad_y
    } else {
        r <- .pgeom_max_poly_radius(x)
        if (is.null(r) || is.na(r) || r <= 0) {
            stop(
                "[calculateOverlap] `x@params$max_poly_radius` is missing or zero — ",
                "padding cannot be derived automatically.\n",
                "Supply `pad_y` directly (in data units) to set the point fetch buffer.",
                call. = FALSE
            )
        }
        r * (1 + poly_buf_factor)
    }

    # resolve columns to fetch — avoids materializing unused attributes
    feat_col_names <- colnames(y)
    extra_cols <- keep_cols %||% character(0L)
    if ("count" %in% feat_col_names &&
            !"count" %in% c(feat_id_col, extra_cols)) {
        extra_cols <- c(extra_cols, "count")
    }
    extra_cols <- intersect(extra_cols, feat_col_names)
    x_sub <- x[, poly_id_col]
    y_sub <- y[, c(feat_id_col, extra_cols, specialCols(y))]

    tile_overlap_fn <- function(poly_sv, feat_sv, .I) {
        if (is.null(poly_sv) || nrow(poly_sv) == 0L) return(NULL)
        if (is.null(feat_sv) || nrow(feat_sv) == 0L) return(NULL)

        extracted <- terra::extract(poly_sv, feat_sv)
        na_mask <- is.na(extracted[[2L]])
        if (all(na_mask)) return(NULL)
        extracted <- extracted[!na_mask, , drop = FALSE]

        pt_idx <- extracted[[1L]]
        # omit_internals = FALSE (via get_params_y): tile_index + row_index present
        pt_vals <- terra::values(feat_sv)

        result_df <- .build_overlap_df(
            poly_id_vals = extracted[[poly_id_col]],
            pt_vals = pt_vals,
            pt_idx = pt_idx,
            feat_id_col = feat_id_col,
            keep_cols = keep_cols
        )

        if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
        arrow::write_parquet(
            result_df,
            file.path(write_dir, sprintf("tile_%04d.parquet", .I))
        )
        NULL
    }

    tilework::tileApply(
        x_sub, y_sub,
        tiles = tile_sel,
        FUN = tile_overlap_fn,
        pad_y = poly_buffer,
        get_params_x = list(output = "terra"),
        get_params_y = list(output = "terra", contiguous = TRUE, omit_internals = FALSE)
    )

    if (!dir.exists(write_dir)) return(result_store)
    initialize(result_store)
}

.calculate_overlap_terra <- function(x, y,
    dir = file.path(tempdir(), .make_uid()),
    threshold = NULL,
    tiles = NULL,
    pad_y = 500,
    poly_id_col = "poly_ID",
    feat_id_col = "feat_ID",
    keep_cols = NULL
) {
    result_store <- parquetStore(path = dir)
    write_dir <- file.path(dir, paste0("source_id=", result_store@uid))

    poly_atab <- storeRead(x, output = "query")
    n_poly <- .dplyr_nrow(poly_atab)
    if (n_poly == 0L) {
        warning("[calculateOverlap] no polygon data found", call. = FALSE)
        return(result_store)
    }

    threshold <- threshold %||% .auto_threshold(n_poly, type = "polygons")
    if (is.null(tiles)) {
        tiles <- tilework::tilePlan("spatial")
        data_ext <- .dplyr_ext(poly_atab, sdimx = "x_index", sdimy = "y_index")
        terra::ext(tiles) <- data_ext
        erange <- range(data_ext)
        length(tiles) <- round(max(erange) / min(erange)) * 4L
    }

    fp <- tilework::quadtreePlan(x,
        tiles = tiles,
        threshold = threshold
    )

    nonempty <- which(fp@metadata$n_records > 0L)
    if (length(nonempty) == 0L) {
        warning("[calculateOverlap] no polygon data found", call. = FALSE)
        return(result_store)
    }
    tile_sel <- fp[i = as.integer(nonempty), drop = FALSE]

    # resolve point columns to fetch upfront
    feat_col_names <- colnames(y)
    extra_cols <- keep_cols %||% character(0L)
    if ("count" %in% feat_col_names &&
            !"count" %in% c(feat_id_col, extra_cols)) {
        extra_cols <- c(extra_cols, "count")
    }
    extra_cols <- intersect(extra_cols, feat_col_names)

    x_sub <- x[, poly_id_col]
    y_sub <- y[, c(feat_id_col, extra_cols, specialCols(y))]

    tile_overlap_fn <- function(poly_sv, .I) {
        if (is.null(poly_sv) || nrow(poly_sv) == 0L) return(NULL)
        # omit_internals = FALSE: need tile_index + row_index from point values
        feat_sv <- getBoundedData(y_sub, terra::ext(poly_sv) + pad_y,
            output = "terra", omit_internals = FALSE)
        if (is.null(feat_sv) || nrow(feat_sv) == 0L) return(NULL)

        extracted <- terra::extract(poly_sv, feat_sv)
        na_mask <- is.na(extracted[[2L]])
        if (all(na_mask)) return(NULL)
        extracted <- extracted[!na_mask, , drop = FALSE]

        pt_idx <- extracted[[1L]]
        pt_vals <- terra::values(feat_sv)

        result_df <- .build_overlap_df(
            poly_id_vals = extracted[[poly_id_col]],
            pt_vals = pt_vals,
            pt_idx = pt_idx,
            feat_id_col = feat_id_col,
            keep_cols = keep_cols
        )

        if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE)
        arrow::write_parquet(
            result_df,
            file.path(write_dir, sprintf("tile_%04d.parquet", .I))
        )
        NULL
    }

    tilework::tileApply(x_sub,
        tiles = tile_sel,
        FUN = tile_overlap_fn,
        get_params_x = list(output = "terra")
        # omit_internals = TRUE (default): poly_sv only needs poly_id_col
    )

    # nothing was written — no overlaps found
    if (!dir.exists(write_dir)) return(result_store)
    initialize(result_store)
}

# Build the per-tile overlap data.frame.
# poly_id_vals: polygon ID values for each overlap row (from terra::extract)
# pt_vals: terra::values(pt_sv) — must include tile_index + row_index
#   (point store fetched with omit_internals = FALSE)
# pt_idx: point row indices from extracted[[1L]]
# feat_id_col, keep_cols: column specs
.build_overlap_df <- function(poly_id_vals, pt_vals, pt_idx,
    feat_id_col, keep_cols) {
    extra_cols <- keep_cols %||% character(0L)
    # auto-include "count" if present and not already requested
    if ("count" %in% names(pt_vals) &&
            !"count" %in% c(feat_id_col, extra_cols)) {
        extra_cols <- c(extra_cols, "count")
    }
    extra_cols <- intersect(extra_cols, names(pt_vals))

    result_df <- data.frame(
        poly_ID = poly_id_vals,
        feat_ID = pt_vals[[feat_id_col]][pt_idx],
        stringsAsFactors = FALSE
    )
    if (length(extra_cols) > 0L) {
        result_df <- cbind(result_df,
            pt_vals[pt_idx, extra_cols, drop = FALSE])
    }
    result_df$pt_tile_index <- as.integer(pt_vals$tile_index[pt_idx])
    result_df$pt_row_index  <- as.integer(pt_vals$row_index[pt_idx])
    result_df$row_index     <- seq_len(nrow(result_df))
    row.names(result_df) <- NULL
    result_df
}

# Shared DuckDB engine path for both dispatch signatures.
# Reads both stores via glob + hive_partitioning, performs a full-dataset
# spatial join (ST_Intersects on WKB geometries), and writes the unified
# result schema to a single parquet. DuckDB handles internal parallelism;
# tiling params (threshold, tiles, pad_y, poly_buf_factor) are not needed.
.calculate_overlap_duckdb <- function(x, y,
    dir = file.path(tempdir(), .make_uid()),
    poly_id_col = "poly_ID",
    feat_id_col = "feat_ID",
    keep_cols = NULL
) {
    result_store <- parquetStore(path = dir)
    write_dir <- file.path(dir, paste0("source_id=", result_store@uid))
    if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE)

    conn <- .duckdb_connect()
    on.exit(duckdb::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    .duckdb_load_spatial(conn)

    # DuckDB requires forward-slash paths in glob expressions
    .fwd <- function(p) gsub("\\\\", "/", p)
    poly_glob <- .fwd(file.path(x@path, "**", "*.parquet"))
    pt_glob   <- .fwd(file.path(y@path, "**", "*.parquet"))

    # Inspect point store schema to resolve keep_cols / auto-include count
    pt_schema <- DBI::dbGetQuery(conn, sprintf(
        "DESCRIBE SELECT * FROM read_parquet('%s', hive_partitioning = true)",
        pt_glob
    ))
    pt_col_names <- pt_schema$column_name

    extra_cols <- keep_cols %||% character(0L)
    if ("count" %in% pt_col_names &&
            !"count" %in% c(feat_id_col, extra_cols)) {
        extra_cols <- c(extra_cols, "count")
    }
    extra_cols <- intersect(extra_cols, pt_col_names)

    # Build SELECT clause — quote identifiers to handle arbitrary column names
    q <- function(tbl, col) sprintf('%s."%s"', tbl, col)

    sel <- c(
        sprintf('%s AS poly_ID', q("poly", poly_id_col)),
        sprintf('%s AS feat_ID', q("pt", feat_id_col))
    )
    if (length(extra_cols) > 0L) {
        sel <- c(sel, vapply(extra_cols, q, character(1L), tbl = "pt"))
    }
    sel <- c(sel,
        'CAST(pt.tile_index AS INTEGER) AS pt_tile_index',
        'CAST(pt.row_index  AS INTEGER) AS pt_row_index',
        'CAST(row_number() OVER () AS INTEGER) AS row_index'
    )
    select_sql <- paste(sel, collapse = ",\n        ")

    out_file <- .fwd(file.path(write_dir, "overlap.parquet"))

    sql <- sprintf(
        "COPY (
    SELECT
        %s
    FROM read_parquet('%s', hive_partitioning = true) AS poly
    JOIN read_parquet('%s', hive_partitioning = true) AS pt
      ON ST_Intersects(poly.geom, pt.geom)
) TO '%s' (FORMAT PARQUET)",
        select_sql, poly_glob, pt_glob, out_file
    )

    DBI::dbExecute(conn, sql)
    initialize(result_store)
}

# overlapToMatrix ####

#' @rdname overlapToMatrix
#' @export
setMethod("overlapToMatrix", signature("parquetStore"),
    function(x,
        path = .dump_tempfile(),
        feat_id_col = "feat_ID",
        poly_id_col = "poly_ID",
        count_col = NULL,
        store_type = getOption("giotto.gdsrc_matrix_format", "bpcells"),
        ...
    ) {
    GiottoUtils::package_check("arrow")
    checkmate::assert_string(feat_id_col)
    checkmate::assert_string(poly_id_col)
    checkmate::assert_string(count_col, null.ok = TRUE)

    # aggregate via Arrow — only the COO count table is collected
    atab <- storeRead(x, output = "query")
    if (!is.null(count_col)) {
        agg <- atab |>
            dplyr::group_by(
                !!as.name(feat_id_col),
                !!as.name(poly_id_col)
            ) |>
            dplyr::summarize(
                n = sum(!!as.name(count_col), na.rm = TRUE),
                .groups = "drop"
            ) |>
            dplyr::collect()
    } else {
        agg <- atab |>
            dplyr::count(
                !!as.name(feat_id_col),
                !!as.name(poly_id_col)
            ) |>
            dplyr::collect()
    }

    feat_ids <- sort(unique(agg[[feat_id_col]]))
    cell_ids <- sort(unique(agg[[poly_id_col]]))
    i <- match(agg[[feat_id_col]], feat_ids)
    j <- match(agg[[poly_id_col]], cell_ids)
    vals <- as.integer(agg$n)

    mtx_dir <- .write_overlap_mtx(i, j, vals, feat_ids, cell_ids, path)
    .mtx_to_store(mtx_dir, store_type = store_type,
        feat_ids = feat_ids, cell_ids = cell_ids)
})

## internals ####

.write_overlap_mtx <- function(i, j, vals, feat_ids, cell_ids, path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    mtx_path <- file.path(path, "matrix.mtx")

    # header — nnz known upfront from aggregated COO
    con <- file(mtx_path, "w")
    writeLines("%%MatrixMarket matrix coordinate integer general", con)
    writeLines(
        sprintf("%d %d %d", length(feat_ids), length(cell_ids), length(vals)),
        con
    )
    close(con)

    # stream triplets via fwrite (no in-memory sparse matrix object)
    data.table::fwrite(
        data.table::data.table(i = i, j = j, x = vals),
        file = mtx_path,
        append = TRUE,
        sep = " ",
        col.names = FALSE
    )

    # 10x-compatible sidecar files
    writeLines(cell_ids, file.path(path, "barcodes.tsv"))
    writeLines(feat_ids, file.path(path, "features.tsv"))

    invisible(path)
}

# convert a Matrix Market directory to a store of the requested type.
# feat_ids and cell_ids are passed directly to avoid re-reading sidecar files.
.mtx_to_store <- function(path, store_type, feat_ids, cell_ids) {
    store_type <- tolower(store_type)
    mtx_path <- file.path(path, "matrix.mtx")
    switch(store_type,
        "bpcells" = {
            GiottoUtils::package_check("BPCells")
            bp_path <- file.path(path, "bpcells")
            # import_matrix_market streams .mtx to a lazy BPCells matrix;
            # write_matrix_dir streams that to the binary CSC format on disk
            mat <- BPCells::import_matrix_market(mtx_path)
            rownames(mat) <- feat_ids
            colnames(mat) <- cell_ids
            BPCells::write_matrix_dir(mat, dir = bp_path)
            bpcMatrixStore(path = bp_path)
        },
        "h5" = {
            GiottoUtils::package_check("HDF5Array")
            mat <- Matrix::readMM(mtx_path)
            rownames(mat) <- feat_ids
            colnames(mat) <- cell_ids
            h5_path <- file.path(path, "matrix.h5")
            store <- h5ArrayStore(path = h5_path)
            storeWrite(store, mat)
        },
        stop(
            sprintf("[overlapToMatrix] unsupported store_type: '%s'", store_type),
            call. = FALSE
        )
    )
}
