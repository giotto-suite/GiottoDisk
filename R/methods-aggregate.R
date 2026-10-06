#' @name calculateOverlap
#' @title Calculate Overlap
#' @description
#' Calculate features `y` that are overlapped by polygons `x`. GiottoDisk
#' provides methods for [GiottoClass::calculateOverlap()] for operating on
#' disk backed stores.
#'
#' Both dispatch paths return an `overlapPointDisk` wrapping a
#' `queryableStore` over parquet files with a fixed schema: `poly_ID`,
#' `feat_ID`, any `keep_cols`, auto-included `count` (if present),
#' `pt_tile_index`, `pt_row_index`. A row is keyed by `(pt_tile_index,
#' pt_row_index, poly_ID)`: the point, and the polygon it falls in. The
#' result is self-contained for [overlapToMatrix()] and keeps the point
#' keys for joins back to the point store.
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
#'   `NULL` (default) processes all tiles. Ignored unless `engine = "terra"`.
#' @param prune_tiles `logical` (default `FALSE`). When `y` is a
#'   `parquetGeomTileStore`, run a parallel pre-filter pass to skip point tiles
#'   that contain no polygon centroids. Enable when polygons are coarsely
#'   distributed relative to the point tile plan to avoid spawning empty workers.
#'   (For example when coarsely binning)
#' @param engine `character` one of `"terra"`, `"duckdb"` or `"sedona"`, or
#'   `NULL` (default) to use `getOption("giottodisk.spatial_query_engine")`,
#'   the same option [spatRelate()] reads. Unset or `"auto"`, that picks the
#'   first installed of sedona, duckdb and terra, except that supplying any
#'   tiling param below picks terra. `"terra"` iterates over adaptive polygon
#'   tiles (or point tiles for `parquetGeomTileStore` `y`). `"duckdb"` and
#'   `"sedona"` perform one full-dataset spatial join over the stores'
#'   [storeRead()] scans, so filters pending on `x` and `y` apply; tiling
#'   params (`threshold`, `tiles`, `pad_y`, `poly_buf_factor`, `tile_idx`,
#'   `prune_tiles`) are ignored. `"duckdb"` needs DuckDB's spatial extension
#'   (`INSTALL spatial`); `"sedona"` needs sedonadb >= 0.4.
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
#' Aggregates the output of [calculateOverlap()] into a feature x cell sparse
#' count matrix store. `store_type = "parquetExpr"` (the default) writes a
#' `parquetExprStore`, built in one DuckDB query by default (see `engine`).
#' `"bpcells"` and `"h5"` are built in Arrow through a Matrix Market
#' intermediate.
#'
#' The matrix covers the overlap's feature and cell ID universes
#' (`@feat_ids` / `@spat_ids`, or `all_feat_ids` / `all_cell_ids`): IDs with
#' no overlaps get empty rows or columns, and overlap rows outside them are
#' dropped.
#' @param x `overlapPointDisk` output from [calculateOverlap()], or a
#'   `queryableStore` of overlap rows
#' @param path `character` output directory for the matrix store
#' @param feat_id_col `character` feature ID column name (default `"feat_ID"`)
#' @param poly_id_col `character` polygon ID column name (default `"poly_ID"`)
#' @param count_col `character` (optional) column to sum instead of counting
#'   rows. Useful when feature detections carry a `count` field.
#' @param engine (`overlapPointDisk` only) `character` `"duckdb"` or
#'   `"arrow"`, or `NULL` (default): duckdb when it is installed and
#'   `store_type = "parquetExpr"`, arrow otherwise. duckdb builds only
#'   `"parquetExpr"`.
#' @param ... additional params to pass
#' @returns the matrix store, or an `exprObj` wrapping it when
#'   `output = "exprObj"`
NULL

# overlapPointDisk class ####

setClass("overlapPointDisk",
    contains = "overlapInfo",
    slots = list(
        poly_id_col = "character",
        feat_id_col = "character",
        spat_ids = "character",  # all polygon IDs (including zero-overlap)
        feat_ids = "character",  # all feature IDs (including zero-overlap)
        poly_uids = "character", # uid(s) of polygon source store -- provenance + future depends tracking
        feat_uids = "character"  # uid(s) of feature source store -- provenance + future depends tracking
    )
)

# calculateOverlap #####

#' @rdname calculateOverlap
#' @export
setMethod("calculateOverlap", signature("parquetGeomStore", "parquetGeomStore"),
    function(x, y,
        method = c("vector", "raster"),
        threshold = NULL,
        tiles = NULL,
        pad_y = 500,
        engine = NULL,
        path = .dump_tempfile(),
        poly_id_col = "poly_ID",
        feat_id_col = "feat_ID",
        keep_cols = NULL,
        ...
    ) {
    method <- match.arg(method, c("vector", "raster"))
    engine <- .resolve_overlap_engine(engine,
        tiling = !missing(threshold) || !missing(tiles) || !missing(pad_y))
    if (inherits(x, "unionParquetStore") || inherits(y, "unionParquetStore")) {
        stop(
            "[calculateOverlap] union stores are not supported; ",
            "calculate per substore",
            call. = FALSE
        )
    }
    spat_ids <- .collect_ids(x, poly_id_col)
    feat_ids <- .collect_ids(y, feat_id_col)
    result <- if (engine != "terra") {
        .calculate_overlap_sql(x, y,
            engine = engine,
            dir = path,
            poly_id_col = poly_id_col,
            feat_id_col = feat_id_col,
            keep_cols = keep_cols
        )
    } else {
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
    }
    .wrap_overlap(result, poly_id_col, feat_id_col,
        spat_ids = spat_ids, feat_ids = feat_ids,
        poly_uids = storeUID(x), feat_uids = storeUID(y))
})

#' @rdname calculateOverlap
#' @export
setMethod("calculateOverlap", signature("parquetGeomStore", "parquetGeomTileStore"),
    function(x, y,
        method = c("vector", "raster"),
        poly_buf_factor = 0.15,
        pad_y = NULL,
        tile_idx = NULL,
        prune_tiles = FALSE,
        engine = NULL,
        path = .dump_tempfile(),
        poly_id_col = "poly_ID",
        feat_id_col = "feat_ID",
        keep_cols = NULL,
        ...
    ) {
    method <- match.arg(method, c("vector", "raster"))
    engine <- .resolve_overlap_engine(engine,
        tiling = !missing(poly_buf_factor) || !missing(pad_y) ||
            !missing(tile_idx) || !missing(prune_tiles))
    spat_ids <- .collect_ids(x, poly_id_col)
    feat_ids <- .collect_ids(y, feat_id_col)
    result <- if (engine != "terra") {
        .calculate_overlap_sql(x, y,
            engine = engine,
            dir = path,
            poly_id_col = poly_id_col,
            feat_id_col = feat_id_col,
            keep_cols = keep_cols
        )
    } else {
        if (method == "raster") {
            stop("[calculateOverlap] raster method not yet implemented", call. = FALSE)
        }
        .calculate_overlap_terra_tiled(x, y,
            dir = path,
            poly_buf_factor = poly_buf_factor,
            pad_y = pad_y,
            tile_idx = tile_idx,
            prune_tiles = prune_tiles,
            poly_id_col = poly_id_col,
            feat_id_col = feat_id_col,
            keep_cols = keep_cols
        )
    }
    .wrap_overlap(result, poly_id_col, feat_id_col,
        spat_ids = spat_ids, feat_ids = feat_ids,
        poly_uids = storeUID(x), feat_uids = storeUID(y))
})

## internals ####

# Engine for calculateOverlap: the shared spatial-engine resolver, except that
# "auto" with a terra tiling param supplied means terra -- the SQL engines
# ignore those params, and auto must not drop them silently. An engine named
# explicitly (argument or option) is honoured as given.
.resolve_overlap_engine <- function(engine = NULL, tiling = FALSE) {
    auto <- is.null(engine) &&
        identical(getOption("giottodisk.spatial_query_engine", "auto"), "auto")
    resolved <- if (auto && isTRUE(tiling)) "terra" else
        .resolve_spatial_engine(engine, verb = "calculateOverlap")
    match.arg(resolved, c("terra", "duckdb", "sedona"))
}

.calculate_overlap_terra_tiled <- function(x, y,
    dir = file.path(tempdir(), .make_uid()),
    poly_buf_factor = 0.15,
    pad_y = NULL,
    tile_idx = NULL,
    prune_tiles = FALSE,
    poly_id_col = "poly_ID",
    feat_id_col = "feat_ID",
    keep_cols = NULL
) {
    write_dir <- dir

    tile_sel <- y@tiles
    if (length(tile_sel) == 0L) {
        warning("[calculateOverlap] no point tiles found", call. = FALSE)
        return(.overlap_store(dir))
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
    # Select after expanding: `$bounds` is plan geometry, which a
    # tileSelection does not expose (its `$` reads per-tile metadata), and
    # "outermost" means outermost in the plan, not in the selection.
    if (!is.null(tile_idx)) {
        tile_sel <- tile_sel[i = as.integer(tile_idx), drop = FALSE]
    }

    # polygon buffer: explicit pad_y takes priority; otherwise derive from
    # max_poly_radius. Error loudly if neither is available.
    poly_buffer <- if (!is.null(pad_y)) {
        pad_y
    } else {
        r <- .pgeom_max_poly_radius(x)
        if (is.null(r) || is.na(r) || r <= 0) {
            stop(
                "[calculateOverlap] `x@params$max_poly_radius` is missing or zero -- ",
                "padding cannot be derived automatically.\n",
                "Supply `pad_y` directly (in data units) to set the point fetch buffer.",
                call. = FALSE
            )
        }
        r * (1 + poly_buf_factor)
    }

    # resolve columns to fetch -- avoids materializing unused attributes
    feat_col_names <- colnames(y)
    extra_cols <- keep_cols %||% character(0L)
    if ("count" %in% feat_col_names &&
            !"count" %in% c(feat_id_col, extra_cols)) {
        extra_cols <- c(extra_cols, "count")
    }
    extra_cols <- intersect(extra_cols, feat_col_names)
    x_sub <- x[, poly_id_col]
    y_sub <- y[, c(feat_id_col, extra_cols, specialCols(y))]

    if (isTRUE(prune_tiles)) {
        # Pre-filter: skip point tiles with no polygon centroids.
        # Polygons are assigned to exactly one tile by centroid, so tiles with
        # zero centroids will always produce empty overlap results. Avoids
        # spawning workers for empty tiles -- significant when polygon binning is
        # coarse relative to the point tile plan.
        # Run as tileApply so counts are collected in parallel. FUN returns the
        # flat plan index .I for non-empty tiles and NULL otherwise. @indices is
        # then replaced directly, which is ordering-safe regardless of whether
        # the future backend preserves result order.
        tile_counts <- tilework::tileApply(
            x_sub,
            tiles = tile_sel,
            FUN = function(poly_q, .I) if (.dplyr_nrow(poly_q) > 0L) .I else NULL,
            get_params_x = list(output = "query")
        )
        nonempty_I <- sort(unique(as.integer(unlist(tile_counts))))
        if (length(nonempty_I) == 0L) {
            warning("[calculateOverlap] no polygon data in any point tile", call. = FALSE)
            return(.overlap_store(dir))
        }
        # Coerce to tileSelection if needed so @indices can be set directly
        if (!is(tile_sel, "tileSelection")) {
            tile_sel <- tile_sel[seq_len(length(tile_sel)), drop = FALSE]
        }
        tile_sel@indices <- nonempty_I
    }

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
        .write_parquet_file(
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

    .overlap_store(dir)
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
    write_dir <- dir

    poly_atab <- storeRead(x, output = "query")
    n_poly <- .dplyr_nrow(poly_atab)
    if (n_poly == 0L) {
        warning("[calculateOverlap] no polygon data found", call. = FALSE)
        return(.overlap_store(dir))
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
        return(.overlap_store(dir))
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
        .write_parquet_file(
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

    # nothing was written -- no overlaps found
    .overlap_store(dir)
}

.collect_ids <- function(store, col, uniques = TRUE) {
    q <- storeRead(store, output = "query") |>
        dplyr::select(!!as.name(col))
    if (uniques) q <- dplyr::distinct(q)
    dplyr::pull(q, col, as_vector = TRUE) |>
        as.character()
}


.wrap_overlap <- function(store, poly_id_col, feat_id_col,
        spat_ids = character(), feat_ids = character(),
        poly_uids = character(), feat_uids = character()) {
    new("overlapPointDisk",
        data = store,
        poly_id_col = poly_id_col,
        feat_id_col = feat_id_col,
        spat_ids = spat_ids,
        feat_ids = feat_ids,
        poly_uids = poly_uids,
        feat_uids = feat_uids
    )
}

# Build the per-tile overlap data.frame.
# poly_id_vals: polygon ID values for each overlap row (from terra::extract)
# pt_vals: terra::values(pt_sv) -- must include tile_index + row_index
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
    row.names(result_df) <- NULL
    result_df
}

# Shared SQL engine path (duckdb / sedona) for both dispatch signatures: one
# full-dataset ST_Intersects join over the two stores' storeRead() scans, so
# pending ops (subset filters, crop, spat_relate) on either store apply. The
# engine parallelises the join itself; tiling params are not needed.
#
# No row_index: the overlap is keyed by (pt_tile_index, pt_row_index,
# poly_ID). Numbering rows with row_number() OVER () forces the join onto one
# thread. adr/0019.
.calculate_overlap_sql <- function(x, y,
    engine = c("duckdb", "sedona"),
    dir = file.path(tempdir(), .make_uid()),
    poly_id_col = "poly_ID",
    feat_id_col = "feat_ID",
    keep_cols = NULL
) {
    engine <- match.arg(engine)
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

    # `count` rides along when the points carry it
    pt_cols <- colnames(y)
    extra_cols <- keep_cols %||% character(0L)
    if ("count" %in% pt_cols && !"count" %in% c(feat_id_col, extra_cols)) {
        extra_cols <- c(extra_cols, "count")
    }
    extra_cols <- intersect(extra_cols, pt_cols)

    q <- function(tbl, col) sprintf('%s."%s"', tbl, col)
    # aliases quoted: DataFusion lowercases unquoted identifiers
    sel <- c(
        sprintf('%s AS "poly_ID"', q("poly", poly_id_col)),
        sprintf('%s AS "feat_ID"', q("pt", feat_id_col)),
        vapply(extra_cols, q, character(1L), tbl = "pt", USE.NAMES = FALSE),
        'CAST(pt.tile_index AS INTEGER) AS "pt_tile_index"',
        'CAST(pt.row_index AS INTEGER) AS "pt_row_index"'
    )
    join_sql <- function(poly_from, pt_from) sprintf(
        "SELECT %s FROM %s AS poly JOIN %s AS pt ON ST_Intersects(poly.geom, pt.geom)",
        paste(sel, collapse = ", "), poly_from, pt_from)

    if (engine == "duckdb") {
        conn <- .duckdb_connect()
        on.exit(duckdb::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
        .duckdb_load_spatial(conn)
        scan <- function(store) sprintf("(%s)", dbplyr::sql_render(storeRead(
            store, output = "duckdb", duckdb_params = list(conn = conn))))
        # forward slashes: DuckDB on Windows rejects backslashed paths
        out_file <- gsub("\\\\", "/", file.path(dir, "overlap.parquet"))
        DBI::dbExecute(conn, sprintf("COPY (%s) TO '%s' (FORMAT PARQUET)",
            join_sql(scan(x), scan(y)), gsub("'", "''", out_file, fixed = TRUE)))
    } else {
        views <- tolower(paste0("gd_ov_", .make_uid(), c("_poly", "_pt", "_res")))
        on.exit(for (v in views) try(sedonadb::sd_drop_view(v), silent = TRUE),
            add = TRUE)
        sedonadb::sd_to_view(storeRead(x, output = "sedona"), views[[1L]],
            overwrite = TRUE)
        sedonadb::sd_to_view(storeRead(y, output = "sedona"), views[[2L]],
            overwrite = TRUE)
        sedonadb::sd_to_view(sedonadb::sd_sql(join_sql(
            sprintf('"%s"', views[[1L]]), sprintf('"%s"', views[[2L]]))),
            views[[3L]], overwrite = TRUE)
        # sedona hands strings back as Utf8View, which R's arrow cannot read
        types <- sedonadb::sd_collect(sedonadb::sd_sql(
            sprintf('DESCRIBE "%s"', views[[3L]])))
        cols <- types$column_name
        sel_out <- ifelse(types$data_type == "Utf8View",
            sprintf("arrow_cast(\"%s\", 'Utf8') AS \"%s\"", cols, cols),
            sprintf('"%s"', cols))
        sedonadb::sd_write_parquet(sedonadb::sd_sql(sprintf(
            'SELECT %s FROM "%s"', paste(sel_out, collapse = ", "), views[[3L]])),
            dir)
    }
    .overlap_store(dir)
}

# The overlap carrier: a queryableStore over the parquet files in `dir`.
# Overlaps are generated internally with a fixed schema and read by one
# consumer, so they get no parquetStore contract (row_index, op chain).
# adr/0019.
# A run that found nothing still leaves one zero-row file, so the store
# opens and reads as empty.
.overlap_store <- function(dir) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    if (length(list.files(dir, pattern = "[.]parquet$", recursive = TRUE)) == 0L) {
        .write_parquet_file(
            arrow::arrow_table(poly_ID = character(), feat_ID = character(),
                pt_tile_index = integer(), pt_row_index = integer()),
            file.path(dir, "overlap.parquet"))
    }
    as(fileStore(path = dir, read_fun = .overlap_read_fun), "queryableStore")
}

#' @keywords internal
#' @noRd
.overlap_read_fun <- function(path, ...) arrow::open_dataset(path, ...)

# overlapToMatrix ####

#' @rdname overlapToMatrix
#' @export
setMethod("overlapToMatrix", signature("overlapPointDisk"),
    function(x,
        name = "raw",
        sort = TRUE,
        count_col = NULL,
        store_type = getOption("giotto.gdsrc_sparsematrix_format", "parquetExpr"),
        path = .dump_tempfile(),
        output = c("store", "exprObj"),
        engine = NULL,
        ...
    ) {
    output <- match.arg(tolower(output), choices = c("store", "exprobj"))
    engine <- .resolve_matrix_engine(engine, store_type)
    # The overlap files always carry the fixed `poly_ID` / `feat_ID` schema;
    # @poly_id_col / @feat_id_col name the INPUT stores' columns.
    mat_store <- if (engine == "duckdb") {
        feat_ids <- x@feat_ids
        cell_ids <- x@spat_ids
        if (isTRUE(sort)) {
            feat_ids <- GiottoUtils::mixedsort(feat_ids)
            cell_ids <- GiottoUtils::mixedsort(cell_ids)
        }
        .overlap_to_pestore_duckdb(x@data, feat_ids, cell_ids,
            count_col = count_col, path = path)
    } else {
        overlapToMatrix(x@data,
            path = path,
            feat_id_col = "feat_ID",
            poly_id_col = "poly_ID",
            count_col = count_col,
            store_type = store_type,
            sort = sort,
            all_feat_ids = x@feat_ids,
            all_cell_ids = x@spat_ids
        )
    }
    switch(output,
        "store" = mat_store,
        "exprobj" = createExprObj(
            expression_data = if (inherits(mat_store, "parquetExprStore")) {
                mat_store
            } else {
                storeRead(mat_store)
            },
            name = name,
            spat_unit = spatUnit(x),
            feat_type = featType(x),
            provenance = prov(x)
        )
    )
})

#' @rdname overlapToMatrix
#' @export
setMethod("overlapToMatrix", signature("queryableStore"),
    function(x,
        path = .dump_tempfile(),
        feat_id_col = "feat_ID",
        poly_id_col = "poly_ID",
        count_col = NULL,
        store_type = getOption("giotto.gdsrc_sparsematrix_format", "parquetExpr"),
        sort = TRUE,
        all_feat_ids = NULL,
        all_cell_ids = NULL,
        verbose = NULL,
        ...
    ) {
    GiottoUtils::package_check("arrow")
    checkmate::assert_string(feat_id_col)
    checkmate::assert_string(poly_id_col)
    checkmate::assert_string(count_col, null.ok = TRUE)

    atab <- storeRead(x, output = "query")

    # resolve ID universes; distinct() uses plain parquet scans (safe)
    feat_ids <- if (length(all_feat_ids) > 0L) all_feat_ids else {
        atab |>
            dplyr::distinct(!!as.name(feat_id_col)) |>
            dplyr::pull(feat_id_col, as_vector = TRUE) |>
            as.character()
    }
    cell_ids <- if (length(all_cell_ids) > 0L) all_cell_ids else {
        atab |>
            dplyr::distinct(!!as.name(poly_id_col)) |>
            dplyr::pull(poly_id_col, as_vector = TRUE) |>
            as.character()
    }
    if (isTRUE(sort)) {
        feat_ids <- GiottoUtils::mixedsort(feat_ids)
        cell_ids <- GiottoUtils::mixedsort(cell_ids)
    }

    vmsg(.v = verbose, sprintf(
        "[overlapToMatrix] building COO: %d features x %d cells",
        length(feat_ids), length(cell_ids)
    ))

    # Build string -> integer lookup tables from the clean ID vectors.
    feat_df <- data.frame(i = seq_along(feat_ids))
    feat_df[[feat_id_col]] <- feat_ids
    cell_df <- data.frame(j = seq_along(cell_ids))
    cell_df[[poly_id_col]] <- cell_ids
    feat_lut <- arrow::as_arrow_table(feat_df)
    cell_lut <- arrow::as_arrow_table(cell_df)

    # Join string -> integer on raw parquet batches (string buffers are live
    # per-batch — no dangling pointer risk), then aggregate on integer keys.
    # Aggregating strings first and joining after is unsafe: Arrow's hash-
    # aggregate stores keys as utf8_view; at scale those view buffers can be
    # freed before all output batches are consumed, causing join misses.
    # Inner joins: the ID universes narrow as well as number, so overlap rows
    # outside them (a subset object's dropped cells or features) fall away.
    if (!is.null(count_col)) {
        coo_query <- atab |>
            dplyr::inner_join(feat_lut, by = feat_id_col) |>
            dplyr::inner_join(cell_lut, by = poly_id_col) |>
            dplyr::group_by(i, j) |>
            dplyr::summarize(
                n = sum(!!as.name(count_col), na.rm = TRUE),
                .groups = "drop"
            )
    } else {
        coo_query <- atab |>
            dplyr::inner_join(feat_lut, by = feat_id_col) |>
            dplyr::inner_join(cell_lut, by = poly_id_col) |>
            dplyr::count(i, j)
    }

    # parquetexpr destination: skip the Matrix Market intermediate entirely.
    # The COO is already triplets — a transmute + arrange + streaming write
    # produces a sorted parquetExprStore directly.
    if (tolower(store_type) == "parquetexpr") {
        vmsg(.v = verbose, "[overlapToMatrix] writing parquetExprStore...")
        return(.coo_to_parquetexpr(coo_query, feat_ids, cell_ids,
                                   path, verbose = verbose))
    }

    vmsg(.v = verbose, "[overlapToMatrix] writing Matrix Market...")
    mtx_dir <- .write_overlap_mtx(coo_query, feat_ids, cell_ids, path)
    .mtx_to_store(mtx_dir, store_type = store_type,
        feat_ids = feat_ids, cell_ids = cell_ids, verbose = verbose)
})

## internals ####

# Engine for the overlap -> matrix build. NULL takes duckdb and degrades to
# arrow when duckdb is not installed or the destination is not a
# parquetExprStore (the only one the duckdb build writes). sedona was measured
# and not added: on Atera, 40 s and 32-36 GB peak against duckdb's 7 s and
# 10-11 GB, same output.
.resolve_matrix_engine <- function(engine = NULL, store_type) {
    pe_dest <- tolower(store_type) == "parquetexpr"
    if (is.null(engine)) {
        return(if (pe_dest && requireNamespace("duckdb", quietly = TRUE)) {
            "duckdb"
        } else {
            "arrow"
        })
    }
    engine <- match.arg(engine, c("duckdb", "arrow"))
    if (engine == "duckdb" && !pe_dest) {
        stop("[overlapToMatrix] engine = \"duckdb\" builds only ",
            "store_type = \"parquetExpr\"; use engine = \"arrow\" for '",
            store_type, "'", call. = FALSE)
    }
    if (engine == "duckdb") GiottoUtils::package_check("duckdb")
    engine
}

# Overlap -> parquetExprStore in one DuckDB query: join the ID universes as
# integer LUTs, aggregate, sort cell-major, write. DuckDB runs the aggregate
# and the sort in parallel and spills past its memory limit; Arrow runs both
# mostly on one thread and holds the whole COO (adr/0019). The inner joins
# make the universes a filter: a narrowed @spat_ids / @feat_ids drops the
# rest of the overlap.
.overlap_to_pestore_duckdb <- function(store, feat_ids, cell_ids,
                                       count_col = NULL, path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    pe <- parquetExprStore(
        path     = normalizePath(path),
        cell_ids = cell_ids,
        feat_ids = feat_ids
    )
    partition_dir <- .idpath(pe@path, pe@uid)
    dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)

    conn <- .duckdb_connect()
    on.exit(duckdb::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
    # an in-memory database otherwise spills to .tmp/ in the working directory
    spill_dir <- file.path(tempdir(), "giottodisk_duckdb_spill")
    dir.create(spill_dir, showWarnings = FALSE)
    fwd <- function(p) gsub("'", "''", gsub("\\\\", "/", p))
    DBI::dbExecute(conn, sprintf("SET temp_directory = '%s'", fwd(spill_dir)))
    duckdb::duckdb_register(conn, "gd_feat_lut",
        data.frame(feat_ID = feat_ids, i = seq_along(feat_ids)))
    duckdb::duckdb_register(conn, "gd_cell_lut",
        data.frame(poly_ID = cell_ids, j = seq_along(cell_ids)))

    value_sql <- if (is.null(count_col)) "count(*)" else
        sprintf('sum(o."%s")', count_col)
    DBI::dbExecute(conn, sprintf(
        "COPY (
            SELECT c.j::INTEGER AS row_id, f.i::INTEGER AS col_id,
                   CAST(%s AS DOUBLE) AS value
            FROM read_parquet('%s') AS o
            JOIN gd_feat_lut AS f ON o.feat_ID = f.feat_ID
            JOIN gd_cell_lut AS c ON o.poly_ID = c.poly_ID
            GROUP BY ALL
            ORDER BY row_id, col_id
        ) TO '%s' (FORMAT PARQUET, COMPRESSION %s, ROW_GROUP_SIZE 1048576)",
        value_sql,
        fwd(file.path(store@path, "**", "*.parquet")),
        fwd(file.path(partition_dir, "part-0.parquet")),
        toupper(.parquet_compression())))

    .pestore_finalize_stats(pe)
}

# coo_query: lazy Arrow query with columns i (int), j (int), n (int)
# Streams record batches to a temp file to avoid peak memory, then prepends
# the MTX header (which requires nnz, only known after streaming completes).
.write_overlap_mtx <- function(coo_query, feat_ids, cell_ids, path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    mtx_path <- file.path(path, "matrix.mtx")
    tmp_path  <- paste0(mtx_path, ".tmp")
    on.exit(unlink(tmp_path), add = TRUE)

    reader <- arrow::as_record_batch_reader(coo_query)
    nnz <- 0L
    repeat {
        batch <- reader$read_next_batch()
        if (is.null(batch)) break
        data.table::fwrite(
            data.table::as.data.table(batch),
            file      = tmp_path,
            append    = TRUE,
            sep       = " ",
            col.names = FALSE
        )
        nnz <- nnz + nrow(batch)
    }

    # write header now that nnz is known, then append data
    con <- file(mtx_path, "w")
    writeLines("%%MatrixMarket matrix coordinate integer general", con)
    writeLines(sprintf("%d %d %d", length(feat_ids), length(cell_ids), nnz), con)
    close(con)
    file.append(mtx_path, tmp_path)

    # 10x-compatible sidecar files
    writeLines(cell_ids, file.path(path, "barcodes.tsv"))
    writeLines(feat_ids, file.path(path, "features.tsv"))

    invisible(path)
}

# convert a Matrix Market directory to a store of the requested type.
# feat_ids and cell_ids are passed directly to avoid re-reading sidecar files.
.mtx_to_store <- function(path, store_type, feat_ids, cell_ids, verbose = NULL) {
    store_type <- tolower(store_type)
    mtx_path <- file.path(path, "matrix.mtx")
    switch(store_type,
        "bpcells" = {
            GiottoUtils::package_check("BPCells")
            bp_path <- file.path(path, "bpcells")
            vmsg(.v = verbose, "[overlapToMatrix] importing MTX to BPCells...")
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

# Stream a lazy (i, j, n) COO query into a sorted parquetExprStore. Skips
# the Matrix Market intermediate — the COO is already triplets, so the
# conversion is column-rename + arrange. The arrange is required by
# parquetExprStore's sorted-by-row_id contract, which enables Arrow
# row-group skipping on cell-banded streaming reads downstream.
#
# Column mapping: parquetExprStore stores (row_id = cell_idx,
# col_id = gene_idx, value), so j -> row_id (cells), i -> col_id (feats).
.coo_to_parquetexpr <- function(coo_query, feat_ids, cell_ids, path,
                                verbose = NULL) {
    # NSE bindings
    i <- j <- n <- row_id <- col_id <- value <- NULL

    if (!dir.exists(path)) dir.create(path, recursive = TRUE)

    # Construct the parquetExprStore upfront so its auto-generated uid
    # is available for the source_id partition path (.idpath helper).
    pe <- parquetExprStore(
        path     = normalizePath(path),
        cell_ids = cell_ids,
        feat_ids = feat_ids
    )
    partition_dir <- .idpath(pe@path, pe@uid)
    dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)

    pe_query <- coo_query |>
        dplyr::transmute(
            row_id = as.integer(j),
            col_id = as.integer(i),
            value  = as.double(n)
        ) |>
        dplyr::arrange(row_id, col_id)

    reader <- arrow::as_record_batch_reader(pe_query)
    batch_idx <- 0L
    repeat {
        batch <- reader$read_next_batch()
        if (is.null(batch)) break
        if (batch$num_rows == 0L) next
        batch_idx <- batch_idx + 1L
        .write_parquet_file(
            batch,
            file.path(partition_dir,
                      sprintf("part-%d.parquet", batch_idx - 1L))
        )
    }

    # Empty overlap: write a single zero-row chunk so the dataset is
    # openable. parquetExprStore's read_fun (open_dataset) needs at
    # least one file present to resolve the schema.
    if (batch_idx == 0L) {
        empty <- data.table::data.table(
            row_id = integer(0L),
            col_id = integer(0L),
            value  = double(0L)
        )
        .write_parquet_file(empty,
            file.path(partition_dir, "part-0.parquet"))
    }

    # Parquet is written directly here rather than through `storeWrite()`, so
    # the marginal cache has to be filled explicitly.
    .pestore_finalize_stats(pe)
}
