#' @name plot
#' @title Visualize a Store
#' @description
#' Plot a parquet-backed geometry store. Routing depends on `geomtype(x)`:
#'
#' * `"points"`: default path rasterizes points into a bin matrix via
#'   [rasterize()] (one Arrow `GROUP BY` over spatial bins) and plots the
#'   resulting raster. Mirrors the `giottoPoints` plot API.
#' * `"polygons"` / `"lines"`: always uses the vector path — collects via
#'   `storeRead(., output = "terra")` (with `sample_max` regular-stride
#'   downsample) and plots the resulting `SpatVector`. `rasterize` is
#'   points-only at present, and polygon/line plots are usually wanted as
#'   geometry anyway.
#'
#' Signature is `parquetGeomBase` so the method covers both single-store
#' (`parquetGeomStore`, `parquetGeomTileStore`) and
#' `unionParquetGeomStore` — both branches inherit `parquetGeomBase` and
#' both have `storeRead(., output = "terra")` methods.
#' @param feats `character` (optional, points only). Restrict to one or
#'   more `feat_ID`s before binning. Composed as a `subset(x, feat_ID
#'   %in% feats)` filter op (queued, executed lazily at Arrow read time).
#' @param raster `logical`. Points only. `TRUE` (default) bins points to
#'   a raster and plots that; `FALSE` collects points as a `SpatVector`
#'   (subject to `sample_max`) and plots that. Ignored for non-points
#'   geomtypes.
#' @param raster_size `integer`. Points raster path: major-axis pixel
#'   count for the raster template. Default `600`.
#' @param count `logical`. Points raster path: `TRUE` (default) plots
#'   per-cell counts; `FALSE` plots binary presence.
#' @param sample_max `integer`-like. Vector path only. Maximum number of
#'   geometries to plot (regularly sampled). `NULL` disables sampling.
#'   Default `getOption("giottodisk.plot_sample_max", 1e5)`.
#' @param ... forwarded to `terra::plot()`.
#' @export
setMethod("plot", signature("parquetGeomBase", "missing"),
    function(x,
        feats = NULL,
        raster = TRUE,
        raster_size = 600L,
        count = TRUE,
        sample_max = getOption("giottodisk.plot_sample_max", 1e5),
        ...) {

    is_points <- identical(geomtype(x), "points")

    if (!is.null(feats)) {
        if (!is_points) {
            stop("[plot(parquetGeomBase)] `feats` filter applies to ",
                "points stores only (current geomtype: \"",
                geomtype(x), "\"). Use subset(x, ...) for ",
                "polygon/line filtering.", call. = FALSE)
        }
        filt_expr <- bquote(feat_ID %in% .(feats))
        x <- subset(x, filt_expr, quote = FALSE)
    }

    if (isTRUE(raster) && is_points) {
        # Build a template raster matching the store's extent and aspect.
        ex <- ext(x)
        xrange <- ex$xmax - ex$xmin
        yrange <- ex$ymax - ex$ymin
        if (xrange >= yrange) {
            ncols <- as.integer(raster_size)
            nrows <- max(1L,
                as.integer(round(raster_size * yrange / xrange)))
        } else {
            nrows <- as.integer(raster_size)
            ncols <- max(1L,
                as.integer(round(raster_size * xrange / yrange)))
        }
        tmpl <- terra::rast(extent = ex, ncols = ncols, nrows = nrows)
        r <- rasterize(x, tmpl, fun = "count")
        if (!isTRUE(count)) r <- r > 0
        terra::plot(r, ...)
        return(invisible())
    }

    # Vector path: either user opted out of raster (raster = FALSE) or
    # geomtype isn't points (rasterize is points-only).
    sample_callback <- NULL
    if (!is.null(sample_max)) {
        sample_callback <- function(atab) {
            .arrow_sample_max_rows(atab, sample_max)
        }
    }
    sv <- storeRead(x, output = "terra", callback = sample_callback)
    plot(sv, ...)
})

# Register sedonadb_dataframe as a known S3 class so S4 plot dispatch finds it.
setOldClass("sedonadb_dataframe")

#' @rdname plot
#' @param x A `sedonadb_dataframe` from `storeRead(output = "sedona")`.
#' @param values `character` (optional). Column name(s) to include as
#'   SpatVector attributes for coloring. Mirrors terra's `plot(sv, y)`.
#' @param n `integer`. Target number of rows to display (default from
#'   `getOption("giottodisk.plot_sample_max")`). Implemented as systematic
#'   stride sampling (`row_index %% k == 0`) after a `COUNT(*)` pass —
#'   proportional coverage across write order.
#' @export
setMethod("plot", signature("sedonadb_dataframe", "missing"),
    function(x, values = NULL,
        n = getOption("giottodisk.plot_sample_max", 1e5), ...) {

    cols_sql <- if (!is.null(values)) {
        paste(c('"geom"', sprintf('"%s"', values)), collapse = ", ")
    } else {
        "geom"
    }
    ref <- sd_view_ref(x)
    total <- sedonadb::sd_collect(
        sedonadb::sd_sql(sprintf("SELECT COUNT(*) AS n FROM %s", ref))
    )$n
    k <- max(1L, as.integer(ceiling(total / n)))
    sample_sql <- sprintf(
        "SELECT %s FROM %s WHERE (row_index - 1) %% %d = 0",
        cols_sql, ref, k)
    df <- sedonadb::sd_collect(sedonadb::sd_sql(sample_sql))
    wkb <- unclass(wk::as_wkb(df$geom))
    df$geom <- NULL
    sv <- terra::vect(wkb)
    if (ncol(df) > 0L) terra::values(sv) <- df
    if (!is.null(values)) plot(sv, values, ...) else plot(sv, ...)
})
