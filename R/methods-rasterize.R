#' @name rasterize
#' @title Rasterize Parquet Point Data
#' @description
#' Rasterize a `parquetGeomStore` of points onto a template `SpatRaster`.
#' Rather than materializing the full geometry as a `SpatVector`, the
#' template raster's extent and pixel dimensions are used to define spatial
#' bins. Points are assigned to bins in Arrow and aggregated there — only
#' the small summary table is pulled into R.
#'
#' Currently supports point geometries only.
#'
#' @param x `parquetGeomStore` (points)
#' @param y `SpatRaster` template — defines extent, resolution, and CRS.
#' @param field `character` (optional). Column in `x` to aggregate. Not
#'   required when `fun = "count"`.
#' @param fun `character`. Aggregation function. One of `"count"` (default),
#'   `"sum"`, `"mean"`, `"min"`, `"max"`.
#' @param ... additional params passed to [storeRead()]
#' @returns `SpatRaster` with the same extent and CRS as `y`. Cells with no
#'   points are `NA`.
#' @examples
#' # create 500 random points
#' set.seed(42)
#' x <- rnorm(500, 0, 100)
#' y <- rnorm(500, 0, 100)
#' pts <- terra::vect(
#'     data.frame(
#'         x = x,
#'         y = y,
#'         value = x + y),
#'     geom = c("x", "y")
#' )
#'
#' # write to parquetGeomStore
#' store <- storeCreate(type = "parquetgeom")
#' store <- storeWrite(store, pts)
#'
#' # 20 x 20 template raster over the same extent
#' template <- terra::rast(
#'     nrows = 20, ncols = 20,
#'     xmin = -300, xmax = 300, ymin = -300, ymax = 300
#' )
#'
#' # point density per cell
#' r_count <- rasterize(store, template)
#' plot(r_count)
#'
#' # mean of `value` field per cell
#' r_mean <- rasterize(store, template, field = "value", fun = "mean")
#' plot(r_mean)
#' @export
setMethod("rasterize", signature("parquetGeomStore", "SpatRaster"),
    function(x, y, field = NULL, fun = "count", ...) {
    fun <- match.arg(fun, c("count", "sum", "mean", "min", "max"))
    if (fun != "count") {
        checkmate::assert_string(field)
    }

    # template grid parameters
    ext_y  <- terra::ext(y)
    xmin   <- ext_y$xmin
    xmax   <- ext_y$xmax
    ymin   <- ext_y$ymin
    ymax   <- ext_y$ymax
    ncols  <- terra::ncol(y)
    nrows  <- terra::nrow(y)
    xres   <- (xmax - xmin) / ncols
    yres   <- (ymax - ymin) / nrows

    # query, bin, aggregate entirely in Arrow
    # half-open intervals [xmin, xmax) x [ymin, ymax) — no boundary clamping needed:
    #   col_bin in [1, ncols], row_bin in [1, nrows]
    q <- storeRead(x, output = "query", ...) |>
        dplyr::filter(
            x_index >= xmin, x_index < xmax,
            y_index >= ymin, y_index < ymax
        ) |>
        dplyr::mutate(
            col_bin = as.integer(floor((x_index - !!xmin) / !!xres)) + 1L,
            row_bin = !!nrows - as.integer(floor((y_index - !!ymin) / !!yres))
        )

    agg <- if (fun == "count") {
        q |>
            dplyr::group_by(row_bin, col_bin) |>
            dplyr::summarise(value = dplyr::n(), .groups = "drop") |>
            dplyr::collect()
    } else {
        fld <- as.name(field)
        q |>
            dplyr::group_by(row_bin, col_bin) |>
            dplyr::summarise(
                value = switch(fun,
                    sum  = sum(!!fld,  na.rm = TRUE),
                    mean = mean(!!fld, na.rm = TRUE),
                    min  = min(!!fld,  na.rm = TRUE),
                    max  = max(!!fld,  na.rm = TRUE)
                ),
                .groups = "drop"
            ) |>
            dplyr::collect()
    }

    # build matrix from (row, col, value) triplets
    m <- matrix(NA_real_, nrow = nrows, ncol = ncols)
    if (nrow(agg) > 0L) {
        m[cbind(agg$row_bin, agg$col_bin)] <- agg$value
    }

    # wrap as SpatRaster with template georeferencing
    r <- terra::rast(m)
    terra::ext(r) <- ext_y
    crs_y <- terra::crs(y)
    if (nzchar(crs_y)) terra::crs(r) <- crs_y
    r
})
