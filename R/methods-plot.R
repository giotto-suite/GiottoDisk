#' @name plot
#' @title Visualize a Store
#' @description
#' Plot a store's contents. Currently only for geometry
#' stores.
#' @param sample_max `integer`-like. Maximum number of geometries to plot
#'   (regularly sampled). If `NULL`, no sampling is performed. Defaults to
#'   `getOption("giottodisk.plot_sample_max")`.
#' @export
setMethod("plot", signature("parquetGeomStore", "missing"),
    function(x, sample_max = getOption("giottodisk.plot_sample_max", 1e5), ...) {

    sample_callback <- NULL
    if (!is.null(sample_max)) {
        sample_callback <- function(atab) {
            .arrow_sample_max_rows(atab, sample_max)
        }
    }

    sv <- storeRead(x,
        output = "terra",
        callback = sample_callback
    )
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
