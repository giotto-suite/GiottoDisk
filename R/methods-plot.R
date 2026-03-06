#' @name plot
#' @title Visualize a Store
#' @description
#' Plot a store's contents. Currently only for geometry
#' stores.
#' @param sample_max `integer`-like (default = 1e4). Maximum
#'   number of geometries to plot (regularly sampled). If
#'   `NULL`, no sampling is performed.
setMethod("plot", signature("parquetGeomStore", "missing"),
    function(x, sample_max = 1e5, ...) {
    
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
