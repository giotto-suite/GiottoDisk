# docs ----------------------------------------------------------- #
#' @title Get the area of individual polygons
#' @name expanse
#' @description Compute the area covered by polygons in a `parquetGeomStore`
#' or `parquetGeomTileStore`. The tiled terra path parallelizes across tiles
#' via `tileApply`; the sedona path issues a single `ST_Area` query regardless
#' of tiling.
#' @param x `parquetGeomStore` or `parquetGeomTileStore`
#' @param output one of `"data.table"` (default), `"named"`, or `"vector"`.
#' `"data.table"` returns a `data.table` with columns `cell_ID` and `area`.
#' `"named"` returns a named numeric vector with `cell_ID` as names.
#' `"vector"` returns a plain unnamed numeric vector — polygon identity is not
#' preserved.
#' @param engine one of `"terra"` (default) or `"sedona"`. `"sedona"` issues a
#' single `ST_Area(geom)` SQL query via SedonaDB/DataFusion without
#' materializing geometry in R. Pending affine transforms are applied via
#' `ST_Affine` before area computation, so results are correct on lazily
#' scaled or transformed stores. Requires the `sedonadb` package.
#' @param poly_id_col name of the polygon ID column (default `"poly_ID"`)
#' @inheritDotParams terra::expanse
#' @returns depends on `output`: a `data.table`, named `numeric`, or `numeric`
#' @examples
#' # create simple polygons
#' polys <- terra::vect(
#'     c("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
#'       "POLYGON ((1 0, 3 0, 3 2, 1 2, 1 0))",
#'       "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))")
#' )
#' terra::values(polys) <- data.frame(poly_ID = c("a", "b", "c"))
#'
#' store <- storeCreate(type = "parquetgeom")
#' store <- storeWrite(store, polys)
#'
#' # data.table output (default)
#' expanse(store)
#'
#' # named vector
#' expanse(store, output = "named")
#'
#' # sedona engine (requires sedonadb)
#' \dontrun{
#' if (requireNamespace("sedonadb", quietly = TRUE)) {
#'     expanse(store, engine = "sedona")
#' }
#' }
NULL
# ---------------------------------------------------------------- #

#' @rdname expanse
#' @export
setMethod("expanse", signature("parquetGeomStore"), function(
    x,
    output = c("data.table", "named", "vector"),
    engine = c("terra", "sedona"),
    poly_id_col = "poly_ID",
    ...
) {
    output <- match.arg(output)
    engine <- match.arg(engine)
    if (output == "vector") {
        warning("[expanse] output = 'vector' does not preserve polygon IDs",
            call. = FALSE)
    }
    if (engine == "sedona") {
        return(.expanse_sedona(x, output = output, poly_id_col = poly_id_col))
    }
    sv <- storeRead(x, output = "terra")
    terra::crs(sv) <- "local"
    args <- list(sv, ...)
    args$transform <- args$transform %null% FALSE
    areas <- do.call(terra::expanse, args)
    .format_expanse(areas, ids = terra::values(sv)[[poly_id_col]], output = output)
})

#' @rdname expanse
#' @export
setMethod("expanse", signature("parquetGeomTileStore"), function(
    x,
    output = c("data.table", "named", "vector"),
    engine = c("terra", "sedona"),
    poly_id_col = "poly_ID",
    ...
) {
    output <- match.arg(output)
    engine <- match.arg(engine)
    if (output == "vector") {
        warning("[expanse] output = 'vector' does not preserve polygon IDs",
            call. = FALSE)
    }
    if (engine == "sedona") {
        return(.expanse_sedona(x, output = output, poly_id_col = poly_id_col))
    }
    dots <- list(...)
    tile_results <- tilework::tileApply(
        x,
        tiles = x@tiles,
        FUN = function(sv, .I) {
            if (is.null(sv) || nrow(sv) == 0L) return(NULL)
            terra::crs(sv) <- "local"
            tile_args <- c(list(sv), dots)
            tile_args$transform <- tile_args$transform %null% FALSE
            areas <- do.call(terra::expanse, tile_args)
            data.table::data.table(
                cell_ID = terra::values(sv)[[poly_id_col]],
                area = areas
            )
        },
        get_params_x = list(output = "terra")
    )
    result <- data.table::rbindlist(tile_results)
    .format_expanse(result$area, ids = result$cell_ID, output = output)
})

# ---------------------------------------------------------------- #

.expanse_sedona <- function(x, output, poly_id_col) {
    GiottoUtils::package_check("sedonadb",
        repository = "github:apache/sedona-db/r/sedonadb")
    sdf <- storeRead(x, output = "sedona")
    # Register the full query (including any pending ST_Affine) as a new view
    # so that ST_Area operates on the transformed geometry.
    # sd_view_ref() only points to the base parquet view.
    exp_view <- tolower(paste0("gd_exp_", .make_uid()))
    sedonadb::sd_to_view(sdf, exp_view, overwrite = TRUE)
    sql <- sprintf(
        'SELECT "%s" AS cell_ID, ST_Area(geom) AS area FROM "%s"',
        poly_id_col, exp_view
    )
    result <- data.table::as.data.table(
        sedonadb::sd_collect(sedonadb::sd_sql(sql))
    )
    # DataFusion lowercases all identifiers on collection
    .format_expanse(result$area, ids = result$cell_id, output = output)
}

.format_expanse <- function(areas, ids, output) {
    switch(output,
        vector = areas,
        named  = stats::setNames(areas, ids),
        data.table = data.table::data.table(cell_ID = ids, area = areas)
    )
}
