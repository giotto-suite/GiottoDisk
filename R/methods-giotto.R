# Methods for creating GiottoClass objects from disk-backed stores.
NULL

#' @name createGiottoPolygon
#' @title Create a giottoPolygon from a Parquet Geometry Store
#' @description
#' Create a [GiottoClass::giottoPolygon] object from a `parquetGeomBase`
#' store. The store is kept in the `spatVector` slot -- nothing is
#' materialized. The `unique_ID_cache` is populated via a single Arrow
#' `DISTINCT` query at construction time.
#' @param x `parquetGeomBase` store (polygons)
#' @param name `character`. Polygon set name. Default `"cell"`.
#' @param poly_ID_colname `character`. Column in `x` containing polygon IDs.
#'   Default `"poly_ID"`.
#' @param ... unused
#' @returns `giottoPolygon`
#' @export
setMethod("createGiottoPolygon", signature("parquetGeomBase"),
    function(x, name = "cell", poly_ID_colname = "poly_ID", ...) {
    poly_ids <- .collect_ids(x, poly_ID_colname)
    new("giottoPolygon",
        name                = name,
        spatVector          = x,
        spatVectorCentroids = NULL,
        overlaps            = NULL,
        unique_ID_cache     = poly_ids
    )
})

#' @name createGiottoPoints
#' @title Create a giottoPoints from a Parquet Geometry Store
#' @description
#' Create a [GiottoClass::giottoPoints] object from a `parquetGeomBase`
#' store. The store is kept in the `spatVector` slot -- nothing is
#' materialized. The `unique_ID_cache` is populated via a single Arrow
#' `DISTINCT` query at construction time.
#'
#' When `split_keyword` is provided, a named list of `giottoPoints` is
#' returned, each backed by a lazily filtered subset of the store.
#' @param x `parquetGeomBase` store (points)
#' @param feat_type `character`. Feature type label(s). When using
#'   `split_keyword`, provide one value per split plus one for the remainder.
#' @param feat_ID_colname `character`. Column in `x` containing feature IDs.
#'   Default `"feat_ID"`.
#' @param split_keyword `list` of character vectors of keywords. Each vector
#'   is `grepl()`-matched against feat IDs to define a split. Unmatched
#'   features form the first group. See [GiottoClass::createGiottoPoints()].
#' @param ... unused
#' @returns `giottoPoints`, or named list of `giottoPoints` when
#'   `split_keyword` is provided
#' @export
setMethod("createGiottoPoints", signature("parquetGeomBase"),
    function(x, feat_type = "rna", feat_ID_colname = "feat_ID",
            split_keyword = NULL, ...) {
    feat_ids <- .collect_ids(x, feat_ID_colname)

    if (is.null(split_keyword)) {
        return(new("giottoPoints",
            feat_type       = feat_type,
            spatVector      = x,
            networks        = NULL,
            unique_ID_cache = feat_ids
        ))
    }

    checkmate::assert_list(split_keyword)
    checkmate::assert_character(feat_type, len = length(split_keyword) + 1L)

    split_bools <- lapply(split_keyword, function(keywords) {
        grepl(paste(keywords, collapse = "|"), feat_ids)
    })
    default_bool <- !Reduce("|", split_bools)
    split_bools <- c(list(default_bool), split_bools)
    names(split_bools) <- feat_type

    fld <- as.name(feat_ID_colname)
    mapply(function(bool, ft) {
        ids <- feat_ids[bool]
        store <- x
        store@ops <- c(store@ops, list(list(
            type = "filter",
            expr = call("%in%", fld, ids)
        )))
        new("giottoPoints",
            feat_type       = ft,
            spatVector      = store,
            networks        = NULL,
            unique_ID_cache = ids
        )
    }, split_bools, feat_type, SIMPLIFY = FALSE)
})
