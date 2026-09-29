#' @include pkg_imports.R
NULL

# =============================================================================
# parquetCoordinator — view + space coordinator for parquet-backed subobjects
# =============================================================================
#
# Extends GiottoClass's `dataTableCoordinator`, so in-memory subobjects in an
# otherwise-backed gobject fall through to the in-memory leaves.
#
# A view resolves to one surviving cell_ID set, queued on every backed
# cell-keyed target as a single `id_filter`, so a recipe means the same thing
# whichever slot it is read through. `geom` picks the geometry standing for a
# cell and `engine` the evaluator; neither depends on storage kind
# (adr/0015). Only the cell axis is modelled: a crop reaching a transcript
# points store clips its own geometry, as in memory.
#
# See `R/methods-resolveSubobject.R` for the methods.
# =============================================================================


#' @title parquetCoordinator
#' @description View + space recipe coordinator for gobjects whose
#' subobjects use `parquetStore`-inheriting backings.
#'
#' A recipe step resolves to a surviving `cell_ID` set, queued on the
#' target store as a single `id_filter`. Every cell-keyed target — cell
#' metadata, spatial locations, expression, and a cell-polygon store —
#' narrows by that same set, so a recipe cannot mean different things
#' depending on which slot it is read through.
#'
#' Which engine evaluates a spatial predicate is set by `engine` on
#' [spatRelate()] (or the `giottodisk.spatial_query_engine` option), and
#' is independent of which carrier delivers the rows — it is not a
#' function of `storeRead(output =)`.
#'
#' This version resolves the **cell** axis. Feature subsets and
#' subcellular-point subsets are separate axes and are not modelled yet;
#' a crop reaching a transcript points store is applied as a geometric
#' clip, matching the in-memory path.
#'
#' Inherits from [GiottoClass::dataTableCoordinator-class] so that in-memory
#' subobjects within an otherwise-backed gobject fall through to the
#' in-memory dispatch naturally.
#'
#' @returns `parquetCoordinator`
#' @examples
#' parquetCoordinator()
#' @export
#' @exportClass parquetCoordinator
setClass("parquetCoordinator",
    contains = "dataTableCoordinator")

#' @rdname parquetCoordinator-class
#' @export
parquetCoordinator <- function() new("parquetCoordinator")
