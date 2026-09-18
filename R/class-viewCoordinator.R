#' @include pkg_imports.R
NULL

# =============================================================================
# parquetCoordinator — view + space coordinator for parquet-backed subobjects
# =============================================================================
#
# `parquetCoordinator` brokers view + space recipe resolution against
# gobjects whose subobjects use `parquetStore`-inheriting backings. It
# inherits from GiottoClass's `dataTableCoordinator` so that in-memory
# subobjects within an otherwise-backed gobject naturally fall through to
# the in-memory dispatch — `parquetCoordinator IS-A dataTableCoordinator`
# with the added ability to push recipe steps onto store lazy-op queues.
#
# A step resolves to a surviving cell_ID set, queued on the target store
# as a single `id_filter`. Every cell-keyed target — cell metadata,
# spatial locations, expression, and a cell-polygon store — narrows by
# that same set, so a recipe cannot mean different things depending on
# which slot it is read through. `geom` selects the geometry representing
# a cell; `engine` selects the evaluator; neither is a function of the
# target's storage kind. **adr/0015** has the argument, the alternatives,
# and the cost.
#
# Only the CELL axis is modelled. Features and subcellular points are
# separate axes and are not resolved yet, so a crop reaching a transcript
# points store is applied as a geometric clip on its own geometry,
# matching the in-memory path. That is a missing axis, not a property of
# points — adr/0015 "Consequences".
#
# Cross-storage bridging (an in-memory column needed to filter a backed
# store) is handled in `.narrow_store_by_predicate` / `.narrow_dt_via_arrow`:
# a backed owner narrows lazily and joins, an in-memory owner narrows
# eagerly and inlines its ids.
#
# Why inheritance from dataTableCoordinator: a mixed gobject (some slots
# backed, some in-memory) hits the parquetCoordinator's methods for backed
# subobjects and inherits dataTableCoordinator's methods for in-memory
# ones automatically. Same coordinator instance handles both.
#
# See `R/methods-resolveSubobject.R` for the resolveSubobject methods
# and the `defaultViewCoordinator` dispatch registration on `gsource`.
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
