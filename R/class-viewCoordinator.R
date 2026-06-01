#' @include pkg_imports.R
NULL

# =============================================================================
# parquetCoordinator — view + space coordinator for parquet-backed subobjects
# =============================================================================
#
# DESIGN NOTES (sketch, 2026-05-28)
# ---------------------------------
# `parquetCoordinator` brokers view + space recipe resolution against
# gobjects whose subobjects use `parquetStore`-inheriting backings. It
# inherits from GiottoClass's `dataTableCoordinator` so that in-memory
# subobjects within an otherwise-backed gobject naturally fall through to
# the in-memory dispatch — `parquetCoordinator IS-A dataTableCoordinator`
# with the added ability to push view steps onto store lazy-op queues.
#
# v1 scope (engine-agnostic, pure lazy queue manipulation):
#   * resolveSubobject methods push view steps onto the relevant store
#     using the existing dispatch (`subset()`, `crop()`, pending `affine()`)
#     on `parquetBase`-inheriting classes. No I/O at coordinator time.
#   * Engine choice (arrow / duckdb / sedona) is deferred to `storeRead()`
#     consume time via its `output =` argument.
#   * Cross-storage bridging (in-memory column needed to filter a backed
#     store) is also deferred — when it appears, the parquetize-via-
#     `storeWrite()` happens at storeRead-time, not coordinator-time.
#     Session-scoped temp parquets, prunable via `sourcePrune`.
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
#' subobjects use `parquetStore`-inheriting backings. Pushes recipe steps
#' onto each store's lazy-op slots via the existing `subset()`, `crop()`,
#' and pending `affine()` dispatch on `parquetBase`. No I/O is performed at
#' coordinator time; materialisation happens when the user (or a downstream
#' getter) calls [storeRead()] with their preferred `output =` engine.
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
