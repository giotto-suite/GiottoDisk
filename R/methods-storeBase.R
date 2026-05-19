#' @include class-dataStore.R class-parquetStore.R class-parquetExprStore.R
NULL

# storeBase ####
# Returns the identity-bearing projection of a store: the same object
# stripped of all view-state (lazy @ops, crop/window/post_ops, subset
# @cell_idx / @gene_idx, tile_filter, etc.) and of @params (treated as
# opaque scratch space — never identity-bearing). The result represents
# what a fresh-from-disk handle to the same on-disk artifact would be.
#
# Use sites: .hash, sourceAdopt, sourceContains, snapshotSave artifact
# tagging — anywhere "is this the same on-disk thing" matters.
#
# Dispatch is on concrete classes only — the parquetBase / parquetGeomBase
# mixins are pure-virtual (extend "VIRTUAL", not "dataStore"), so a
# method on them would have nowhere to callNextMethod() to. Each concrete
# class chains through its own contains() hierarchy.

#' @name storeBase
#' @title Identity-bearing projection of a store
#' @description
#' Returns the store with all view-state (lazy ops, crop / window, subset
#' indices, tile filters) and `@params` reset to their prototypes. The
#' result represents the same store as if it were freshly opened from
#' disk: no pending filters, no cached metadata. Used by `.hash` and the
#' artifact-tracking machinery (sourceAdopt, snapshotSave) so two
#' references to the same on-disk thing always hash identically,
#' regardless of any lazy filters applied at call time.
#' @param x a `dataStore` inheriting object
#' @param ... unused
#' @return the same object class with view-state slots zeroed
#' @export
NULL

# Default — identity. Concrete classes call up to here via the
# fileStore chain.
setMethod("storeBase", "dataStore", function(x, ...) x)

# fileStore strips @params.
setMethod("storeBase", "fileStore", function(x, ...) {
    x <- callNextMethod()
    x@params <- list()
    x
})

# parquetStore (queryableStore -> fileStore): adds @ops strip.
setMethod("storeBase", "parquetStore", function(x, ...) {
    x <- callNextMethod()
    x@ops <- list()
    x
})

# parquetGeomStore (parquetStore + parquetGeomBase mixin): adds spatial
# view-state strip on top of parquetStore.
setMethod("storeBase", "parquetGeomStore", function(x, ...) {
    x <- callNextMethod()
    x@crop     <- numeric(0L)
    x@window   <- numeric(0L)
    x@post_ops <- list()
    x
})

# parquetGeomTileStore: adds @tile_filter strip.
setMethod("storeBase", "parquetGeomTileStore", function(x, ...) {
    x <- callNextMethod()
    x@tile_filter <- integer(0L)
    x
})

# parquetExprStore (queryableStore -> fileStore): strips subset state.
# @cell_ids / @feat_ids are left narrowed — they're not read by the
# read_fun (which only consults @path), so they don't affect hashing.
setMethod("storeBase", "parquetExprStore", function(x, ...) {
    x <- callNextMethod()
    x@cell_idx <- integer(0L)
    x@gene_idx <- integer(0L)
    x
})

# unionParquetStore extends parquetBase only (not fileStore), so no
# callNextMethod chain to fileStore. Handle @params + @ops directly and
# recurse storeBase into substores so their view-state is also stripped.
setMethod("storeBase", "unionParquetStore", function(x, ...) {
    x@ops    <- list()
    x@params <- list()
    x@stores <- lapply(x@stores, storeBase)
    x
})

# unionParquetGeomStore extends c("unionParquetStore", "parquetGeomBase").
setMethod("storeBase", "unionParquetGeomStore", function(x, ...) {
    x <- callNextMethod()
    x@crop     <- numeric(0L)
    x@window   <- numeric(0L)
    x@post_ops <- list()
    x
})
