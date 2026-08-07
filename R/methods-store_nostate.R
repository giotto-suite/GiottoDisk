#' @include class-dataStore.R class-parquetStore.R class-parquetExprStore.R
NULL

# .store_nostate ####
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

# Internal — returns the store with all view-state (lazy ops, crop /
# window, subset indices, tile filters) and @params reset to their
# prototypes. Used by .hash so two references to the same on-disk thing
# always hash identically regardless of any pending filters.
setGeneric(".store_nostate",
    function(x, ...) standardGeneric(".store_nostate")
)

# Default — identity. Concrete classes call up to here via the
# fileStore chain.
setMethod(".store_nostate", "dataStore", function(x, ...) x)

# fileStore strips @params.
setMethod(".store_nostate", "fileStore", function(x, ...) {
    x <- callNextMethod()
    x@params <- list()
    x
})

# parquetStore (queryableStore -> fileStore): adds @ops strip.
setMethod(".store_nostate", "parquetStore", function(x, ...) {
    x <- callNextMethod()
    x@ops <- list()
    x
})

# parquetGeomStore (parquetStore + parquetGeomBase mixin): adds spatial
# view-state strip on top of parquetStore.
setMethod(".store_nostate", "parquetGeomStore", function(x, ...) {
    x <- callNextMethod()
    x@crop     <- numeric(0L)
    x@window   <- numeric(0L)
    x@post_ops <- list()
    x
})

# parquetGeomTileStore: adds @tile_filter strip.
setMethod(".store_nostate", "parquetGeomTileStore", function(x, ...) {
    x <- callNextMethod()
    x@tile_filter <- integer(0L)
    x
})

# parquetExprStore (queryableStore -> fileStore): strips subset state +
# both halves of the op chain (@ops prefix, @post_ops suffix). @cell_ids /
# @feat_ids are left narrowed — they're not read by the read_fun (which
# only consults @path), so they don't affect hashing.
setMethod(".store_nostate", "parquetExprStore", function(x, ...) {
    x <- callNextMethod()
    x@cell_idx <- integer(0L)
    x@gene_idx <- integer(0L)
    x@ops      <- list()
    x@post_ops <- list()
    x
})

# unionParquetStore extends parquetBase only (not fileStore), so no
# callNextMethod chain to fileStore. Handle @params + @ops directly and
# recurse .store_nostate into substores so their view-state is also stripped.
setMethod(".store_nostate", "unionParquetStore", function(x, ...) {
    x@ops    <- list()
    x@params <- list()
    x@stores <- lapply(x@stores, .store_nostate)
    x
})

# unionParquetGeomStore extends c("unionParquetStore", "parquetGeomBase").
setMethod(".store_nostate", "unionParquetGeomStore", function(x, ...) {
    x <- callNextMethod()
    x@crop     <- numeric(0L)
    x@window   <- numeric(0L)
    x@post_ops <- list()
    x
})

# unionParquetExprStore: contains "dataStore" only — strip @params,
# both phase chains, and recurse into substores.
setMethod(".store_nostate", "unionParquetExprStore", function(x, ...) {
    x@params   <- list()
    x@ops      <- list()
    x@post_ops <- list()
    x@stores   <- lapply(x@stores, .store_nostate)
    x
})
