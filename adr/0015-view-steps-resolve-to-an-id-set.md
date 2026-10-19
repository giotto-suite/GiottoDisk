# 0015. A view step resolves to an ID set; the coordinator models the cell axis

- **Status:** Accepted
- **Date:** 2026-09-09
- **Supersedes:** —
- **Superseded by:** —

*(Takes 0015 rather than 0014: `refactor/points-drop-wkb` already holds 0014
locally.)*

## Context

`parquetCoordinator` resolves view + space recipes against gobjects whose
subobjects are `parquetStore`-backed. A crop step declares `geom = "centroid" |
"poly"` — which geometry represents a cell — and the coordinator has to turn
that into an actual narrowing on each subobject.

Two properties of the machinery bear on how:

- `spat_relate` is **already an eager checkpoint, not a lowered predicate**.
  `.spat_relate_narrow()` picks an engine, evaluates against a trimmed store,
  and returns surviving `(source_id, row_index[, tile_index])` as an arrow
  Table, which the carrier `semi_join`s (design.Rmd §"Spatial Predicates").
- `engine` therefore selects *who computes geometry* independently of *who
  carries rows*. It is not a function of `storeRead(output =)`.

The tempting design is per-carrier pushdown: let each store evaluate whatever
part of a recipe it can. A backed cell-polygon store holds the cell polygon, so
a `geom = "poly"` crop looks like free pushdown onto its own `geom` column.

That was built and reversed (`e5cf9db`). It cost two things. The answer landed
in `sr_cache`, a `list()` local to one `storeRead` call, invisible to every
other subobject — so a `materialize()` containing a poly crop evaluated the same
predicate **twice, on two different code paths**, and two paths for one
predicate is the divergence class this design exists to remove. It also masked a
live error: the other path ran GiottoClass's
`spatIDs(spatRelate(polys, region, relation))`, which unwraps a `giottoPolygon`
to its geometry and calls `terra::relate` — no method for a store — so any
tabular target with a poly crop errored on a backed gobject. No test caught it
because no fixture paired a backed polygon source with a tabular target.

Two stale design notes had by then accumulated saying the opposite of the
settled position: `class-viewCoordinator.R`'s 2026-05-28 sketch ("pure lazy
queue manipulation", "No I/O at coordinator time", engine choice "deferred to
`storeRead()` consume time via its `output =` argument"), written before
`engine` existed; and `methods-spatRelate.R`'s header, claiming the arrow
backend "errors loudly" on a spatial predicate and demanding `output =
"sedona"`, written before the checkpoint model.

Separately: there are **three** subset axes — cells (`cell_ID`), features
(`feat_ID`), and subcellular points (transcript id). Only cells is modelled.

## Decision

A recipe step resolves to a surviving `cell_ID` set and is queued on the target
as a single `id_filter`. Every cell-keyed target — cell metadata, spatial
locations, expression, and a cell-polygon store — narrows by that same set.
`geom` selects the geometry representing a cell; `engine` selects the evaluator;
neither is a function of the target's storage kind.

**The coordinator models the cell axis only.** Features and subcellular points
are deferred. A crop reaching a transcript points store is applied as a
geometric clip on the store's own geometry — parity with the in-memory path,
which clips points the same way via `.apply_crops_geometrically()`.

Recorded steps are walked **in order** and intersected, rather than gathered by
type. Order is information a read-time collapse needs, and it is also what
removes the last reason anything wanted to select steps by type — which is why
the coordinator needs no step-filtering accessor from GiottoClass.

## Consequences

One evaluation per crop step per `materialize()`, shared through the resolution
cache, and one usage layer per predicate: a recipe cannot mean different things
depending on which slot it is read through.

Both crop arms reduce to one public expression —
`spatIDs(spatRelate(<carrier>, region, relation))`, with `spatRelate()`
returning the narrowed `cell_ID`s directly on the centroid carrier. There is
nothing left for the coordinator to reach into GiottoClass for: `grep -rc
"GiottoClass:::" R/` is zero, and the only `@` access left is `view@space`.

The cost is real and worth naming. A poly crop now materializes the polygon
source's surviving ids even when the only target is the polygon store itself,
where they could have stayed inside a single read. Correctness bought that;
revisit if a workflow appears that resolves a poly crop against the polygon
store alone and nothing else.

The constraint on future code: do not add a per-carrier pushdown arm to the
coordinator unless its answer lands in the shared resolution cache. A pushdown
whose result is private to one `storeRead` is a second evaluation path by
construction.

The points path is an **exception owed to a missing axis, not a property of
points**. Points can carry tracked transcript IDs and are subject to feature
subsets, so once those axes exist a points crop resolves to an ID set like any
other. `IMPLEMENTATION_viewspace.md` §4 calls the geom path "reserved for
`giottoPoints` (not cell-aggregatable)", which reads as permanent; it is not.
Revisit this ADR when the feature and subcellular axes land.

## Alternatives considered

- **Per-carrier pushdown** — what `e5cf9db` did. Rejected: double evaluation on
  two code paths, a private cache the other subobjects cannot see, and it hid
  the dispatch error above.
- **Three cache slots** (`filter_ids`, `crop_ids:centroid`, `crop_ids:poly`)
  instead of one. That split existed only so the polygon store could take the
  filter arm eagerly while pushing its own crops down lazily. With per-carrier
  pushdown gone, every cell-keyed target consumes the same set, so one
  target-independent slot is both sufficient and the thing that makes a single
  cache safe to share across every subobject in one `materialize()`.
- **A new `relateIDs()` generic in GiottoClass**, giving the poly arm a
  carrier-agnostic ID path. Rejected as unnecessary: `spatIDs(giottoPolygon)`
  already has a non-terra branch (`36ad8eab`), and the gap was only that
  `spatRelate` did not delegate to the carrier.
- **Resolving the poly arm by forcing single columns off the narrowed store**
  (`[` then `as.vector`) rather than adding methods. Rejected: a
  `spatIDs(parquetGeomBase)` / `featIDs(parquetGeomBase)` pair routed through
  the existing `.collect_ids()` is what `spatIDs(giottoPolygon)`'s non-terra
  branch is *for*, it matches `spatIDs(parquetEdgeStore)`, and it keeps the
  ID question one dispatch rather than a column-forcing idiom repeated per
  call site. The methods were added.
- **Positionally indexing the narrowed store** to reuse GiottoClass's in-memory
  `spatRelate` tail verbatim. Rejected by ADR 0001.

## References

- `R/class-viewCoordinator.R`, `R/methods-resolveSubobject.R`
  (`.push_view_to_pstore`, `.surviving_cell_ids_arrow`, `.crop_step_ids`)
- `R/methods-storeRead.R:313-319` — the checkpoint, and `sr_cache`'s scope
- `vignettes/articles/design.Rmd` §"Spatial Predicates (`spat_relate`)"
- ADR 0001 (no positional row indexing), ADR 0012 (one classifier, many
  carriers — the same "do not build a second path" rule, one layer down)
- GiottoClass `IMPLEMENTATION_viewspace.md` §4, `R/methods-IDs.R:167-175`
- `e5cf9db` — the reversed attempt