# 0021. A crop selects by centroid; a cull is an optimizer step and never decides membership

- **Status:** Accepted
- **Date:** 2026-10-27
- **Supersedes:** —
- **Superseded by:** —

## Context

What a spatial crop *means* for a dataset had never been written down, and the
code had drifted into using one meaning as an optimization for another.

**What the code did.** `crop()` / `window()` on a `parquetGeomBase` record an
extent, and every read path applies it as a plain range on the centroid columns
`x_index` / `y_index` (arrow `.dplyr_crop()`, the duckdb/sedona `WHERE`).
Tile stores hold whole cells, each written into the tile containing its centroid
(`getBoundedData()`, inclusive bounds), and `crop()` on a tile store also
narrows `@tile_filter` eagerly. `spatRelate()` pre-culls by calling
`crop(trim_store, <region AABB>)` before its exact predicate.

**The failure that follows.** That pre-cull is centroid membership placed in
front of an `intersects` test, so a polygon that overlaps the region with its
centroid outside the region's bounding box is dropped before the exact test
sees it. Reproduced with a square `[3,7]²` against the region `[0,4]²`:
`terra::relate()` keeps it; the backed store drops it on all three engines
(terra, duckdb, sedona).

**Forces on the meaning of a crop.**

- Vector GIS and spatial SQL default to *intersects* and keep clipping as a
  separate verb; raster tools default to the pixel centre with an opt-in
  `ALL_TOUCHED`. Cells, spots and bins behave like pixels: each is a measurement
  unit whose expression, area and metadata were computed over the whole
  footprint.
- Clipping a cell leaves those attributes stale (the case sf's `st_agr` warns
  about). Intersects membership puts a boundary cell in both neighbouring
  regions, so per-region counts stop summing to the total and tiled work
  double-counts.
- The area centroid commutes with affine transforms, which pending transforms
  in `@post_ops` rely on; a bbox centre does not under rotation. Spatial
  locations derived from polygons are already area centroids, and view crops
  test spatial locations.
- A per-store reach statistic already exists: `params$max_poly_radius`
  (`sqrt(max area / pi)`), written at `storeWrite()`, carried through tiling,
  and used only by `calculateOverlap()`.
- Plotting needs a margin around the viewport so edge cells are drawn and
  clipped. fastplotR built a per-cell extent index for this, keyed on
  `store@ops`, so every distinct view rebuilds it in a full WKB pass.
- A 2M-polygon test (tile-sorted, 5k-row row groups) pruned as well on
  centroids padded by the maximum half-width as on a per-polygon bbox column:
  0.046 s vs 0.057 s, 14 extra cells of ~1340. The bbox column cost +57% file
  size.

## Decision

1. **A crop selects measurement units by membership; the default rule is the
   area centroid in the region.** The point tested is the cell's spatial
   location, of which `x_index` / `y_index` is the store's copy. Geometry of
   cell-keyed units is never altered by a crop or a view. Points and images
   are clipped (ADR 0015's exception). Region-annotation polygons, once they
   carry a role, are clipped.
2. **Other rules are named modes of `crop()`, not separate verbs or flags:**
   `"centroid"` (default), `"intersects"` (whole units touching the region),
   `"clip"` (geometry altered; guarded on cell-keyed units). This replaces the
   view step's `geom = "centroid" | "poly"`.
3. **A cull is an optimizer step.** It returns a superset of the final answer,
   and only the refinement step decides membership. Every mode is one plan:
   cull on centroids in the region's bbox widened by the mode's pad (0 for
   `"centroid"`, cell reach otherwise), then refine (point in region,
   geometry intersects, or intersect-and-clip). No public operation is used as
   a cull.
4. **Planning happens at two levels.** `resolveRecipe()` plans logically:
   which steps, which slots, the surviving ID set once (ADR 0015), one composed
   transform per slot; it records constraints and does no physical work.
   `storeRead()` plans physically, right before execution: tile pruning,
   centroid pushdown, per-mode pads, predicate order. Public store ops append
   constraints instead of mutating `@tile_filter` as they are called.
5. **A consumer may only widen a cull**, e.g. a plot's draw margin. Rows read
   only because of a margin are flagged as context, not selection.
6. **Shared boundaries are half-open** wherever cells are partitioned (tile
   assignment, adjacent regions), so a centroid on an edge lands in one side.

## Consequences

Eager tile pruning is correct under rule 1: tiles hold whole cells by centroid,
so a tile that misses the region holds no qualifying cell. The `spatRelate()`
pre-cull is a bug for `relation = "intersects"` (and for any relation that can
hold for a unit whose centroid is outside the region's bbox) until it uses the
padded cull.

Constraint on future code: do not call a semantic operation (`crop()`,
`window()`, `subset()`) to narrow candidates for another operation. Narrow with
the cull and prove it returns a superset.

`max_poly_radius` is an equal-area radius, not a bound: it under-covers
elongated polygons (a 10×1 cell gives 1.78 against a half-length of 5). Rule 3
is exact by construction only with a recorded maximum half-extent per axis;
until then, the pad is `max_poly_radius * (1 + poly_buf_factor)`, as in
`calculateOverlap()`.

Revisit the store-level statistic if polygons of very different scales share a
store (one outlier widens every read): first a per-tile pad in the tile plan's
metadata, then per-polygon bbox columns.

## Alternatives considered

- **Clip as the default** — invalidates aggregate attributes (counts, area,
  overlaps) of every edge cell; a clipped cell is a new measurement unit, which
  is re-aggregation, not a crop.
- **Intersects as the default** — breaks the partition property, depends on
  segmentation quality, and needs geometry where centroid needs only spatial
  locations.
- **`crop()` meaning clip everywhere, subsets via `spatRelate()`** — aligns with
  terra/sf names, but a dataset-level crop must be membership anyway, so it is a
  large migration that leaves `crop()` on a gobject as a guard rail.
- **A `clip = TRUE/FALSE` argument** — lets a user clip cells by accident, and
  one flag cannot express per-slot behaviour or further membership rules.
- **Per-polygon bbox columns (GeoParquet 1.1 covering)** — same pruning for
  current cells at +57% size, and more state to keep correct under pending
  transforms.
- **Native Parquet `GEOMETRY` for pruning** — helps only engines running a
  spatial predicate; arrow reads, `crop()` and tile selection cannot use it.
- **A per-view extent index (fastplotR's approach)** — rebuilt on every distinct
  view; replaced by a margin on the cull.

## References

- `R/methods-spatRelate.R` — the pre-cull block (`.combined_spatrelate_aabb()`
  then `crop(trim_store, ...)`)
- `R/methods-ops.R` — `crop()` on `parquetGeomBase` / `parquetGeomTileStore`
- `R/methods-storeRead.R` — centroid AABB in the arrow and SQL paths
- `R/tilework.R` — `getBoundedData()`; `R/methods-storeWrite.R` —
  `.pgeom_max_poly_radius()`, `.terra_to_parquet_format()`
- `R/methods-aggregate.R` — `calculateOverlap()` pad
- ADR 0015 (a view step resolves to an ID set)
- GiottoClass `vignettes/articles/design_view_space.Rmd`, "What a crop selects"
- roadmap.Rmd, "Read planning: cull, then refine"
