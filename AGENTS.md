# GiottoDisk

Disk-backed data storage and processing for the Giotto spatial omics suite.
Provides lazy/query-based access to tabular and spatial data in Apache Parquet format.

## Documentation map

This file is the entry point for code work — navigation, constraints,
invariants, dispatch patterns. Architectural rationale and user-facing
walkthroughs live in `vignettes/articles/`:

| Doc | Role |
|---|---|
| `AGENTS.md` (this file) | Code navigation, constraints, invariants. Read first when modifying code. |
| `vignettes/articles/design.Rmd` | Architectural rationale: class hierarchy, lazy op system, transform back-projection, manifest design, extent tracking, scale targets. |
| `vignettes/articles/gsource.Rmd` | `gDirSource` walkthrough (the directory-backed `gsource` — the only backend currently shipped): source verbs (`sourceWrite`/`sourceContains`/`sourceAdopt`/`sourcePrune`), snapshot lifecycle, deployment patterns. |
| `vignettes/articles/roadmap.Rmd` | Public-facing direction. Headline items: `parquetMutableStore`, partition hardlink utility, `gSdataSource`. |
| `vignettes/articles/parquetEdgeStore.Rmd` | Edge-store (graph) specifics. |
| `adr/` | Architecture Decision Records: why a choice was made, what was rejected, what it costs. Dated and immutable — read when you are about to change a decision, not to learn current behaviour. |
| `bench/` | Re-runnable regression benchmark (see *Benchmarks* below). Not part of the package — Rbuildignored, results gitignored. |
| `NEWS.md` | User-visible changes per version. Add an entry when you change behaviour, an argument, or an export. |
| `DESCRIPTION` `Remotes:` | Authoritative upstream branch pins. See *Upstream branch pins* below for why each one exists. |

When in doubt, search AGENTS.md first; for deeper "why" follow the
pointers to the relevant vignette. If a decision looks arbitrary and you
are tempted to undo it, check `adr/` before doing so — the alternatives
were often already tried. Write a new ADR when you make a call a future
reader could reasonably reverse, especially one backed by a measurement;
format, criteria and the numbering convention are in `adr/README.md`.
Each ADR carries a code pointer, which is the intended discovery path —
you should meet the relevant ADR by following it from the code.

## Upstream branch pins

GiottoDisk builds against **development branches** of the suite. `Remotes:` in
`DESCRIPTION` is the authoritative, machine-enforced list; this section is the
*why*, and the exit criterion for each pin.

| Package | Branch | Required because | Drop the pin when |
|---|---|---|---|
| `GiottoClass` | `gsource` | `analyzeData`, `reduceData`, `filterData` generics (all in `NAMESPACE` imports) and `labelProportionsParam` (`R/stream-labelProportions.R`) are gsource-only. | those four are exported on `dev`. |
| `Giotto` | `gsource` | The whole param layer dispatched on: `pcaParam` / `autoPcaParam` / `randomPcaParam` / `irlbaPcaParam` / `exactPcaParam`, `varParam`, `covLoessParam`, `covGroupsParam`, `cellStatsParam`, `featStatsParam`, `logNormParam`, `filterParam`. | `suite_dev` exports them. |
| `GiottoUtils` | `dev` | Suite convention; `dev` carries everything used. | `main` catches up. |
| `tilework` | default | Hard `Imports:` dependency, `drieslab/tilework`, not on CRAN. | it ships to CRAN. |

Two traps worth knowing:

- **Giotto has no `dev` branch.** Its general integration line is `suite_dev`,
  which carries only `libraryNormParam` out of the list above. "Just use dev"
  fails differently here than for GiottoClass.
- `gramEigenPcaParam` looks upstream but is **GiottoDisk-owned** (`R/pca-param.R`).
  Do not go looking for it in Giotto.

Not dependencies, so not pinned: `GiottoVisuals` (not imported; arrives via
Giotto), `GiottoData` (one `skip_if_not_installed()` test — worth adding to
`Suggests:` if more tests come to need it). No branch of any package with gmulti,
`@mapping`, `giottoView`, or `giottoSpace` work is required.

**Keeping this current:** when a pin's exit criterion is met, change
`DESCRIPTION` `Remotes:` and this table in the same commit, and add a `NEWS.md`
entry — a loosened build requirement is user-visible. When you add an import
that only exists on a non-default upstream branch, add the symbol to the
"Required because" cell rather than starting a new row.

## Package Structure

```
R/
  AllGenerics.R          # S4 generic definitions
  pkg_imports.R          # roxygen import/export package-level tags
  zzz.R                  # .onLoad / .onAttach hooks
  class-dataStore.R      # fileStore, queryableStore, h5ArrayStore, etc.
  class-parquetStore.R   # parquetBase / parquetGeomBase virtuals; parquetStore /
                         # parquetGeomStore / parquetGeomTileStore; union variants
  class-parquetExprStore.R  # parquetExprStore, unionParquetExprStore (long-format
                            # expression: row_id, col_id, value, source_id)
  class-parquetEdgeStore.R  # parquetEdgeStore (graph edges) + edgeInput hierarchy
  class-fileInputs.R     # exprInput, mtxInput, tenxH5Input, cellbinGefInput, etc.
  class-gsource.R        # gsource, gDirSource
  methods-accessors.R    # [,j] indexing, [i, on] join, colnames, nrow, dim, ext, window
  methods-ops.R          # subset, rowSample, head, tail, unique, crop, window<-; .ptabular_apply_op()
  methods-transforms.R   # affine, spin, rescale, shear, spatShift, t, flip; transform helpers
  methods-spatRelate.R   # spat_relate op + engine dispatch (sedona/duckdb/terra)
  methods-storeRead.R    # storeRead + .pstore_to_sedona / .pstore_to_duckdb /
                         # .pgstore_to_spatial; shared SQL builders
  methods-storeWrite.R   # storeWrite methods
  methods-storeInspect.R # storeExists, storePaths methods
  methods-store_nostate.R # store-state strip for read-only checks
  methods-combine.R      # rbind2 methods (parquet + union variants)
  methods-show.R         # show methods
  methods-specialCols.R  # specialCols — columns enforced on disk per store type
  methods-expanse.R      # expanse() / centroid extraction
  methods-plot.R         # plot methods (parquetGeomBase)
  methods-aggregate.R    # calculateOverlap, overlapToMatrix, overlapPointDisk class
  methods-rasterize.R    # terra::rasterize, terra::centroids for parquetGeomBase
  methods-giotto.R       # createGiottoPoints, createGiottoPolygon for parquetGeomBase
  methods-parquetExprStore.R  # subset / union / storeWrite / generic dispatch
                              # for parquetExprStore
  methods-fileInputs.R   # readers for the *Input file-input types
  methods-sourceWrite.R  # sourceWrite (allocate uid → write → register)
  methods-sourceAdopt.R  # sourceContains, sourceAdopt
  methods-sourcePrune.R  # sourcePrune
  methods-snapshotSave.R # snapshotSave (giottosave snapshot lifecycle)
  methods-snapshotLoad.R # snapshotLoad
  methods-snapshotDelete.R # snapshotDelete
  stream-filter.R        # filterData(parquetExprStore, ...)
  stream-normalize.R     # processData(parquetExprStore, libraryNormParam/logNormParam)
  stream-hvf.R           # processData(parquetExprStore, varParam) HVF selection
  stream-pca.R           # reduceData(parquetExprStore, randomPcaParam)
  stream-qc.R            # processData(parquetExprStore, cellStatsParam/featStatsParam)
  stream-recommend.R     # streaming recommender utilities
  convenience-cosmx.R    # CosMx import convenience
  convenience-stereoseq.R # Stereo-seq import convenience
  convenience-xenium.R   # Xenium import convenience
  utils.R                # .dplyr_nrow, .dump_tempfile, .move_path, etc.
  utils-arrow.R          # .arrow_sample_max_rows, .dplyr_ext, .dplyr_crop, etc.
  utils-spatial.R        # affine half-plane helpers, AABB, etc.
  utils-parquetExprStore.R # .pestore_* helpers (LUT remap, finalize, etc.)
  tilework.R             # tile plan integration
```

## Class Hierarchy

```
dataStore (VIRTUAL)
└── fileStore                         # path + read_fun + uid + params
    └── queryableStore                # dplyr-queryable fileStore
        └── parquetStore              # extends c("queryableStore", "parquetBase")
            └── parquetGeomStore      # extends c("parquetStore", "parquetGeomBase")
                └── parquetGeomTileStore  # adds @tiles (tilePlan) + @tile_filter (integer) slots
        └── parquetExprStore          # long-format triplet expression
                                      # (row_id, col_id, value, source_id)
        └── parquetEdgeStore          # graph edges + @nodes parquetStore sidecar

parquetBase (VIRTUAL)                 # @fields, @ops — shared by parquet + union
parquetGeomBase (VIRTUAL)             # @window, @crop, @geomtype

unionParquetStore                     # extends "parquetBase" — multi-store union
└── unionParquetGeomStore             # extends c("unionParquetStore", "parquetGeomBase")

unionParquetExprStore                 # multi-store union of parquetExprStores
```

`parquetBase` is used for shared dispatch (colnames, [,j], subset, rowSample, nrow, show ops)
without traversing the fileStore chain. Methods on `parquetBase` must be self-contained
and never use `callNextMethod` — storeRead is NOT defined on parquetBase for this reason.

## Standard Column Schema

Enforced via `specialCols()`:
- `parquetBase`: `row_index`, `source_id`
- `parquetGeomBase`: `x_index`, `y_index`, `geom` (WKB), `tile_index`

Special cols are hidden from `colnames()` but available in queries; injected and dropped
during materialization as needed.

`row_index`: intrinsic row ordering within a source
`source_id`: stable sort key across rbind/union stores — hive partition (`source_id=<uid>/`)
`x_index`, `y_index`: float centroid coordinates for fast spatial extent filtering
`geom`: WKB geometry column (GeoParquet spec)
`tile_index`: tile membership — hive partition (`tile_index=<n>/`). Flat `parquetGeomStore`
  writes use `tile_index=000` so joins always work on a uniform `(tile_index, row_index)` key.

## Key Design Decisions

### No row indexing
Positional row indexing not supported — parquet is not row-addressable, and a 2B-row index
costs ~16GB. Use `subset()` for value-based filtering, `rowSample()` for downsampling.

### nrow() returns numeric (double)
Handles counts up to 2^53. Arrow COUNT(*) returns int64 → `as.numeric()` converts cleanly.
Always queries via COUNT(*) — no caching.

### Lazy ops via @ops slot
Operations recorded lazily as a list of steps. User-facing op types:
`filter`, `head`, `tail`, `sample`, `distinct`, `join`, `spat_relate`. Applied
in order at `storeRead()` time. `arrange` deferred to `.pstore_to_tibble()` to
avoid breaking sample/count ops.

`crop` and `window` are NOT ops — they use dedicated slots on `parquetGeomBase`.

Most ops are arrow-evaluable via `.ptabular_apply_op()` (filter / head / tail / sample /
distinct / join / id_filter). Two ops route through specialized handlers:

- **`spat_relate`** is engine-evaluated (sedona / duckdb / terra) — see
  *spat_relate op + engine dispatch* below. `.pbase_storeread_processing()`
  intercepts spat_relate ops in the op loop and calls `.spat_relate_narrow()`,
  which returns an arrow Table of surviving ids; arrow then `semi_join`s the
  main query. Never reaches `.ptabular_apply_op`.
- **`id_filter`** is an internal op type created ephemerally by the spat_relate
  narrow path to cache surviving ids across chained spat_relate ops. Has both
  an arrow handler (`.ptabular_apply_op` "id_filter" → `semi_join`) and a SQL handler
  (`.pstore_sql_inner` "id_filter" → correlated `EXISTS` subquery).

`@post_ops` holds R-level post-materialization ops with no Arrow equivalent.
Currently one type: `"transform"` (affine2d).

`crop()`/`window<-` may inject a `"filter"` half-plane op into `@ops` when
rotation/shear is pending — a plain Arrow filter op injected by
`.pgeom_resolve_extent()`, not a post-op.

### Lazy spatial transforms (`@post_ops` `"transform"` type on `parquetGeomBase`)

`affine(store, affine2d)` records the transform lazily; `.apply_post_ops()`
applies it at `storeRead(output = "tibble"/"terra"/"sf")`. `output = "query"`
is unaffected. Composition auto-collapses to one `affine2d`; convention is
post-multiply (`xy_out = xy_in %*% A + t`). When a transform is pending,
`crop()`/`window<-` back-project the query extent through the inverse affine
into intrinsic space and inject a half-plane filter op for rotation/shear.

Helpers: `.pgeom_pending_transform`, `.affine_halfplane_expr`, `.apply_post_ops`
(in `utils-spatial.R` / `methods-transforms.R`). See **design.Rmd:
Spatial Transforms are Post-Ops** for the affine composition convention,
`@anchor` semantics, and the full back-projection algorithm.

### @crop vs @window (parquetGeomBase)
- `@crop`: permanent composable spatial subset (numeric(4): xmin, xmax, ymin, ymax)
- `@window`: temporary spatial filter (same format, numeric(0) = unset)
- `storeRead` uses `.pstore_active_extent()`: window takes priority over crop
- `show()` display priority: `ext(x, exact = FALSE)` — fast metadata path, no scan

### ext() semantics and fast path

Three helpers, picked by `exact` and call site:

- `ext(x, exact = TRUE)` (default) — **scans** with all Arrow-phase filters
  applied. `.dplyr_ext` for no-transform; `.dplyr_ext_affine` for pending
  transform (scans affine-projected coords).
- `ext(x, exact = FALSE)` — no scan; returns `.pgeom_ext_estimate` (tightest
  metadata bound from `@crop`/`disk_extent` ∩ `@window`, projected through
  any pending affine). `show()` uses this.
- `.pgeom_ext_intrinsic(x)` — internal helper that always live-scans
  intrinsic (on-disk) space, ignoring any pending transform. Used by
  `crop()`, `window<-`, `affine()` where intrinsic bounds are required.

See **design.Rmd: Spatial Extent Tracking** for the layered cache /
runtime / live-query rationale and edge cases.

### Spatial param provenance (`@params` on `parquetGeomBase`)

| Field | Set where | Consumed where |
|-------|-----------|----------------|
| `@params$disk_extent` | `storeWrite(,SpatVector)` → `ext(data)`; tile stores → tile plan bounds | GeoParquet `bbox` metadata; `show()` fast path |
| `@params$crs` | `storeWrite(,SpatVector)` → `terra::crs(data)`; tile stores → propagated from first written tile | `as.terra`/`as.sf` output; GeoParquet metadata |
| `@params$max_poly_radius` | `storeWrite(,SpatVector)` polygons → `sqrt(max_expanse/pi)`; tile stores → max over tiles | `calculateOverlap` `pad_y` |
| `@params$use_xy_as_geom` | `terra::centroids()` → `TRUE` | `storeRead` geom path selector; `as.data.frame(geom="XY")` guard |

`@geomtype` slot — set by: `storeWrite(,SpatVector)`, `storeWrite(TileStore,*)` from `type` param,
`storeWrite(,parquetGeomStore)` (copied), `terra::centroids()` (forced `"points"`).
Consumed by: `as.data.frame(geom="XY")` guard, GeoParquet metadata, `show()`.

### @datatype — on-disk type overrides

`@datatype` is a named list mapping column names to R type strings (`"character"`, `"raw"`,
`"integer"`, etc.). Intended as a mechanism for callers to signal that certain columns need
to be read as a different Arrow type than what is physically in the file — e.g. when a
writer uses a type that Arrow R cannot execute compute kernels on. Callers set `@datatype`
manually before calling `initialize()`. The storeRead hook that applies the override schema
is not yet wired; see roadmap.

### storeRead shared implementation
`parquetStore::storeRead` and `unionParquetStore::storeRead` share post-dataset logic via
`.pbase_storeread_processing()`, splitting after `atab` is obtained.

### terra output: WKB direct path
`output = "terra"` uses `terra::vect(as.list(geom_col))` directly — avoids sf intermediate
(~6x faster at scale). `output = "sf"` still routes through sf.

When `store@params$use_xy_as_geom = TRUE` (set by `terra::centroids()`), the read path builds
point geometries from `x_index`/`y_index` and drops `geom` entirely — no store rewrite.

### storeRead field injection: `.pstore_fields_requested` + `.pstore_lazy_fields`
Dispatch chain changes `output` to `"query"` via `callNextMethod`, so special column
injection must happen at the **topmost** level before the first `callNextMethod`.

- `.pstore_fields_requested(store, fields)` — resolves user fields; no-ops if already expanded
- `.pstore_lazy_fields(store, fields, output)` — expands to include required special cols;
  returns `NULL` if `fields = NULL` (fetch all)
- Required injections: `x_index`/`y_index` always (geom stores); `tile_index` (tileStore);
  `source_id`/`row_index` for materialized outputs; `geom` for `"terra"`/`"sf"`

### [i, on] join op
`x[y, on = c(...)]` records a join against `y` (`parquetBase`) as a `"join"` op
in `x@ops`. `x` drives the result; `y` provides columns.

- `on`: named character vector — names are `x` columns, values are `y` columns
- `nomatch` follows data.table convention: default (or `NA`) = left join (preserve
  `x`, NA fill on miss); `NULL` = inner (drop unmatched). Any other value errors.
  Stored as `"inner"`/`"left"` strings — `NULL` is lost in R lists.
- `rbind2`/`storeWrite` blocked while join pending

### spat_relate op + engine dispatch
`spatRelate(x, y, relation, engine = NULL)` records a `"spat_relate"` op carrying
the predicate (`"intersects"`/`"within"`/etc.), the query geometry (as `y_wkt` or
`y_store`), and an optional per-call `engine`. The op is evaluated at
`storeRead()` time by one of three engines:

- **`sedona`** — DataFusion via `{sedonadb}`; emits `ST_<pred>(geom, ...)` SQL.
- **`duckdb`** — DuckDB spatial extension via `{duckdb}` + `{dbplyr}`; same SQL
  shape, no SRID (DuckDB has no CRS concept).
- **`terra`** — deps-free fallback. Tile stores stream per-tile via
  `tilework::tileApply`; non-tile stores materialize the trim and run
  `terra::relate`.

Engine selection precedence: per-call `engine` arg > option
`giottodisk.spatial_query_engine` > `"auto"` (sedona > duckdb > terra). When
"auto" falls through to terra, `rlang::inform` fires once per session.

**Narrow path** (`.spat_relate_narrow`): builds a trim store of ops *before* this
spat_relate (any earlier spat_relate ops swapped for internal `id_filter` ops
carrying cached surviving ids), runs the engine to collect surviving
`(source_id, row_index[, tile_index])` as an arrow Table, then `semi_join`s the
main arrow query — arrow stays lazy past the spatial step.

**`id_filter` is an internal op type** (`.ptabular_apply_op` "id_filter" + `.pstore_sql_inner`
"id_filter") used only inside the trim store during narrow eval; never reaches
the user's `@ops`. SQL form is a correlated `EXISTS` subquery against a temp
view registered from the cached id arrow Table — DataFusion doesn't support
tuple-IN subqueries.

### subset() NSE
Uses `rlang::enquo()`. `.inline_local_vars()` walks AST inlining non-column symbols.
When calling `subset()` programmatically (e.g. from within `mapply`), S4 dispatch may
fall through to `base::subset.default` — append to `@ops` directly instead:
```r
store@ops <- c(store@ops, list(list(type = "filter", expr = my_call)))
```

### GeoParquet metadata
`.arrow_meta_add_geoparquet()` attaches GeoParquet-compliant `"geo"` metadata. Uses
`store@params$disk_extent` as bounding box, `store@params$crs %||% ""` for CRS.
Tiled stores write the top-level extent to each tile file (intentional — files are internal).

## Output Formats
- `"query"`: Arrow lazy dataset (default)
- `"tibble"`: collected data.table, arranged by source_id/tile_index/row_index
- `"duckdb"`: lazy `tbl_dbi` over a duckdb `TEMP VIEW` of the parquet dataset.
  Native compile path via `.pstore_to_duckdb` — `read_parquet` SQL with
  per-tile UNION ALL, `@ops` translated to WHERE/SELECT/LIMIT/EXISTS, spatial
  extension's `ST_*` for any pending transforms or spat_relate ops. User-
  supplied connection honoured via `duckdb_params$conn`; otherwise an
  ephemeral in-memory connection is created and kept alive by the returned
  `tbl_dbi`.
- `"sedona"`: lazy `sedonadb_dataframe` (DataFusion via sedonadb) — same
  shape as duckdb path: per-tile UNION ALL, `@ops` translated, `ST_*` for
  spatial. Built by `.pstore_to_sedona`. Shares the @ops translation
  builder `.pstore_sql_inner` with the duckdb path.
- `"terra"`: SpatVector (`parquetGeomBase` only)
- `"sf"`: sf object (`parquetGeomBase` only)

`omit_internals = TRUE` (default): strips `row_index`/`source_id` from materialized output.

`tile_idx` on `storeRead(parquetGeomTileStore)`: explicit tile partition filter (integer vector).
When `tile_idx = NULL`, falls through to `@tile_filter` if non-empty. Explicit `tile_idx`
always takes precedence.

`@tile_filter` on `parquetGeomTileStore`: set by `crop()` to the indices of tiles whose
padded bounds intersect `@crop` (computed via `tilework::intersect(x@tiles, ext(x@crop))`).
Applied as an Arrow hive-partition filter at `storeRead` time — eliminates irrelevant tile
directories before any row-level filtering. Overwritten (not accumulated) on each `crop()`
call; `integer(0L)` means no tile pruning. For rotated crops the tile set may be slightly
over-inclusive (AABB overestimate); exact filtering is handled by the half-plane op in `@ops`.

## GiottoClass Integration

### as.data.frame / as.terra
`as.data.frame(parquetBase)` → `storeRead(output = "tibble")` → `as.data.frame()`.
`as.data.frame(parquetGeomBase, geom = "XY")` → tibble with `x`/`y` centroid columns (points only).
`as.terra(parquetGeomBase)` → `storeRead(output = "terra")`.

### unique / as.vector
`unique(parquetBase)` records a `"distinct"` lazy op on `@ops` (deduplicates by all user
columns at read time, ignoring special cols). Compose with `[, j]` for column-scoped distinct.

`as.vector(parquetBase)` materializes a single-column store as `list(col = c(...))` — matches
base R's `as.vector()` on a data.frame, keeping pipelines interchangeable with in-memory objects.
Errors if more than one column is selected.

```r
unique(x[, "feat_ID"]) |> as.vector()  # list(feat_ID = c(...))
```

### createGiottoPoints / createGiottoPolygon
Store is placed directly in the `spatVector` slot (`ANY`) — nothing materialized.
`unique_ID_cache` populated via single Arrow `DISTINCT` query at construction (avoids broken
behavior when GiottoClass code calls `featIDs()`/`spatIDs()`).

`createGiottoPoints` supports `split_keyword`: IDs matched via `grepl` on the cached IDs,
each split gets a lazily filtered store via direct `@ops` append. Cache reused for both
grepl matching and per-split caching — no extra file reads.

`feat_ID_colname` / `poly_ID_colname` default to `"feat_ID"` / `"poly_ID"`. No rename applied.

### terra::rasterize (parquetGeomStore)
Bins points into raster cells using `x_index`/`y_index` — no geometry materialization.
Aggregation (`"count"`, `"sum"`, `"mean"`, `"min"`, `"max"`) in Arrow; only summary table
pulled into R. Returns `SpatRaster` with template extent and CRS.

### terra::centroids (parquetGeomBase)
Sets `@geomtype = "points"` and `@params$use_xy_as_geom = TRUE`. No data read or rewritten.
At read time, `x_index`/`y_index` are used as point coordinates instead of parsing `geom` WKB.

## gsource / Project Directory System

`gDirSource` manages a project directory with vault, manifest, and snapshot
lifecycle. See `vignettes/articles/gsource.Rmd` for the user-facing walkthrough
of the verbs (`sourceWrite`/`sourceContains`/`sourceAdopt`/`sourcePrune`,
`snapshotSave`/`snapshotLoad`/`snapshotDelete`) and deployment patterns;
`vignettes/articles/design.Rmd` "Project Management" for architectural
rationale (manifest design, dump lifecycle, multi-analysis guarantees).

```
<project_dir>/
  giottodir.json       # manifest (consolidated from _pending/ on read)
  artifacts/<uid>/...  # vault
  _pending/            # WAL-style pending manifest edits
  giottosave/          # .rds/.qs snapshots
```

Manifest access: `src["uid"]`, `src["uid", "field"]`, `src[, "field"]`,
`as.data.frame(src)`. Consolidation uses `.json_atomic_write` (temp →
rename); pending writes are also atomic. `filelock` (Suggests) provides
advisory locking for the consolidation step.

### Implementation constraints (for code work)

- **`sourceContains` dispatch**: uid-based for `fileStore`; path-prefix
  for `SpatRaster`; all-substores for union; path-prefix over
  `.im_leaf_dirs(x)` for `IterableMatrix` (ANY fallback).
- **`sourceAdopt` dispatch**: `fileStore` preserves `@uid` + moves
  files; `SpatRaster` fresh uid (`file.copy` on-disk, COG write
  in-memory); union delegates per substore; `IterableMatrix` walks
  leaves via `.im_map_leaves(x, f)`.
- **External-path guard**: adoption only proceeds for dump-resident
  paths or when `giottodisk.adopt_external = TRUE`.
- **Session cache** (`.adopt_session_map`): records old → new path per
  leaf within one `snapshotSave` so shared leaves (`raw`/`normalized`
  with same `@dir`) only move once. Reset at the start of each
  `snapshotSave` via `.adopt_session_reset()`.
- **Vault-resident unregistered**: when `setArtifactDumpDir(src)` is
  active, artifacts land in the vault directly. `sourceAdopt` extracts
  uid via `basename(dirname(path))`, checks `src[uid]`, writes the
  manifest entry if missing — no file movement.
- **`snapshotSave` lifecycle**: session-reset → adopt external images
  → adopt external expression matrices → atomic snapshot write
  (temp → rename) → tag artifacts (`name` appended to `giottosave`
  manifest field; hash-based detection strips lazy `@ops`).
- **`sourcePrune`**: removes untagged artifacts; protection transitive
  via BFS over `depends`.

### Tiled geometry write planning
`storeWrite(parquetGeomTileStore)` uses `tilework::quadtreePlan()`. `dry_run = TRUE` plots
layout and returns `freeTilePlan` for reuse. Auto-threshold via `.auto_threshold(n, type)`.

### parquetGeomStore params
- `.pgeom_max_poly_radius`: max polygon radius for `pad_y` in `calculateOverlap`
- `.gdsrc_allocate_artifact_dir(p, uid, create)`: returns **named** character (name = uid)
- `.move_path(from, to)`: `file.rename` with cross-device fallback

## calculateOverlap / overlapToMatrix

### overlapPointDisk class
Inherits `overlapInfo`. Wraps overlap `parquetStore` with metadata:
- `poly_id_col`, `feat_id_col`: column names for polygon/feature IDs
- `spat_ids`, `feat_ids`: full ID universes (Arrow `DISTINCT` at `calculateOverlap` time)
  — ensures zero-overlap rows/cols appear in the output matrix
- `poly_uids`, `feat_uids`: UID(s) of the polygon and feature source stores captured at
  `calculateOverlap` time — provenance tracking; also intended for future `depends` BFS
  walk in `sourcePrune` (not yet wired). Use `storeUID(overlap)` to access as a named list.

`storeUID()` dispatch:
- `fileStore` → `store@uid` (single string)
- `unionParquetStore` → character vector, one per substore
- `overlapPointDisk` → `list(poly = ..., feat = ...)` named list of UID vectors

### Unified overlap result schema
Both dispatch paths produce a flat `parquetStore`:
- `poly_ID`, `feat_ID` (character) — polygon and feature IDs
- `keep_cols` (optional) — extra columns from point store
- `count` (if present in point store)
- `pt_tile_index`, `pt_row_index` — join keys back to point store

### calculateOverlap dispatches
- `(parquetGeomStore, parquetGeomTileStore)`: `tileApply` with `pad_y` padding
- `(parquetGeomStore, parquetGeomStore)`: `quadtreePlan` + `tileApply`

Both paths return `overlapPointDisk`. Internals return plain `parquetStore`.

### overlapToMatrix
Groups by `(feat_ID, poly_ID)`, counts/sums, builds COO on integer keys (string→int LUTs
joined onto raw batches *before* aggregation — aggregating strings first leaves dangling
utf8_view buffers at scale). IDs sorted via `GiottoUtils::mixedsort()` before COO
construction. `all_feat_ids`/`all_cell_ids` params preserve zero-overlap entries.

**Destination branches on `store_type`** (default: `getOption("giotto.gdsrc_matrix_format", "bpcells")`):
- `"bpcells"` / `"hdf5"`: COO → `.write_overlap_mtx()` (streamed Matrix Market) →
  `.mtx_to_store()` → BPCells or HDF5Array.
- `"parquetexpr"`: COO → `.coo_to_parquetexpr()` directly. Transmutes
  `(i, j, n)` → `(col_id, row_id, value)`, arranges by `row_id` for row-group skipping,
  streams record batches to a `parquetExprStore` — **no MTX intermediate, no dgCMatrix
  materialization**.

## Metadata

Cell/feature metadata is currently stored as in-memory `data.table`. This is sufficient
for typical cell metadata sizes (up to ~1M cells × many cols). At transcript/feature scale
(~50M rows) it becomes heavier but remains workable.

## Roadmap

High-level public-facing direction is in `vignettes/articles/roadmap.Rmd`. Headline items:

- **`parquetMutableStore`**: disk-backed metadata with main file + named sidecar files per
  column group. Sidecars joined lazily on `[, col]`; consolidated at `snapshotSave()`.
  Join key: `(source_id, row_index)` — stable for the lifetime of the main file.
- **Partition hardlink utility**: creates hardlinks to existing parquet files under a new
  hive partition layout (e.g. `sample_id=s1/`) without copying data. Enables partition
  pushdown for logical groupings (e.g. per-sample blanket ops) added after initial write.
- **`gSdataSource` (SpatialData adapter)**: read-only adapter over a published sdata Zarr
  store. Implementation notes when starting work:
  - sdata points/shapes are GeoParquet — `parquetGeomStore` already reads/writes
    GeoParquet (`.arrow_meta_add_geoparquet()`). The points/shapes adapter is likely a
    thin schema-conformance + column-rename layer, not a new store class.
  - sdata coordinate transforms are 2D affine (translate/scale/rotate/affine). Map them
    onto `parquetGeomBase @post_ops "transform"` via `affine2d` — the existing
    composition path covers it; no new transform machinery needed.
  - sdata tables (AnnData/Zarr) and OME-NGFF image pyramids are the heavier mappings;
    points/shapes are the cheap wins to land first.

## Benchmarks

`bench/` is a re-runnable regression benchmark, not part of the package
(Rbuildignored; `bench/results/` gitignored).

```sh
Rscript bench/regress.R                        # every dataset, every case
Rscript bench/regress.R --cases=PCA --reps=5   # one slice (regex)
Rscript bench/regress.R --ref=A --now=B        # pin both ends
```

It reports the **ratio** between two trees measured back to back, never absolute
times — **below 1.00 means slower than the ref**. `--ref` compares against your
working tree; `--now` replaces the working-tree side with a second worktree, so
pin both ends for any number that lands in an issue, commit message or ADR,
since otherwise the comparison drifts as the branch moves. `--data=atera` needs
`GD_BENCH_H5` pointing at a `cell_feature_matrix.h5`; without it the default run
does synthetic only and says so. Reasoning and the full flag list are in
`bench/README.md`.

## Conventions
- Internal helpers: `.` prefix, snake_case
- `parquetStore` internals: `.pstore_*` / `parquetBase` internals: `.pbase_*`
- Arrow utilities: `.arrow_*`, `.dplyr_*` / gsource internals: `.gdsrc_*`
- Public store verbs: `storeCreate`, `storeRead`, `storeWrite`, `storeExists`, `storePaths`
- Public source verbs: `sourceWrite`, `sourceContains`, `sourceAdopt`, `sourcePrune`
