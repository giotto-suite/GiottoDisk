# GiottoDisk

Disk-backed data storage and processing for the Giotto spatial omics suite.
Provides lazy/query-based access to tabular and spatial data in Apache Parquet format.

## Package Structure

```
R/
  AllGenerics.R          # S4 generic definitions
  class-dataStore.R      # fileStore, queryableStore, h5ArrayStore, etc.
  class-parquetStore.R   # parquetStore, parquetGeomStore, unionParquetStore hierarchy
  class-gsource.R        # gsource class
  methods-accessors.R    # [,j] indexing, [i, on] join, colnames, nrow, dim, ext, window
  methods-ops.R          # subset, rowSample, head, tail, crop, window<-; .do_op()
  methods-transforms.R   # affine, spin, rescale, shear, spatShift, t, flip; transform helpers
  methods-storeRead.R    # storeRead; as.data.frame, as.terra
  methods-storeWrite.R   # storeWrite methods
  methods-storeInspect.R # storeExists, storePaths methods
  methods-sourceAdopt.R  # sourceContains, sourceAdopt; .move_path in utils.R
  methods-combine.R      # rbind2 methods
  methods-show.R         # show methods
  methods-specialCols.R  # specialCols — columns enforced on disk per store type
  methods-aggregate.R    # calculateOverlap, overlapToMatrix, overlapPointDisk class
  methods-giotto.R       # createGiottoPoints, createGiottoPolygon for parquetGeomBase
  methods-rasterize.R    # terra::rasterize, terra::centroids for parquetGeomBase
  utils.R                # .dplyr_nrow, .dump_tempfile, etc.
  utils-arrow.R          # .arrow_sample_max_rows, .dplyr_ext, .dplyr_crop, etc.
  utils-spatial.R        # spatial helpers
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

parquetBase (VIRTUAL)                 # @fields, @ops — shared by parquet + union
parquetGeomBase (VIRTUAL)             # @window, @crop, @geomtype

unionParquetStore                     # extends "parquetBase" — multi-store union
└── unionParquetGeomStore             # extends c("unionParquetStore", "parquetGeomBase")
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
Operations recorded lazily as a list of steps: `filter`, `head`, `tail`, `sample`, `join`.
Applied in order at `storeRead()` time via `.do_op()`. `arrange` deferred to
`.pstore_to_tibble()` to avoid breaking sample/count ops.

`crop` and `window` are NOT ops — they use dedicated slots on `parquetGeomBase`.

`@ops` contains Arrow-executable ops only. `.do_op()` dispatches each type to its dplyr/Arrow equivalent.

`@post_ops` holds R-level post-materialization ops with no Arrow equivalent. Currently one type: `"transform"` (affine2d). Planned: `"geom_filter"` (exact polygon mask; AABB pre-filter goes to `@ops`, exact polygon test in `@post_ops`).

`crop()`/`window<-` may inject a `"filter"` half-plane op into `@ops` when rotation/shear is pending — a plain Arrow filter op injected by `.pgeom_resolve_extent()`, not a post-op.

### Lazy spatial transforms (`@post_ops` `"transform"` type on `parquetGeomBase`)

`affine(store, affine2d)` records the transform lazily. At `storeRead(output = "tibble"/"terra"/"sf")`, `.apply_post_ops()` applies it. `output = "query"` is unaffected.

- `affine2d@affine` is a **3×3 homogeneous matrix** (GiottoClass `initialize()` converts 2×2 input → 3×3). `@anchor` = intrinsic-space extent (numeric(4)). `@translate` = decomposed translation.
- Transform composition: `affine(store, aff2)` after `affine(store, aff1)` auto-collapses — `affine(existing_affine2d, aff2@affine)` dispatches to `affine(affine2d, matrix)` which composes via centroid-shift method.
- Convention: **post-multiply** (`xy_out = xy_in %*% A + t`), matching GiottoClass's `pre_multiply = FALSE` default.
- `@anchor` is a **composition scaffold only** — used by `affine(affine2d, matrix)` for centroid-shift calculation when composing transforms. NOT read when applying a transform to actual coordinates (`affine(SpatVector, affine2d)` uses only `@affine`). Only `spin()` and `rescale()` need the centroid; callers build the `affine2d` with the correct extent set externally before passing to `affine(store, aff)`. `affine(parquetGeomBase, affine2d)` sets the anchor via `.pgeom_ext_intrinsic(x)` only when `y@anchor` is still at the prototype default — avoids redundant scans on composition chains while ensuring a valid pivot on first use.

**crop/window back-projection when transform pending:**
1. Back-project query extent to intrinsic space: `affine(as.polygons(ext(y)), aff, inv = TRUE)`
2. If rotation/shear: append exact parallelogram half-plane filter to `@ops` (Arrow-native: `a*x_index + b*y_index >= c`)
3. AABB of back-projected corners stored in `@crop`/`@window` for display and AABB pre-filter
4. Axis-aligned transforms (scale/translate only): back-projected extent is already a rectangle, no extra filter needed

Utilities: `.affine_has_rotation()`, `.affine_halfplane_expr()`, `.affine_aabb()`, `.apply_post_ops()`, `.aff_linear_2d()` — all in `utils-spatial.R`.
Transform helpers: `.pgeom_pending_transform()`, `.pgeom_get_transform()`, `.pgeom_set_transform()` — in `methods-transforms.R`.

### @crop vs @window (parquetGeomBase)
- `@crop`: permanent composable spatial subset (numeric(4): xmin, xmax, ymin, ymax)
- `@window`: temporary spatial filter (same format, numeric(0) = unset)
- `storeRead` uses `.pstore_active_extent()`: window takes priority over crop
- `show()` display priority: `ext(x, exact = FALSE)` — fast metadata path, no scan

### ext() semantics and fast path

`ext(x, exact = TRUE)` (default): **always scans** coordinate columns with all
Arrow-phase filter ops applied (crop AABB pre-filters, half-plane filters). Reflects
the true data extent including effects of any spatial subset ops.
- **No transform pending**: `.dplyr_ext(q)` — scan `x_index`/`y_index` min/max
- **Transform pending**: `.dplyr_ext_affine(q, aff)` — scan affine-projected coordinates;
  works for all transform types (axis-aligned cross-terms collapse to 0)

`ext(x, exact = FALSE)`: always returns `.pgeom_ext_estimate()` without scanning.
Used for display — `show()` always calls `ext(x, exact = FALSE)`. Half-plane filter
effects are NOT reflected (uses `@crop`/`disk_extent` AABB).

`.pgeom_ext_estimate(x, aff)`: tightest metadata upper bound.
- Intrinsic base: `@crop` (already ⊆ disk_extent) > `disk_extent` > live scan fallback
- Applies `@window` intersection on top
- Projects 4 AABB corners through affine if pending — exact for axis-aligned, AABB
  overestimate for rotation/shear

`.pgeom_ext_intrinsic(x)`: always live scan in intrinsic (on-disk) space, ignoring
any pending transform. Used internally by `crop()`, `window<-`, `affine()` where
exact intrinsic bounds are required.

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
`x[y, on = c(...), nomatch = NULL]` records a join against `y` (`parquetBase`) as a `"join"`
op in `x@ops`. `x` drives the result; `y` provides columns.

- `on`: named character vector — names are `x` columns, values are `y` columns
- `nomatch = NULL` (inner) stored as string `"inner"` — `NULL` is lost in R lists
- `rbind2`/`storeWrite` blocked while join pending

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
- `"duckdb"`: DuckDB connection
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

`gDirSource` manages a project directory. See `vignettes/articles/design.Rmd` for
full architectural rationale (manifest design, dump lifecycle, deployment, multi-analysis).

```
<project_dir>/
  giottodir.json       # manifest of all tracked artifacts
  artifacts/<uid>/data # vault — artifact subdirs named by uid
  _pending/            # WAL-style pending manifest edits
  giottosave/          # .rds/.qs snapshots of giotto objects
  artifacts/           # vault — also the dump when setArtifactDumpDir(src) is used
```

Artifacts written via `sourceWrite()`: allocate uid → write → hash → queue pending edit.
Pending edits consolidated into `giottodir.json` on-read (optionally locked via `filelock`).
Consolidation is atomic: `.json_atomic_write` writes to temp then renames.

Manifest access: `src["uid"]`, `src["uid", "field"]`, `src[, "field"]`, `as.data.frame(src)`.

### sourceContains / sourceAdopt

`sourceContains(src, store)`: uid-based for `fileStore`; path-prefix for `SpatRaster`;
all-substores for union; path-prefix over all leaves for `IterableMatrix`.

`sourceAdopt(src, store)`: moves unmanaged store into vault and registers it.
- `fileStore`: preserves `@uid`, moves files
- `SpatRaster`: fresh uid; `file.copy` for on-disk, COG write for in-memory
- Union: delegates per substore
- `IterableMatrix` (ANY fallback, soft dep): uses `.im_map_leaves` to adopt each leaf

**IterableMatrix leaf traversal** (`methods-sourceAdopt.R`):
- `.im_leaf_dirs(x)`: all leaf `@dir` paths — for vault checks and hashing
- `.im_map_leaves(x, f)`: applies `f` to each leaf, rebuilds compound structure

**External path guard**: adoption proceeds only for dump-resident paths or when
`giottodisk.adopt_external = TRUE`.

**Session cache** (`.adopt_session_map`): records old → new path per leaf within
one `snapshotSave`. Handles shared-leaf compounds (`raw`/`normalized` same `@dir`).
Reset via `.adopt_session_reset()` at the start of each `snapshotSave`.

**Vault-resident unregistered**: uid extracted via `basename(dirname(path))`;
`src[uid]` checked; manifest entry written if missing — no file movement.

### snapshotSave lifecycle

1. `.adopt_session_reset()`
2. Register external images (`SpatRaster` not in vault)
3. Register external expression matrices (`IterableMatrix`/`HDF5Array` not in vault)
4. Write snapshot atomically: temp file → `file.rename` to `.rds`/`.qs`
5. Tag artifact uids: hash-based detection via `.ss_hash_expr_base` (strips lazy ops)

`sourcePrune`: removes untagged artifacts (protection transitive via `depends` BFS walk).

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

Future plans are documented in `vignettes/articles/roadmap.Rmd`. Key items:

- **`parquetMutableStore`**: disk-backed metadata with main file + named sidecar files per
  column group. Sidecars joined lazily on `[, col]`; consolidated at `snapshotSave()`.
  Join key: `(source_id, row_index)` — stable for the lifetime of the main file.
- **Partition hardlink utility**: creates hardlinks to existing parquet files under a new
  hive partition layout (e.g. `sample_id=s1/`) without copying data. Enables partition
  pushdown for logical groupings (e.g. per-sample blanket ops) added after initial write.

## Conventions
- Internal helpers: `.` prefix, snake_case
- `parquetStore` internals: `.pstore_*` / `parquetBase` internals: `.pbase_*`
- Arrow utilities: `.arrow_*`, `.dplyr_*` / gsource internals: `.gdsrc_*`
- Public store verbs: `storeCreate`, `storeRead`, `storeWrite`, `storeExists`, `storePaths`
- Public source verbs: `sourceWrite`, `sourceContains`, `sourceAdopt`, `sourcePrune`
