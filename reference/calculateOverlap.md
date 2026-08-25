# Calculate Overlap

Calculate features `y` that are overlapped by polygons `x`. GiottoDisk
provides methods for
[`GiottoClass::calculateOverlap()`](https://giotto-suite.github.io/GiottoClass/reference/calculateOverlap.html)
for operating on disk backed stores.

Both dispatch paths return an `overlapPointDisk` wrapping a flat parquet
store with a unified schema: `poly_ID`, `feat_ID`, any `keep_cols`,
auto-included `count` (if present), `pt_tile_index`, `pt_row_index`,
`row_index`. The result is self-contained for
[`overlapToMatrix()`](https://giotto-suite.github.io/GiottoDisk/reference/overlapToMatrix.md)
and retains integer join keys for optional joins back to the point
store.

When `y` is a `parquetGeomTileStore`, iteration is driven by the point
tiles. Each point belongs to exactly one tile (no deduplication needed).

When `y` is a flat `parquetGeomStore`, iteration is driven by adaptive
polygon tiles built via `quadtreePlan`.

## Usage

``` r
# S4 method for class 'parquetGeomStore,parquetGeomStore'
calculateOverlap(
  x,
  y,
  method = c("vector", "raster"),
  threshold = NULL,
  tiles = NULL,
  pad_y = 500,
  engine = c("terra", "duckdb"),
  path = .dump_tempfile(),
  poly_id_col = "poly_ID",
  feat_id_col = "feat_ID",
  keep_cols = NULL,
  ...
)

# S4 method for class 'parquetGeomStore,parquetGeomTileStore'
calculateOverlap(
  x,
  y,
  method = c("vector", "raster"),
  poly_buf_factor = 0.15,
  pad_y = NULL,
  tile_idx = NULL,
  prune_tiles = FALSE,
  engine = c("terra", "duckdb"),
  path = .dump_tempfile(),
  poly_id_col = "poly_ID",
  feat_id_col = "feat_ID",
  keep_cols = NULL,
  ...
)
```

## Arguments

- x:

  `parquetGeomStore`-inheriting object containing polygons

- y:

  `parquetGeomStore`-inheriting object containing features (points)

- method:

  (`engine = "terra"`) `character`. One of `"vector"` or `"raster"`.
  Method of overlap calculation. See
  [`GiottoClass::calculateOverlap()`](https://giotto-suite.github.io/GiottoClass/reference/calculateOverlap.html)
  for details.

- threshold:

  (optional) `numeric` maximum number of polygons per tile when `y` is a
  flat store. `NULL` auto-computes via `.auto_threshold()`.

- tiles:

  (optional) seed `tilePlan` for quadtree planning when `y` is a flat
  store. `NULL` builds a default grid from the data extent and aspect
  ratio. A `freeTilePlan` (e.g. from `dry_run`) is used directly.

- pad_y:

  (optional) `numeric`. Fixed spatial padding (in data units) around
  each tile extent when fetching points.

  - `parquetGeomStore` dispatch: padding around each polygon tile.
    Default `500`.

  - `parquetGeomTileStore` dispatch: padding around each point tile.
    When `NULL` (default), derived from
    `x@params$max_poly_radius * (1 + poly_buf_factor)`. Supply a value
    directly when `max_poly_radius` was not recorded at write time.

- engine:

  `character` one of `"terra"` or `"duckdb"` (default = `"terra"`).
  `"terra"` iterates over adaptive polygon tiles (or point tiles for
  `parquetGeomTileStore` `y`). `"duckdb"` performs a single full-dataset
  spatial join via DuckDB's spatial extension – tiling params
  (`threshold`, `tiles`, `pad_y`, `poly_buf_factor`) are ignored.
  Requires the DuckDB spatial extension (`INSTALL spatial`).

- path:

  `character` filepath for the result store

- poly_id_col:

  `character` column in `x` holding polygon IDs. Written as `poly_ID` in
  the result (default = `"poly_ID"`).

- feat_id_col:

  `character` column in `y` holding feature IDs. Written as `feat_ID` in
  the result (default = `"feat_ID"`).

- keep_cols:

  `character` (optional) additional columns from `y` to carry into the
  result. `"count"` is auto-included if present in `y`.

- ...:

  additional params to pass

- poly_buf_factor:

  (optional) `numeric`. Fractional buffer beyond
  `x@params$max_poly_radius` (auto-computed at write time) applied
  around each point tile extent when fetching points. Default `0.15`.
  Ignored when `pad_y` is supplied directly.

- tile_idx:

  (optional) `integerlike`. When `y` is a `parquetGeomTileStore`,
  restrict processing to these tile indices only. `NULL` (default)
  processes all tiles. Ignored for `engine = "duckdb"`.

- prune_tiles:

  `logical` (default `FALSE`). When `y` is a `parquetGeomTileStore`, run
  a parallel pre-filter pass to skip point tiles that contain no polygon
  centroids. Enable when polygons are coarsely distributed relative to
  the point tile plan to avoid spawning empty workers. (For example when
  coarsely binning)
