# Write to a Parquet Geometry Tiled Storage

Write geometry data to hive-partitoned spatially tiled parquets. The
partition column is `"tile_index"`.

### Accepted `data` inputs

- `queryableStore`-inheriting - Any dplyr lazy queryable store. Requires
  selection of `type` param.

  - `"points"` - written directly.

  - `"polygons"` - must already be written to
    `"parquetStore"`-inheriting class. Computes geometries from vertices
    which requires row ordering.

- `parquetGeomStore`-inheriting - Simple write to a new tiled store.

## Usage

``` r
# S4 method for class 'parquetGeomTileStore,queryableStore'
storeWrite(
  store,
  data,
  threshold = NULL,
  tiles = NULL,
  i = NULL,
  id_col,
  sdimx,
  sdimy,
  type = c("points", "polygons"),
  contiguous = TRUE,
  dry_run = FALSE,
  write_param = list(),
  verbose = NULL,
  ...
)

# S4 method for class 'parquetGeomTileStore,parquetStore'
storeWrite(
  store,
  data,
  threshold = NULL,
  tiles = NULL,
  i = NULL,
  id_col,
  sdimx,
  sdimy,
  part_col = NULL,
  group_col = id_col,
  type = c("points", "polygons"),
  contiguous = TRUE,
  dry_run = FALSE,
  split_geom = FALSE,
  split_geom_fmt = "poly_%d",
  split_geom_sourcename = NULL,
  flip_vertical = FALSE,
  write_param = list(),
  verbose = NULL,
  ...
)

# S4 method for class 'parquetGeomTileStore,parquetGeomStore'
storeWrite(
  store,
  data,
  threshold = NULL,
  tiles = NULL,
  i = NULL,
  contiguous = TRUE,
  dry_run = FALSE,
  write_param = list(),
  verbose = NULL,
  ...
)
```

## Arguments

- store:

  `dataStore` inheriting class

- data:

  data to write

- threshold:

  `numeric` or `NULL`. Maximum number of records per tile. `NULL`
  (default) auto-selects based on total row count and `type`:

  - `"points"`: 1M/tile up to 1B rows; clamped ramp `n/1000` (1M–100M)
    above.

  - `"polygons"` (vertex rows): 500k/tile up to 500M rows; clamped ramp
    `n/1000` (500k–5M) above that.

- tiles:

  `freeTilePlan`, `tilePlan`, or `NULL`. When a `freeTilePlan` is
  provided (e.g. the output of a `dry_run`), it is used directly without
  replanning. When a `tilePlan` is provided, it is used as the seed grid
  for
  [`tilework::quadtreePlan()`](https://drieslab.github.io/tilework/reference/quadtreePlan.html).
  When `NULL`, a default seed grid is built from the data extent's
  aspect ratio.

- i:

  `integerlike` specific tile(s) to write. Leave as NULL (default) to
  write all.

- id_col:

  `character`. Column name in `data` that holds the geometry identifier.
  For polygons this groups vertices into individual geometries; for
  points this is the feature name. Regardless of the input name, the
  column is written to disk as `poly_ID` (polygons) or `feat_ID`
  (points) so downstream consumers (e.g.
  [`calculateOverlap()`](https://giotto-suite.github.io/GiottoDisk/reference/calculateOverlap.md))
  can assume standardized names.

- sdimx, sdimy:

  `character`. Names of columns containing x and y spatial coordinate
  values, respectively. Used for spatial tile planning and geometry
  construction.

- type:

  `character`. Either `"points"` or `"polygons"`. What type of geometry
  to write.

- contiguous:

  `logical` (default = TRUE) When `TRUE`, tile inclusivity rules (see
  [getBoundedData](https://giotto-suite.github.io/GiottoDisk/reference/getBoundedData.md))
  prevent duplication of data at tile boundaries (assuming no padding).

- dry_run:

  `logical` (default = FALSE). When `TRUE`, runs tile planning and plots
  the adaptive tile layout with record counts, then returns the
  `freeTilePlan` invisibly without writing for inspection. This can then
  be passed to `tiles` to skip tile planning.

- write_param:

  named `list` (optional). Additional params to pass to the writing
  function.

- ...:

  additional params to pass to `tileApply`

- part_col:

  `character` (optional, polygons only). Column name for multi-part
  polygon grouping. See
  [GiottoClass::createGiottoPolygon](https://giotto-suite.github.io/GiottoClass/reference/createGiottoPolygon.html).

- group_col:

  `character`. If `type = "polygons"`, name of column used for spatial
  tile planning via envelope centroids. Defaults to `id_col`.

- split_geom:

  `logical` (polygons only, default `FALSE`). When `TRUE`, multipart
  polygons are split into single-part polygons via
  [`GiottoClass::splitGeom()`](https://giotto-suite.github.io/GiottoClass/reference/combine_split_geoms.html)
  before write. Caller is responsible for producing globally unique IDs
  across all writes; tile orchestration layers (the
  `parquetGeomTileStore` methods) bake the tile index into
  `split_geom_fmt` before dispatch.

- split_geom_fmt:

  `character`. `sprintf` format with a single `%d` slot for the row
  counter. Default `"poly_%d"`. See
  [`splitGeom()`](https://giotto-suite.github.io/GiottoClass/reference/combine_split_geoms.html).

- split_geom_sourcename:

  `character` or `NULL`. Column name carrying the original `poly_ID`
  values (duplicated per part of a multipart polygon). Default `NULL`.

## See also

Other storeWrite methods:
[`storeWrite,parquetEdgeStore,edgeInput-method`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md),
[`storeWrite-parquetGeomStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomStore.md),
[`storeWrite-parquetStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetStore.md)
