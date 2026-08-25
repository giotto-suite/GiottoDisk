# Write to a Parquet Geometry Storage

Write geometry data to a `parquetGeomStore`. Each row is a single
geometry. Enforces the presence of `x_index` and `y_index` cols which
are the centroid (as calculated by the center of the geometry's spatial
extent) and the `row_index` which preserves row ordering.

May be written either directly from terra `SpatVector` or `data.frame`
inputs that are first parsed via GiottoClass `createGiottoPoints` or
`createGiottoPolygon`, depending on the `type` param.

## Usage

``` r
# S4 method for class 'parquetGeomStore,SpatVector'
storeWrite(
  store,
  data,
  meta = NULL,
  row_offset = 0,
  tile_idx = 0L,
  split_geom = FALSE,
  split_geom_fmt = "poly_%d",
  split_geom_sourcename = NULL,
  .arrow_meta = NULL,
  ...
)

# S4 method for class 'parquetGeomStore,data.frame'
storeWrite(
  store,
  data,
  type = c("points", "polygons"),
  id_col,
  sdimx,
  sdimy,
  part_col = NULL,
  row_offset = 0,
  ...
)
```

## Arguments

- store:

  `dataStore` inheriting class

- data:

  data to write

- meta:

  (optional) `data.frame-like` that contains attributes info in the same
  ordering as the spatvector in `data`. If not provided and attributes
  are in the spatvector, they will be automatically extracted with
  [`terra::values()`](https://rspatial.github.io/terra/reference/values.html)

- row_offset:

  `numeric`, will be coerced to `integer` (default = 0L). Offset to
  apply to row number indexing. For example `row_offset = 0L` means the
  first `row_index` value in this file will be `1`.

- tile_idx:

  `integer(1)`. Tile index to write under as a hive partition
  (`tile_index=<n>/`). Defaults to `0L` so flat geom stores share a
  uniform `(tile_index, row_index)` join key with tiled stores without
  needing
  [`dplyr::coalesce`](https://dplyr.tidyverse.org/reference/coalesce.html).
  Pass `NULL` to write without a tile subdir (not recommended for
  `parquetGeomStore`).

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

- .arrow_meta:

  `named list` of `character` strings to add as parquet metadata

- ...:

  addtional params to pass to `parquetStore` method and
  [`arrow::write_dataset()`](https://arrow.apache.org/docs/r/reference/write_dataset.html)

- type:

  `character`. Either `"points"` or `"polygons"`. What type of geometry
  to write.

- id_col:

  `character`. Column name in `data` that holds the geometry identifier.
  For polygons this groups vertices into individual geometries; for
  points this is the feature name. Regardless of the input name, the
  column is written to disk as `poly_ID` (polygons) or `feat_ID`
  (points) so downstream consumers (e.g.
  [`calculateOverlap()`](https://giotto-suite.github.io/GiottoDisk/reference/calculateOverlap.md))
  can assume standardized names.

- sdimx:

  `character`. Column name containing x spatial coordinates.

- sdimy:

  `character`. Column name containing y spatial coordinates.

- part_col:

  `character` (optional, polygons only). Column name for multi-part
  polygon grouping. See
  [GiottoClass::createGiottoPolygon](https://giotto-suite.github.io/GiottoClass/reference/createGiottoPolygon.html).

## See also

Other storeWrite methods:
[`storeWrite,parquetEdgeStore,edgeInput-method`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md),
[`storeWrite-parquetGeomTileStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomTileStore.md),
[`storeWrite-parquetStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetStore.md)
