# Read a `dataStore`

Read from a `dataStore` inheriting object. The output should be a useful
representation of the contained data.

## Usage

``` r
# S4 method for class 'parquetEdgeStore'
storeRead(store, output = c("arrow", "tibble", "igraph"), minimal = TRUE, ...)

# S4 method for class 'mtxInput'
storeRead(store, ...)

# S4 method for class 'tenxH5Input'
storeRead(store, ...)

# S4 method for class 'cellbinGefInput'
storeRead(store, ...)

# S4 method for class 'binGefInput'
storeRead(store, ...)

# S4 method for class 'csvWideInput'
storeRead(store, ...)

# S4 method for class 'parquetExprStore'
storeRead(
  store,
  output = c("query", "tibble", "duckdb", "dgcmatrix"),
  max_rows = NULL,
  max_cols = NULL,
  ...
)

# S4 method for class 'unionParquetExprStore'
storeRead(
  store,
  fields = NULL,
  output = c("query", "tibble", "duckdb", "dgcmatrix"),
  callback = NULL,
  duckdb_params = list(),
  max_rows = NULL,
  max_cols = NULL,
  ...
)

# S4 method for class 'fileStore'
storeRead(store, ...)

# S4 method for class 'queryableStore'
storeRead(
  store,
  fields = NULL,
  output = c("query", "tibble", "duckdb"),
  callback = NULL,
  duckdb_params = list(),
  ...
)

# S4 method for class 'parquetStore'
storeRead(
  store,
  fields = NULL,
  output = c("query", "tibble", "duckdb", "sedona"),
  callback = NULL,
  duckdb_params = list(),
  omit_internals = TRUE,
  ...
)

# S4 method for class 'unionParquetStore'
storeRead(
  store,
  fields = NULL,
  output = c("query", "tibble", "duckdb", "sedona"),
  callback = NULL,
  duckdb_params = list(),
  omit_internals = TRUE,
  ...
)

# S4 method for class 'unionParquetGeomStore'
storeRead(
  store,
  extent = NULL,
  fields = NULL,
  output = c("query", "tibble", "terra", "sf", "duckdb", "sedona"),
  callback = NULL,
  duckdb_params = list(),
  omit_internals = TRUE,
  ...
)

# S4 method for class 'parquetGeomStore'
storeRead(
  store,
  extent = NULL,
  fields = NULL,
  output = c("query", "tibble", "terra", "sf", "duckdb", "sedona"),
  callback = NULL,
  duckdb_params = list(),
  omit_internals = TRUE,
  ...
)

# S4 method for class 'parquetGeomTileStore'
storeRead(
  store,
  extent = NULL,
  tile_idx = NULL,
  fields = NULL,
  output = c("query", "tibble", "terra", "sf", "duckdb", "sedona"),
  callback = NULL,
  duckdb_params = list(),
  omit_internals = TRUE,
  ...
)

# S4 method for class 'h5ArrayStore'
storeRead(store, ...)

# S4 method for class 'tileDBArrayStore'
storeRead(store, ...)

# S4 method for class 'bpcMatrixStore'
storeRead(store, ...)
```

## Arguments

- store:

  `dataStore` inheriting object

- output:

  `character` (default = "query"). Format to get values in:

  - "query" - produces an arrow lazy query

  - "tibble" - materialized dplyr tibble

  - "terra" - materialized `SpatVector`

  - "sf" - materialized `sf` object

  - "duckdb" - (requires duckdb and dbplyr) produces a `tbl_dbi` lazy
    query. **Note:** should not be used in a parallelized context as
    duckdb handles parallelization internally.

  - "sedona" - (requires sedonadb) produces a `sedonadb_dataframe`. All
    pending Arrow-phase ops (`@ops`, `@crop`, `@window`) are applied via
    the Arrow pipeline before handing off. Pending affine transforms are
    applied via `ST_Affine` on the `geom` column. **Note:**
    `x_index`/`y_index` are not updated after a transform — use SedonaDB
    spatial functions for geometry-based filtering.

- ...:

  additional params to pass (if any implemented)

- max_rows, max_cols:

  integer or `Inf` (optional). Set a dimension guard for
  `output = "dgcmatrix"`. A materialized sparseMatrix must have at least
  one axis narrowed to within the cap — both exceeding errors. Defaults
  via `getOption("giottodisk.dgc_max_rows", 100L)` /
  `getOption("giottodisk.dgc_max_cols", 100L)`. Asymmetric so "narrow
  slice" usage (100 × Inf or Inf × 100) is allowed; the intent is
  `[`-subset along one axis before materializing. Pass `Inf` to disable
  the cap on an axis (`Inf`/`Inf` disables the guard entirely).

- fields:

  `character` (optional) specific fields/columns to read

- callback:

  (optional) `function` where the first param should accept the arrow
  query. A function to apply to the query prior to output and after
  fields or other filters are applied.

  Mostly useful for outputs that require materialization.

- duckdb_params:

  named `list`. Params to pass to `duckdb::duckdb_register_arrow()` if
  `output = "duckdb"`. Key params:

  - `conn` - DBI connection to a duckdb instance.

  - `name` - `character` (optional) If not provided, a random ID for the
    registered table will be generated

- omit_internals:

  `logical` (default `TRUE`). Whether to drop internal special columns
  (`row_index`, `source_id`) from materialized output (`"tibble"`,
  `"terra"`, `"sf"`). Set to `FALSE` to retain them, e.g. when
  downstream code needs these columns for joining. Has no effect for
  `"query"` or `"duckdb"` outputs (those always expose all columns).

- extent:

  `SpatExtent` to filter on (optional)

- tile_idx:

  `integerlike` (optional) specific tile number(s) to read

## Index coordinates by output mode

`parquetExprStore` / `unionParquetExprStore` outputs differ on *two*
independent axes, and switching `output` changes both.

Beyond whether `@post_ops` are applied (see Details), the `row_id` /
`col_id` columns are not in the same coordinate system across modes:

- `output = "dgcmatrix"` returns values indexed by **subset position**,
  with `@feat_ids` / `@cell_ids` attached as dimnames. Rows whose ids
  fall outside the store's current `[`-subset are dropped.

- `output = "tibble"` and `output = "data.table"` return the raw
  **on-disk** `row_id` / `col_id`, without remapping and without
  unmappable-row filtering, even though `@post_ops` *are* applied. A
  caller that needs subset positions applies `.pe_remap_row()` /
  `.pe_remap_col()` itself, and is usually better off doing so after
  aggregating – remapping is a per-row
  [`match()`](https://rdrr.io/r/base/match.html) against `@cell_idx` /
  `@gene_idx`, so it is far cheaper on a per-gene aggregate than on a
  full triplet stream.

- `output = "query"` and `output = "duckdb"` never materialize, so ids
  are likewise raw on-disk values.

The practical consequence: moving from `"dgcmatrix"` to a tabular output
does not error, but `row_id` / `col_id` silently change meaning. The
tabular modes are expected to converge on the remapped coordinates in a
later release; until then, treat their id columns as on-disk.
