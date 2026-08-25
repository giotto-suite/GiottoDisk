# Write to a Parquet Storage Spec

Write tabular data to a parquet. Enforces the presence of an integer
`row_index` column since parquet does not intrinsically encode row
order.

The `ANY` signature method is provided for writing lazily accessed data
in batches with possible pre-processing via `callback`.

## Usage

``` r
# S4 method for class 'parquetStore,data.frame'
storeWrite(
  store,
  data,
  callback = NULL,
  row_offset = 0L,
  uid_partition = TRUE,
  tile_idx = NULL,
  .arrow_meta = NULL,
  ...
)

# S4 method for class 'parquetStore,ANY'
storeWrite(
  store,
  data,
  callback = NULL,
  row_offset = 0L,
  uid_partition = TRUE,
  ...
)
```

## Arguments

- store:

  `dataStore` inheriting class

- data:

  data to write

- callback:

  `function` (optional). Function to apply to `data` before writing. The
  first param of this function should accept the `data`.

- row_offset:

  `numeric`, will be coerced to `integer` (default = 0L). Offset to
  apply to row number indexing. For example `row_offset = 0L` means the
  first `row_index` value in this file will be `1`.

- uid_partition:

  `logical` When `TRUE`, tags the store uid as a column called
  `"source_id"` using filepath hive partitioning rules.

- .arrow_meta:

  `named list` of `character` strings to add as parquet metadata

- ...:

  additional params to pass to
  [`arrow::write_dataset()`](https://arrow.apache.org/docs/r/reference/write_dataset.html)

## Value

A `parquetStore` inheriting class

## See also

Other storeWrite methods:
[`storeWrite,parquetEdgeStore,edgeInput-method`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md),
[`storeWrite-parquetGeomStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomStore.md),
[`storeWrite-parquetGeomTileStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomTileStore.md)
