# Visualize a Store

Plot a store's contents. Currently only for geometry stores.

## Usage

``` r
# S4 method for class 'parquetGeomStore,missing'
plot(x, sample_max = getOption("giottodisk.plot_sample_max", 1e+05), ...)

# S4 method for class 'sedonadb_dataframe,missing'
plot(x, values = NULL, n = getOption("giottodisk.plot_sample_max", 1e+05), ...)
```

## Arguments

- x:

  A `sedonadb_dataframe` from `storeRead(output = "sedona")`.

- sample_max:

  `integer`-like. Maximum number of geometries to plot (regularly
  sampled). If `NULL`, no sampling is performed. Defaults to
  `getOption("giottodisk.plot_sample_max")`.

- values:

  `character` (optional). Column name(s) to include as SpatVector
  attributes for coloring. Mirrors terra's `plot(sv, y)`.

- n:

  `integer`. Target number of rows to display (default from
  `getOption("giottodisk.plot_sample_max")`). Implemented as systematic
  stride sampling (`row_index %% k == 0`) after a `COUNT(*)` pass —
  proportional coverage across write order.
