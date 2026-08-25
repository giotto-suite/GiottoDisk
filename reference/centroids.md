# Get Polygon Centroids from a Parquet Geometry Store

Returns a modified copy of `x` flagged to read as points using the
pre-computed `x_index`/`y_index` centroid columns instead of the `geom`
WKB column. No data is read or rewritten – the flag is applied lazily at
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
time.

`@geomtype` is also updated to `"points"` so downstream dispatch is
consistent with a point store.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
centroids(x, ...)
```

## Arguments

- x:

  `parquetGeomBase` store

- ...:

  unused

## Value

modified `parquetGeomBase` store
