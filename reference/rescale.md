# Scale a parquetGeomBase store

Lazily record a scaling transform on a `parquetGeomBase`-inheriting
store.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
rescale(x, fx = 1, fy = fx, x0, y0, ...)
```

## Arguments

- x:

  `parquetGeomBase`-inheriting store

- fx, fy:

  scale factors for x and y axes. `fy` defaults to `fx`.

- x0, y0:

  pivot coordinates. If omitted, the centroid of the current data extent
  is used.

- ...:

  additional arguments (ignored)

## Value

`x` with the scaling recorded in `@ops`
