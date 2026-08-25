# Shear a parquetGeomBase store

Lazily record a shear transform on a `parquetGeomBase`-inheriting store.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
shear(x, fx = 0, fy = 0, x0, y0, ...)
```

## Arguments

- x:

  `parquetGeomBase`-inheriting store

- fx, fy:

  shear factors for x and y axes (default 0 = no shear)

- x0, y0:

  shear centre coordinates. If omitted, the centroid of the current data
  extent is used.

- ...:

  additional arguments (ignored)

## Value

`x` with the shear recorded in `@ops`
