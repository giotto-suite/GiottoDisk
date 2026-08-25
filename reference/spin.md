# Rotate a parquetGeomBase store

Lazily record a rotation on a `parquetGeomBase`-inheriting store. The
rotation is composed with any existing pending transform and applied at
read time. The centroid anchor is live-scanned so that any preceding
crops or filters are reflected in the rotation pivot.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
spin(x, angle, x0 = NULL, y0 = NULL, ...)
```

## Arguments

- x:

  `parquetGeomBase`-inheriting store

- angle:

  rotation angle in degrees

- x0, y0:

  pivot coordinates. If `NULL` (default), the centroid of the current
  data extent is used.

- ...:

  additional arguments (ignored)

## Value

`x` with the rotation recorded in `@ops`
