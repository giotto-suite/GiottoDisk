# Flip a parquetGeomBase store

Lazily record a reflection on a `parquetGeomBase`-inheriting store.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
flip(x, direction = "vertical", x0 = 0, y0 = 0, ...)
```

## Arguments

- x:

  `parquetGeomBase`-inheriting store

- direction:

  `"vertical"` (default, reflects over a horizontal axis) or
  `"horizontal"` (reflects over a vertical axis)

- x0, y0:

  coordinate of the reflection axis (default 0). For a vertical flip,
  `y0` is the y-position of the axis; for horizontal, `x0`.

- ...:

  additional arguments (ignored)

## Value

`x` with the reflection recorded in `@ops`
