# Translate a parquetGeomBase store

Lazily record a translation on a `parquetGeomBase`-inheriting store.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
spatShift(x, dx = 0, dy = 0, ...)
```

## Arguments

- x:

  `parquetGeomBase`-inheriting store

- dx, dy:

  translation distances along x and y axes (default 0)

- ...:

  additional arguments (ignored)

## Value

`x` with the translation recorded in `@ops`
