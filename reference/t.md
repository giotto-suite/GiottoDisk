# Transpose a parquetGeomBase store

Lazily record a transpose (swap x and y coordinates) on a
`parquetGeomBase`-inheriting store.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
t(x)
```

## Arguments

- x:

  `parquetGeomBase`-inheriting store

## Value

`x` with the transpose recorded in `@post_ops`
