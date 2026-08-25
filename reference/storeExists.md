# Store Existence

Test if a store has been written and is ready to be used.

## Usage

``` r
# S4 method for class 'fileStore'
storeExists(x, ...)

# S4 method for class 'unionParquetStore'
storeExists(x, all = TRUE, ...)
```

## Arguments

- x:

  `store` object

- ...:

  additional params to pass

- all:

  `logical` return single [`all()`](https://rdrr.io/r/base/all.html)
  value instead of a vector of `logical` for each substore

## Value

`logical`
