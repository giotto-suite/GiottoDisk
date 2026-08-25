# Head and tail

Queue a `head` or `tail` op on a lazy parquet store. The op is applied
at
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
time.

## Usage

``` r
# S3 method for class 'parquetBase'
head(x, n = 6L, ...)

# S3 method for class 'parquetBase'
tail(x, n = 6L, ...)
```

## Arguments

- x:

  parquetBase-inheriting store

- n:

  integer. Number of rows to keep.

- ...:

  additional arguments (ignored)

## Value

the input store with a head / tail op queued on `@ops`
