# Subsample rows of a parquetStore

Record a row subsampling operation on a store. The operation is lazy and
applied when reading via
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md).
If the store already has fewer rows than `size`, the store is returned
unchanged.

Sampling is performed by systematic (evenly-spaced) selection rather
than random sampling. Given a target of `size` rows from `n` total rows,
every `k`th row is selected where `k = ceiling(n / size)`. This is
deterministic and reproducible across sessions without a seed.

The actual number of rows returned is approximate since the rows to be
kept are determined based the value of the internal `row_index` col.
Depending on any preceding
[`subset()`](https://giotto-suite.github.io/GiottoDisk/reference/subset.md)
operations, the number of rows surviving the sample filter may be
slightly above or below `size`.

## Usage

``` r
# S4 method for class 'parquetBase'
rowSample(x, size, ...)
```

## Arguments

- x:

  `parquetBase`-inheriting object

- size:

  `numeric(1)`. Maximum number of rows to return

- ...:

  additional params (none implemented)

## Value

`x` with sampling op appended to `@ops`

## See also

[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md),
[`subset()`](https://giotto-suite.github.io/GiottoDisk/reference/subset.md)

## Examples

``` r
ps <- parquetStore()
ps <- storeWrite(ps, mtcars)
ps |> rowSample(10)         # lazy, no disk read
ps |> rowSample(10) |> storeRead(output = "tibble")
```
