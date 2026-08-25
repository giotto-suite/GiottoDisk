# Streaming per-feature statistics

[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
method for a
[Giotto::featStatsParam](https://rdrr.io/pkg/Giotto/man/analyze_param.html)
on a disk-backed expression store. One streamed pass over the triplet
stream, either over every cell or partitioned by a per-cell grouping.

The grouped form is reusable beyond QC: per-cluster mean and
percent-detected is the input to a dot plot, and the group means are a
pseudobulk matrix.

## Usage

``` r
# S4 method for class 'parquetExprBase,featStatsParam'
analyzeData(x, param, ..., groups = NULL, stats = NULL)
```

## Arguments

- x:

  a `parquetExprBase` store.

- param:

  a
  [Giotto::featStatsParam](https://rdrr.io/pkg/Giotto/man/analyze_param.html).

- ...:

  additional arguments (none used).

- groups:

  optional vector of group assignments, one per cell of the current
  view, `NA` to exclude a cell. When supplied, the statistics are taken
  per (feature, group) instead of over every cell, and the result gains
  `group` and `n_cells` columns.

- stats:

  optional character vector of accumulators to compute, any of `"sum"`,
  `"sumsq"`, `"nnz"`, `"sum_det"`. Grouped path only. Emitted columns
  are whichever the requested accumulators support, so asking for less
  genuinely scans for less. Defaults to all four.

## Value

A `data.table`, one row per feature, or per (feature, group) when
`groups` is supplied.
