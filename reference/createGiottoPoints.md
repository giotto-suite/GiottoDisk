# Create a giottoPoints from a Parquet Geometry Store

Create a
[GiottoClass::giottoPoints](https://giotto-suite.github.io/GiottoClass/reference/giottoPoints-class.html)
object from a `parquetGeomBase` store. The store is kept in the
`spatVector` slot – nothing is materialized. The `unique_ID_cache` is
populated via a single Arrow `DISTINCT` query at construction time.

When `split_keyword` is provided, a named list of `giottoPoints` is
returned, each backed by a lazily filtered subset of the store.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
createGiottoPoints(
  x,
  feat_type = "rna",
  feat_ID_colname = "feat_ID",
  split_keyword = NULL,
  ...
)
```

## Arguments

- x:

  `parquetGeomBase` store (points)

- feat_type:

  `character`. Feature type label(s). When using `split_keyword`,
  provide one value per split plus one for the remainder.

- feat_ID_colname:

  `character`. Column in `x` containing feature IDs. Default
  `"feat_ID"`.

- split_keyword:

  `list` of character vectors of keywords. Each vector is
  [`grepl()`](https://rdrr.io/r/base/grep.html)-matched against feat IDs
  to define a split. Unmatched features form the first group. See
  [`GiottoClass::createGiottoPoints()`](https://giotto-suite.github.io/GiottoClass/reference/createGiottoPoints.html).

- ...:

  unused

## Value

`giottoPoints`, or named list of `giottoPoints` when `split_keyword` is
provided
