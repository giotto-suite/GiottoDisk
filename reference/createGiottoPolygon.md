# Create a giottoPolygon from a Parquet Geometry Store

Create a
[GiottoClass::giottoPolygon](https://giotto-suite.github.io/GiottoClass/reference/giottoPolygon-class.html)
object from a `parquetGeomBase` store. The store is kept in the
`spatVector` slot – nothing is materialized. The `unique_ID_cache` is
populated via a single Arrow `DISTINCT` query at construction time.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
createGiottoPolygon(x, name = "cell", poly_ID_colname = "poly_ID", ...)
```

## Arguments

- x:

  `parquetGeomBase` store (polygons)

- name:

  `character`. Polygon set name. Default `"cell"`.

- poly_ID_colname:

  `character`. Column in `x` containing polygon IDs. Default
  `"poly_ID"`.

- ...:

  unused

## Value

`giottoPolygon`
