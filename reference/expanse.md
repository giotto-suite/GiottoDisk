# Get the area of individual polygons

Compute the area covered by polygons in a `parquetGeomStore` or
`parquetGeomTileStore`. The tiled terra path parallelizes across tiles
via `tileApply`; the sedona path issues a single `ST_Area` query
regardless of tiling.

## Usage

``` r
# S4 method for class 'parquetGeomStore'
expanse(
  x,
  output = c("data.table", "named", "vector"),
  engine = c("terra", "sedona"),
  poly_id_col = "poly_ID",
  ...
)

# S4 method for class 'parquetGeomTileStore'
expanse(
  x,
  output = c("data.table", "named", "vector"),
  engine = c("terra", "sedona"),
  poly_id_col = "poly_ID",
  ...
)
```

## Arguments

- x:

  `parquetGeomStore` or `parquetGeomTileStore`

- output:

  one of `"data.table"` (default), `"named"`, or `"vector"`.
  `"data.table"` returns a `data.table` with columns `cell_ID` and
  `area`. `"named"` returns a named numeric vector with `cell_ID` as
  names. `"vector"` returns a plain unnamed numeric vector — polygon
  identity is not preserved.

- engine:

  one of `"terra"` (default) or `"sedona"`. `"sedona"` issues a single
  `ST_Area(geom)` SQL query via SedonaDB/DataFusion without
  materializing geometry in R. Pending affine transforms are applied via
  `ST_Affine` before area computation, so results are correct on lazily
  scaled or transformed stores. Requires the `sedonadb` package.

- poly_id_col:

  name of the polygon ID column (default `"poly_ID"`)

- ...:

  Arguments passed on to
  [`terra::expanse`](https://rspatial.github.io/terra/reference/expanse.html)

  `unit`

  :   character. Output unit of area. One of "m", "km", or "ha"

  `transform`

  :   logical. If `TRUE`, planar CRS are transformed to lon/lat for
      accuracy

  `byValue`

  :   logical. If `TRUE`, the area for each unique cell value is
      returned

  `zones`

  :   NULL or SpatRaster with the same geometry identifying zones in `x`

  `wide`

  :   logical. Should the results be in "wide" rather than "long"
      format?

  `usenames`

  :   logical. If `TRUE` layers are identified by their names instead of
      their numbers

## Value

depends on `output`: a `data.table`, named `numeric`, or `numeric`

## Examples

``` r
# create simple polygons
polys <- terra::vect(
    c("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
      "POLYGON ((1 0, 3 0, 3 2, 1 2, 1 0))",
      "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0))")
)
terra::values(polys) <- data.frame(poly_ID = c("a", "b", "c"))

store <- storeCreate(type = "parquetgeom")
store <- storeWrite(store, polys)

# data.table output (default)
expanse(store)

# named vector
expanse(store, output = "named")

# sedona engine (requires sedonadb)
if (FALSE) { # \dontrun{
if (requireNamespace("sedonadb", quietly = TRUE)) {
    expanse(store, engine = "sedona")
}
} # }
```
