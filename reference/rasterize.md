# Rasterize Parquet Point Data

Rasterize a `parquetGeomStore` of points onto a template `SpatRaster`.
Rather than materializing the full geometry as a `SpatVector`, the
template raster's extent and pixel dimensions are used to define spatial
bins. Points are assigned to bins in Arrow and aggregated there – only
the small summary table is pulled into R.

Currently supports point geometries only.

## Usage

``` r
# S4 method for class 'parquetGeomStore,SpatRaster'
rasterize(x, y, field = NULL, fun = "count", ...)
```

## Arguments

- x:

  `parquetGeomStore` (points)

- y:

  `SpatRaster` template – defines extent, resolution, and CRS.

- field:

  `character` (optional). Column in `x` to aggregate. Not required when
  `fun = "count"`.

- fun:

  `character`. Aggregation function. One of `"count"` (default),
  `"sum"`, `"mean"`, `"min"`, `"max"`.

- ...:

  additional params passed to
  [`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)

## Value

`SpatRaster` with the same extent and CRS as `y`. Cells with no points
are `NA`.

## Examples

``` r
# create 500 random points
set.seed(42)
x <- rnorm(500, 0, 100)
y <- rnorm(500, 0, 100)
pts <- terra::vect(
    data.frame(
        x = x,
        y = y,
        value = x + y),
    geom = c("x", "y")
)

# write to parquetGeomStore
store <- storeCreate(type = "parquetgeom")
store <- storeWrite(store, pts)

# 20 x 20 template raster over the same extent
template <- terra::rast(
    nrows = 20, ncols = 20,
    xmin = -300, xmax = 300, ymin = -300, ymax = 300
)

# point density per cell
r_count <- rasterize(store, template)
plot(r_count)

# mean of `value` field per cell
r_mean <- rasterize(store, template, field = "value", fun = "mean")
plot(r_mean)
```
