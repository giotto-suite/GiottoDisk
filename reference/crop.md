# Crop a `parquetGeomStore`

Apply a spatial crop on geometry data. The crop is more accurately an xy
selection on centroids. Polygon geometries are not modified.

Crop operations are composable operations that limit where future crops
or
[`window()`](https://giotto-suite.github.io/GiottoDisk/reference/window.md)
can be placed.

## Usage

``` r
# S4 method for class 'parquetGeomBase,ANY'
crop(x, y, ...)
```

## Arguments

- x:

  object to crop

- y:

  spatial extent to crop to. Accepts any object that works with
  [`ext()`](https://giotto-suite.github.io/GiottoDisk/reference/ext.md)

- ...:

  additional params to pass (none implemented)

## Value

`parquetGeomBase`-inheriting object

## See also

[`window()`](https://giotto-suite.github.io/GiottoDisk/reference/window.md)
