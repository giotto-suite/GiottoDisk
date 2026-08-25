# Window a `parquetGeomStore`

Similar to
[`crop()`](https://giotto-suite.github.io/GiottoDisk/reference/crop.md),
but does not apply a permanent spatial subset on the data.

## Usage

``` r
# S4 method for class 'parquetGeomBase'
window(x, ...)

# S4 method for class 'parquetGeomBase'
window(x, ...) <- value
```

## Arguments

- x:

  object to window

- ...:

  additional params to pass (none implemented)

- value:

  extent to apply as a window. Can be any object that works with
  [`ext()`](https://giotto-suite.github.io/GiottoDisk/reference/ext.md)
  or `NULL` to remove the window
