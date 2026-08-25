# Lazily record an affine transform on a parquetGeomBase store

Record an affine transformation on a `parquetGeomBase`-inheriting store.
The transform is stored as a `"transform"` entry in `@post_ops` and
applied to output coordinates at materialization time
(`storeRead(output = "tibble"/"terra"/"sf")`). It is a no-op for the
Arrow query phase.

Subsequent `affine()` calls compose with any existing pending transform
– the chain collapses to a single `affine2d` entry. Spatial filters
([`crop()`](https://giotto-suite.github.io/GiottoDisk/reference/crop.md)/`window<-`)
applied after an `affine()` call back-project their query extents to
intrinsic (on-disk) space for accurate Arrow filtering.

`output = "query"` and `output = "duckdb"` are unaffected.

## Usage

``` r
# S4 method for class 'parquetGeomBase,affine2d'
affine(x, y, ...)
```

## Arguments

- x:

  `parquetGeomBase`-inheriting store

- y:

  `affine2d` transform object (from GiottoClass)

- ...:

  additional arguments (ignored)

## Value

`x` with the transform recorded in `@post_ops`
