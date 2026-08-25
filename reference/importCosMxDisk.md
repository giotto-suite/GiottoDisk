# Import a NanoString CosMx assay (disk-backed)

Disk-backed counterpart to
[`Giotto::importCosMx()`](https://rdrr.io/pkg/Giotto/man/importCosMx.html).
Produces a `CosMxDiskReader` whose `load_expression()` call writes a
`parquetExprStore` into a `gDirSource`-managed project vault.
Transcripts / polys / images / cellmeta remain in-memory via the
inherited `CosmxReader` closures.

## Usage

``` r
importCosMxDisk(
  cosmx_dir = NULL,
  backend,
  slide = 1,
  fovs = NULL,
  version = "default",
  micron = FALSE,
  px2um = 0.12028,
  poly_pref = c("mask", "csv")
)
```

## Arguments

- cosmx_dir:

  CosMx output directory

- backend:

  a `gsource` (typically `gDirSource`) project backend. Naming matches
  [`GiottoClass::createGiottoObject()`](https://giotto-suite.github.io/GiottoClass/reference/create_giotto.html)'s
  `backend` param.

- slide, fovs, version, micron, px2um, poly_pref:

  passed through to the parent `CosmxReader` initializer.

## Value

`CosMxDiskReader` object

## See also

[`Giotto::importCosMx()`](https://rdrr.io/pkg/Giotto/man/importCosMx.html)
for the in-memory variant
