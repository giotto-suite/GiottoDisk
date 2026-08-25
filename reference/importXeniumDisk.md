# Import a 10X Xenium assay (disk-backed)

Disk-backed counterpart to
[`Giotto::importXenium()`](https://rdrr.io/pkg/Giotto/man/importXenium.html).
Produces a `XeniumDiskReader` whose `load_transcripts()` and
`load_polys()` calls write to a `gDirSource`-managed project vault as
`parquetGeomTile` stores. Other modalities (expression, featmeta,
cellmeta, images) remain in-memory via the inherited `XeniumReader`
closures.

## Usage

``` r
importXeniumDisk(xenium_dir = NULL, backend, qv_threshold = 20)
```

## Arguments

- xenium_dir:

  Xenium output directory

- backend:

  a `gsource` (typically `gDirSource`) project backend. Naming matches
  [`GiottoClass::createGiottoObject()`](https://giotto-suite.github.io/GiottoClass/reference/create_giotto.html)'s
  `backend` param.

- qv_threshold:

  minimum Phred-scaled quality score retained

## Value

`XeniumDiskReader` object

## See also

[`Giotto::importXenium()`](https://rdrr.io/pkg/Giotto/man/importXenium.html)
for the in-memory variant
