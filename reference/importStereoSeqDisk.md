# Import a BGI Stereo-seq assay (disk-backed)

Disk-backed counterpart to
[`Giotto::importStereoSeq()`](https://rdrr.io/pkg/Giotto/man/importStereoSeq.html).
Produces a `StereoSeqDiskReader` whose `load_expression()` call streams
the source `.gef` file into a `parquetExprStore` written to the
`gDirSource`-managed project vault. For `type = "bin"`, spatial
locations are derived from that same streaming pass. The remaining
modalities (images, masks, binpoints, polygons, and cellbin spatial
locations) come from the inherited `StereoSeqReader` closures.

## Usage

``` r
importStereoSeqDisk(
  stereoseq_dir = NULL,
  backend,
  type = c("bin", "cell"),
  bin_size = "bin100",
  gene_column = c("geneName", "geneID"),
  negative_y = TRUE,
  gef_type
)
```

## Arguments

- stereoseq_dir:

  Stereo-seq output directory

- backend:

  a `gsource` (typically `gDirSource`) project backend. Naming matches
  [`GiottoClass::createGiottoObject()`](https://giotto-suite.github.io/GiottoClass/reference/create_giotto.html)'s
  `backend` param.

- type, bin_size, gene_column, negative_y, gef_type:

  passed through to the parent `StereoSeqReader` initializer.

## Value

`StereoSeqDiskReader` object

## See also

[`Giotto::importStereoSeq()`](https://rdrr.io/pkg/Giotto/man/importStereoSeq.html)
for the in-memory variant
