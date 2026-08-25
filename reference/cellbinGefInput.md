# Create a Stereo-seq cellbin GEF input

Create a Stereo-seq cellbin GEF input

## Usage

``` r
cellbinGefInput(
  gef_path,
  gene_column = c("geneName", "geneID"),
  batch_genes = 500L
)
```

## Arguments

- gef_path:

  character. Path to a Stereo-seq cellbin `.gef` file.

- gene_column:

  character. `"geneName"` (default) or `"geneID"`.

- batch_genes:

  integer. Approximate raw-gene rows per batch. Default 500.

## Value

A `cellbinGefInput` object.

## See also

Other store constructors:
[`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md),
[`csvWideInput()`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput.md),
[`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md),
[`fileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore.md),
[`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md),
[`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md),
[`nnSearchInput()`](https://giotto-suite.github.io/GiottoDisk/reference/nnSearchInput.md),
[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md),
[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md),
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
