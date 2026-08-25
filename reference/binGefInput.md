# Create a Stereo-seq bin GEF input

Create a Stereo-seq bin GEF input

## Usage

``` r
binGefInput(
  gef_path,
  bin_size,
  gene_column = c("geneName", "geneID"),
  batch_genes = 500L
)
```

## Arguments

- gef_path:

  character.

- bin_size:

  character or integer. Bin size key under `geneExp/` (e.g. `50`,
  `"50"`, `"100"`).

- gene_column:

  character. `"geneName"` (default) or `"geneID"`.

- batch_genes:

  integer. Default 500.

## Value

A `binGefInput` object. `cell_ids` / `n_cells` are empty until
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
has been driven to completion.

## See also

Other store constructors:
[`cellbinGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput.md),
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
