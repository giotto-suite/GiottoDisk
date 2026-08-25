# Create a 10x cell_feature_matrix.h5 input

Opens the h5 briefly to read `barcodes`, `features/{id,name}`, and
`shape`; closes immediately. The handle is reopened by
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
for streaming.

## Usage

``` r
tenxH5Input(h5_path, feature_id_col = 2L, batch_cells = 250000L)
```

## Arguments

- h5_path:

  character. Path to `cell_feature_matrix.h5`.

- feature_id_col:

  integer. `1L` = Ensembl ID, `2L` = gene symbol (default).

- batch_cells:

  integer. Cells per batch. Default 250,000.

## Value

A `tenxH5Input` object.

## See also

Other store constructors:
[`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md),
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
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
