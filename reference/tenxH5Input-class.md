# 10x cell_feature_matrix.h5 Input

Wraps a 10x HDF5 sparse-matrix file (CSC layout with `data`, `indices`,
`indptr`, `barcodes`, `features/{id,name}`, `shape`). Metadata is read
eagerly at construction; the actual sparse data is streamed in
cell-chunks by
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
via hyperslab reads on a long-lived `hdf5r::H5File` handle.

## Slots

- `batch_cells`:

  integer. Cells per batch. Default 250,000.

- `feature_id_col`:

  integer. `1L` for Ensembl ID, `2L` for gene symbol (default).

## See also

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
[`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md),
[`csvWideInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput-class.md),
[`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md),
[`edgeInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/edgeInput-class.md),
[`exprInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/exprInput-class.md),
[`fileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md),
[`mtxInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput-class.md),
[`parquetEdgeStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md),
[`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md),
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
