# Expression Matrix Input (virtual)

Virtual base for raw expression-matrix sources. Subclasses describe a
format-specific on-disk layout and expose a batch iterator via
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md);
[`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
on a
[parquetExprStore](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md)
consumes that iterator to materialise a sorted long-format parquet.

Metadata slots (`cell_ids`, `feat_ids`, `n_cells`, `n_genes`) are
populated eagerly at construction when the format allows (e.g. mtx reads
its sidecars). For formats where cell identity is only known during the
stream (e.g. wide CSV), the iterator returned by
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
is responsible for mutating these as it advances.

## Slots

- `cell_ids`:

  character. Cell barcodes (length `n_cells`).

- `feat_ids`:

  character. Gene / feature IDs (length `n_genes`).

- `n_cells`:

  integer. Total cells.

- `n_genes`:

  integer. Total features.

## See also

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
[`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md),
[`csvWideInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput-class.md),
[`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md),
[`edgeInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/edgeInput-class.md),
[`fileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md),
[`mtxInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput-class.md),
[`parquetEdgeStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md),
[`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md),
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
