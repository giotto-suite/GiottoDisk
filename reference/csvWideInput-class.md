# Wide-format CSV Input

Wraps a wide-format CSV (one row per cell, one column per feature, plus
a cell-ID column). The header is read at construction to establish
`feat_ids` and `n_genes`; cells are streamed row-chunk-wise via
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
and `cell_ids` accumulate on the iterator as the file is consumed.

## Slots

- `cell_id_col`:

  character. CSV column carrying the cell barcode.

- `skip_cols`:

  character. Additional non-feature columns to ignore.

- `row_filter_fun`:

  function. Optional per-chunk row filter (e.g. CosMx `cell_ID != 0`).

- `batch_rows`:

  integer. Cells per batch. Default 5,000.

## See also

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
[`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md),
[`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md),
[`edgeInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/edgeInput-class.md),
[`exprInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/exprInput-class.md),
[`fileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md),
[`mtxInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput-class.md),
[`parquetEdgeStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md),
[`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md),
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
