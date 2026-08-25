# 10x / Xenium MatrixMarket Triple Input

Wraps a `matrix.mtx[.gz]` + `barcodes.tsv[.gz]` + `features.tsv[.gz]`
triple. Sidecars are read at construction to populate `cell_ids` /
`feat_ids`; the mtx itself is opened lazily by
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
and streamed in batches of `batch_lines` triplets via a single
long-lived gzip / file connection.

## Slots

- `batch_lines`:

  integer. Triplets per batch. Default 5,000,000 (~120 MB peak per
  batch).

- `feature_id_col`:

  integer. Column of `features.tsv` to use as the gene identifier. `1L`
  for Ensembl ID, `2L` for gene symbol (default).

## Input layout

`@path` is the matrix.mtx.gz file. The two sidecar paths live in
`@params$barcodes_path` and `@params$features_path` (set by the
constructor; can be overridden manually).

## See also

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
[`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md),
[`csvWideInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput-class.md),
[`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md),
[`edgeInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/edgeInput-class.md),
[`exprInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/exprInput-class.md),
[`fileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md),
[`parquetEdgeStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md),
[`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md),
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
