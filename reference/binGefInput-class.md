# Stereo-seq bin GEF Input

Wraps a Stereo-seq bin `.gef` file (`geneExp/<bin_size>/expression`
compound dataset). Unlike cellbin, the cell identity universe (one per
unique `(x, y)` coord) is not known up-front — `(x, y) -> bin_ID` is
assigned as new coords are encountered during the gene-chunk stream.
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)'s
iterator publishes the accumulated `cell_ids` / `n_cells` via its
accessors after iteration completes.

## Slots

- `bin_size`:

  character. Bin size key under `geneExp/` (e.g. `"50"`).

- `batch_genes`:

  integer. Approximate raw-gene rows per batch.

- `gene_column`:

  character.

## See also

Other store types:
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
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
