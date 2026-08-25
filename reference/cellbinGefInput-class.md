# Stereo-seq cellbin GEF Input

Wraps a Stereo-seq cellbin `.gef` file (HDF5 compound datasets under
`cellBin/`). The `cell` and `gene` tables are read in full at
construction (small); the compound `geneExp` is streamed gene-chunk-
wise via
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
using rhdf5 hyperslab reads, respecting the safe-boundary chunk plan
that keeps duplicate-named genes together.

## Slots

- `batch_genes`:

  integer. Approximate raw-gene rows per batch (default 500). Actual
  boundaries may be expanded to keep duplicate- named gene groups
  intact.

- `gene_column`:

  character. `"geneName"` (default) or `"geneID"`.

## See also

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
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
