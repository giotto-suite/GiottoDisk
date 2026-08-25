# Base class for edge-table inputs

Virtual parent of edgeDTInput, igraphInput, nnSearchInput.

## Slots

- `from_col,to_col`:

  character. Endpoint columns in the source.

- `weight_col,distance_col`:

  character or NULL. Optional attr cols.

- `node_meta`:

  data.frame or NULL. Optional per-node metadata joined into the sidecar
  at write time. Must have a `node_id` column matching the input's
  vertex IDs.

## See also

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
[`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md),
[`csvWideInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput-class.md),
[`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md),
[`exprInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/exprInput-class.md),
[`fileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md),
[`mtxInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput-class.md),
[`parquetEdgeStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md),
[`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md),
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
