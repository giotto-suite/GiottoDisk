# Virtual Union of Expression Stores

Lazy column-wise (cell-wise) concatenation of N
[parquetExprStore](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md)
objects. Substores must share an identical `feat_ids` vector (same
panel, same ordering); `cell_ids` accumulate across substores and must
be globally unique (caller pre-prefixes if needed).

Construction is O(1) — the union is purely virtual via Arrow's
`UnionDataset` over the substores' on-disk hive-partitioned datasets. No
data is rewritten; substores remain independently usable.

## Slots

- `stores`:

  list of
  [parquetExprStore](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md)
  objects

- `cell_ids`:

  character. Concatenated cell barcodes.

- `feat_ids`:

  character. Shared feature IDs.

- `n_cells`:

  numeric. Sum of substore n_cells.

- `n_genes`:

  numeric. Shared feature count.

- `params`:

  list. Reserved for downstream pipeline metadata.

- `ops`:

  list. The pre-materialization prefix of the chain. Mirrors
  `parquetExprStore@ops` — same record schema, same `.pe_apply_op`
  executor. Axis-keyed payloads carry one entry per substore keyed by
  `uid`, so a single union-level record covers every substore. Substores
  must have empty `@ops` at union construction time (see constructor);
  the union's own `@ops` carries any subsequent recipes.

- `post_ops`:

  list. The suffix of the chain, running from the first step that cannot
  execute in Acero onward — applied R-side to the materialized
  `data.table` after the `@ops` plan has run. A lowerable record can sit
  here legitimately: once something has forced materialization,
  everything after it must follow. Same pure-data record schema as
  `@ops`; the executor is `.pe_apply_post_ops_df()`. Axis-keyed payloads
  (e.g. a `multiply` op's factor vectors) carry one entry per substore
  keyed by `uid`. Substores must have empty `@post_ops` at construction
  (see constructor); at read time the union transplants its own `@ops`
  and `@post_ops` onto each substore via `.exprbase_inject_parent_ops()`
  so per-substore chunk reads stay self-sufficient.
  [`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  bakes the chain into on-disk values, leaving the output store with
  empty chains. Once `@post_ops` is non-empty, subsequent op pushes are
  routed here regardless of their natural phase (monotonic phase rule).
  See `R/utils-pestore-ops.R`.

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
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md)
