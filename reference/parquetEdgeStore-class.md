# Parquet Edge Store (streaming)

S4 class for **disk-backed networks** stored as long-format Apache
Parquet. Carries a node-ID sidecar parquet for char\<-\>int translation
so the in-memory class footprint stays bounded regardless of vertex
count.

Edge schema (locked):

|            |          |                                  |
|------------|----------|----------------------------------|
| Column     | Type     | Meaning                          |
| `from_id`  | int32/64 | node int ID; canonical sort key  |
| `to_id`    | int32/64 | node int ID                      |
| `weight`   | float32  | edge weight (e.g. sNN Jaccard)   |
| `distance` | float32  | PCA / spatial Euclidean distance |

On-disk layout:

    <store_root>/                  # @path
    ├── edges/                     # arrow::open_dataset reads here
    │   └── *.parquet              # single file today; hive-partitioned later
    └── nodes/                     # auto-derived sidecar
        └── *.parquet

Both subdirs are read via
[`arrow::open_dataset`](https://arrow.apache.org/docs/r/reference/open_dataset.html),
which handles file-or-dir transparently — extending to `source_id=<uid>`
hive partitions for rbind support is a write-time change only.

Node sidecar schema (lives at `<store_root>/nodes/`):

|           |          |          |                                    |
|-----------|----------|----------|------------------------------------|
| Column    | Type     | Required | Meaning                            |
| `node_id` | string   | yes      | original cell barcode / point ID   |
| `int_id`  | int32/64 | yes      | dense integer enumeration          |
| `*`       | any      | no       | arbitrary writer-supplied metadata |

Optional sidecar columns survive round-trips and can be joined into
results at read time. Spatial coords (`x_index`, `y_index`) are NOT
duplicated here — they live in `parquetGeomStore` and are joined across
stores when needed.

For undirected networks (sNN, spatial) edges are stored once in
canonical form (`from_id <= to_id`). kNN networks are symmetrized at
write time and stored canonical. Algorithms requiring directed traversal
can union with a swap.

## Slots

- `n_cells`:

  numeric. Vertex count (= nrow nodes sidecar).

- `n_edges`:

  numeric. Edge count.

- `nodes`:

  parquetStore handle to the node sidecar.

- `type`:

  character. One of "kNN", "sNN", "spatial".

- `directed`:

  logical. TRUE = edges stored directed as-is; FALSE = canonical
  undirected form.

## See also

[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md)

Other store types:
[`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md),
[`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md),
[`csvWideInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput-class.md),
[`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md),
[`edgeInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/edgeInput-class.md),
[`exprInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/exprInput-class.md),
[`fileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md),
[`mtxInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput-class.md),
[`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md),
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
