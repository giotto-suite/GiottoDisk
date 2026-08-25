# Parquet Expression Matrix Store (streaming)

S4 class for **disk-backed expression matrices** stored as long-format
Apache Parquet, designed for streaming read access. Unlike
[parquetStore](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md)
(which stores arbitrary tabular data with a `row_index` column),
`parquetExprStore` represents a sparse cell x gene expression matrix
using three integer / float columns:

|          |         |                    |
|----------|---------|--------------------|
| Column   | Type    | Meaning            |
| `row_id` | int32   | 1-based cell index |
| `col_id` | int32   | 1-based gene index |
| `value`  | float64 | expression count   |

Parquet files are sorted by `row_id` so that Arrow predicate-pushdown
row-group skipping makes chunked reads fast on large datasets.

Cell barcodes (`cell_IDs`) and gene names (`feat_IDs`) are stored as
character slots on the S4 object – the Parquet payload itself stays
minimal. The slot vectors act as a lookup table from integer index to
character ID.

## Slots

- `path`:

  character. Local file path (single Parquet) or directory (one Parquet
  per chunk for very large datasets). Arrow's `open_dataset` handles
  both transparently.

- `uid`:

  character. Auto-generated unique ID for artifact tracking.

- `read_fun`:

  function. Preset to
  [`arrow::open_dataset()`](https://arrow.apache.org/docs/r/reference/open_dataset.html).

- `params`:

  list. Reserved for downstream pipeline metadata (e.g. HVG indices
  after `sc_hvg`). Not used for normalization recipes anymore — those
  live on `@ops`.

- `ops`:

  list. The part of the op chain that runs **before** materialization.
  `@ops` and `@post_ops` are one ordered sequence split at the point
  where execution leaves Acero — this is the prefix, not a collection of
  whichever steps happen to be lowerable.

  Each entry is a pure-data `list(type, ...params)` record (no
  closures), so the chain survives `saveRDS` cleanly. At
  [`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  time `.pe_apply_op()` translates each into arrow-dplyr steps on the
  lazy query, and the whole prefix compiles into one plan executed once
  at collect. Every record here must therefore lower to arrow — a
  consequence of the position rather than the slot's definition.

  Independent of `@cell_idx` / `@gene_idx`: `[` narrows the window and
  never touches the chain, and the chain never consults the window at
  read time. An op whose meaning depends on the window (library
  normalization, whose factors come from column sums over the features
  in view) freezes that statistic into its payload when the producing
  verb runs. Re-running the verb is how you ask for a statistic over a
  new population; subsetting is not. See `adr/0006`.

  Empty by default; populated by `processData()` methods. See
  `R/utils-pestore-ops.R` for the op type registry.

- `n_cells`:

  numeric. Number of cells in the dataset (length of `cell_ids`).

- `n_genes`:

  numeric. Number of genes / features (length of `feat_ids`).

- `cell_ids`:

  character. Cell barcodes; index `i` corresponds to `row_id == i` in
  the Parquet file.

- `feat_ids`:

  character. Gene / feature IDs; index `j` corresponds to `col_id == j`
  in the Parquet file.

- `stats`:

  list. Cached marginal counts for the Parquet payload, filled by
  [`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md).
  Two integer vectors:

  `col_nnz`

  :   stored-entry count per feature, indexed by on-disk `col_id`.

  `row_nnz`

  :   stored-entry count per cell, indexed by on-disk `row_id`.

  Keyed by **on-disk id**, not by identifier name and not by view
  position, which is what makes them invariant under `[`: subsetting
  only narrows `@cell_idx` / `@gene_idx` against the same file, so
  `sum(stats$col_nnz[gene_idx])` is the exact nonzero count of a
  gene-narrowed view. Names would break under feature renaming; view
  positions would be invalidated by every subset.

  Their lengths are the *file's* dimensions, which after a subset is the
  only place that survives — `@n_genes` / `@n_cells` describe the view.

  A subset on one axis is exact; a subset on both scales the exact axis
  by the other's kept fraction, which assumes nonzeros spread uniformly
  across the cell axis. Empty for a store whose Parquet was not written
  through
  [`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md);
  consumers fall back to counting.

  Note
  [`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  **renumbers** ids when the input is subset, so a written store's
  marginals are computed against the new file rather than inherited from
  its parent.

- `cell_idx`:

  integer. Active cell positions in the original Parquet (length 0 = no
  subset, all cells are active). When non-empty:
  `length(cell_idx) == length(cell_ids) == n_cells`. This is *view*
  state and is independent of the op chain — narrowing it never
  invalidates a queued op, because window-dependent ops froze their
  statistic at push time (`adr/0006`).

- `gene_idx`:

  integer. Active gene positions in the original Parquet (length 0 = no
  subset).

## Use case

This store is the streaming expression backend for datasets too large to
hold in memory as a sparse matrix. It is slotted into `exprObj@exprMat`
like any other disk-backed store; downstream Giotto methods that
recognize the class dispatch to streaming implementations, which read
the triplet payload in cell chunks rather than materializing the whole
matrix.

## Subset semantics

`parquetExprStore` supports lightweight subsetting via the `[` operator.
`pe[i, j]` returns a new store whose `feat_ids` / `cell_ids` are
narrowed to the kept rows / columns and whose `gene_idx` / `cell_idx`
slots record the *original* parquet positions of the kept entries. The
Parquet file on disk is **not** rewritten –
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
filters lazily via Arrow using the recorded indices, so chained subsets
stay cheap.

## See also

[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md)

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
[`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md),
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
