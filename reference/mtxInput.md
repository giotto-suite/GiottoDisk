# Create a 10x / Xenium MatrixMarket triple input

Eagerly reads `barcodes.tsv` and `features.tsv` so `cell_ids`,
`feat_ids`, `n_cells`, `n_genes` are known up-front. The mtx file itself
is not opened until
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md).

Three input modes are supported:

1.  **Directory** (10x / Xenium `cell_feature_matrix/` layout):
    `mtxInput(dir)` — auto-resolves `matrix.mtx[.gz]`,
    `barcodes.tsv[.gz]`, `features.tsv[.gz]` inside the directory.

2.  **Matrix file with sidecars resolvable in its parent**:
    `mtxInput(mtx_path)` — used when the three files live as siblings.

3.  **Bare `.mtx` with explicit IDs**:
    `mtxInput(mtx_path, cell_ids = ..., feat_ids = ...)` — for non-10x
    MatrixMarket files where the barcodes/features sidecars aren't
    present. `cell_ids` and `feat_ids` must both be supplied.

Explicit `cell_ids` / `feat_ids` always override sidecar resolution.

## Usage

``` r
mtxInput(
  mtx_path,
  barcodes_path = NULL,
  features_path = NULL,
  cell_ids = NULL,
  feat_ids = NULL,
  feature_id_col = 2L,
  batch_lines = 5000000L
)
```

## Arguments

- mtx_path:

  character. Path to a directory (10x layout) or a `matrix.mtx[.gz]`
  file.

- barcodes_path:

  character. Path to `barcodes.tsv[.gz]`. Default: auto-resolved
  relative to `mtx_path`. Ignored if `cell_ids` is supplied.

- features_path:

  character. Path to `features.tsv[.gz]`. Default: auto-resolved
  relative to `mtx_path`. Ignored if `feat_ids` is supplied.

- cell_ids:

  character. Explicit cell barcodes. If supplied, `barcodes_path` is not
  consulted.

- feat_ids:

  character. Explicit feature IDs. If supplied, `features_path` and
  `feature_id_col` are not consulted.

- feature_id_col:

  integer. `1L` = Ensembl ID, `2L` = gene symbol (default). Used only
  when reading from `features_path`.

- batch_lines:

  integer. Triplets per batch. Default 5,000,000.

## Value

An `mtxInput` object.

## See also

Other store constructors:
[`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md),
[`cellbinGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput.md),
[`csvWideInput()`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput.md),
[`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md),
[`fileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore.md),
[`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md),
[`nnSearchInput()`](https://giotto-suite.github.io/GiottoDisk/reference/nnSearchInput.md),
[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md),
[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md),
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
