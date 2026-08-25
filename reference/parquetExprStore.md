# Create a Parquet Expression Matrix Store

Construct a
[parquetExprStore](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md)
handle around a long-format Parquet file or directory. The Parquet must
contain `row_id`, `col_id`, `value` columns and be sorted by `row_id`.

To populate the Parquet from a 10x / Xenium MatrixMarket triple, build
an
[`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md)
and pass it to
[`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md).

## Usage

``` r
parquetExprStore(
  path = .dump_tempfile(),
  cell_ids = character(0L),
  feat_ids = character(0L),
  n_cells = length(cell_ids),
  n_genes = length(feat_ids),
  scan_stats = FALSE,
  ...
)
```

## Arguments

- path:

  character. Path to a Parquet file or directory of Parquet chunks.

- cell_ids:

  character. Cell barcodes. Length must equal `n_cells`.

- feat_ids:

  character. Gene / feature IDs. Length must equal `n_genes`.

- n_cells:

  numeric. Total number of cells. Defaults to `length(cell_ids)`.

- n_genes:

  numeric. Total number of genes. Defaults to `length(feat_ids)`.

- scan_stats:

  logical. Scan the Parquet at `path` to cache its marginal nonzero
  counts on `@stats` (default `FALSE`). Only meaningful when `path`
  already holds data: the usual pattern is to construct an empty handle
  and populate it with
  [`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md),
  which caches the marginals itself. Set `TRUE` when attaching to a
  Parquet written elsewhere and you would rather pay the scan now than
  have the first consumer pay it. Leaving it `FALSE` costs correctness
  nothing — consumers that need the counts fall back to counting on
  demand — so this is purely about when the scan happens. Reachable
  through
  [`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
  which forwards `...` here.

- ...:

  additional slots passed through to `new()`.

## Value

A `parquetExprStore` S4 object.

## See also

Other store constructors:
[`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md),
[`cellbinGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput.md),
[`csvWideInput()`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput.md),
[`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md),
[`fileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore.md),
[`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md),
[`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md),
[`nnSearchInput()`](https://giotto-suite.github.io/GiottoDisk/reference/nnSearchInput.md),
[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md),
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
