# Create a wide-format CSV input

Create a wide-format CSV input

## Usage

``` r
csvWideInput(
  csv_path,
  cell_id_col = "cell_ID",
  skip_cols = character(0L),
  row_filter_fun = NULL,
  batch_rows = 5000L
)
```

## Arguments

- csv_path:

  character. Path to `.csv` or `.csv.gz`.

- cell_id_col:

  character. CSV column carrying the cell barcode (default `"cell_ID"`).

- skip_cols:

  character. Additional non-feature columns to ignore.

- row_filter_fun:

  function. Optional per-chunk row filter.

- batch_rows:

  integer. Cells per batch. Default 5,000.

## Value

A `csvWideInput` object. `cell_ids` / `n_cells` are empty until
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
has been driven to completion.

## See also

Other store constructors:
[`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md),
[`cellbinGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput.md),
[`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md),
[`fileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore.md),
[`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md),
[`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md),
[`nnSearchInput()`](https://giotto-suite.github.io/GiottoDisk/reference/nnSearchInput.md),
[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md),
[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md),
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
