# arrow Utility: Add row index

Internal that adds a sequential row index column to arrow data processed
batch-by-batch. The index counter is maintained across batches using a
closure, ensuring continuous sequential numbering even when data is
processed in chunks.

This is particularly useful when writing large datasets to parquet where
row indices need to be assigned without loading the entire dataset into
memory.

## Usage

``` r
.arrow_add_row_index(data, col = "row_index", offset = 0L)
```

## Arguments

- data:

  An object coercible to
  [arrow::RecordBatchReader](https://arrow.apache.org/docs/r/reference/RecordBatchReader.html).
  Common inputs include arrow queries (`arrow_dplyr_query`),
  `FileSystemDataset`, arrow `Table`, or `data.frame`/`tibble`.

- col:

  `character`. Name of the column to create for the row index. Default
  is `"row_index"`.

- offset:

  `integer`. Starting value for the row index. Default is `0`, meaning
  the first row will be assigned index `1`. Use non-zero offsets when
  appending to existing indexed data or processing tiles.

## Value

An arrow
[RecordBatchReader](https://arrow.apache.org/docs/r/reference/RecordBatchReader.html)
with the new index column added. The indices continue sequentially
across all batches.

## Examples

``` r
if (FALSE) { # \dontrun{
# Add row_index starting from 1
ds <- arrow::open_dataset("path/to/data")
indexed_ds <- .arrow_add_row_index(ds)

# Add row_index starting from 1000 (for appending)
indexed_ds <- .arrow_add_row_index(ds, offset = 1000)

# Custom column name
indexed_ds <- .arrow_add_row_index(ds, col = "id")
} # }
```
