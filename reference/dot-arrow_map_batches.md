# arrow Utility: Map batches

Internal that applies a transformation function `FUN` to an arrow
dataset batch by batch rather than all at once, returning a new lazy
`RecordBatchReader` with the transformed data. This can be used to tell
arrow how to perform a data transformation that is not natively
supported in a memory-friendly way by performing it in batches. This
implementation is based on
[`arrow::map_batches()`](https://arrow.apache.org/docs/r/reference/map_batches.html).
This may be superceded in arrow's own API in the future. The code runs
in R instead of arrow C++ code

**Key behaviors:**

- Lazy evaluation: Processes batches on-demand (unless `.lazy = FALSE`)

- Schema inference: If schema not provided, processes first batch to
  infer output schema, then uses that schema for remaining batches

- Streaming: Returns a reader that can be consumed incrementally without
  loading entire dataset into memory

- Preserves batch structure: Each batch is transformed independently

## Usage

``` r
.arrow_map_batches(X, FUN, ..., .schema = NULL, .lazy = TRUE)
```

## Arguments

- X:

  An object coercible to
  [arrow::RecordBatchReader](https://arrow.apache.org/docs/r/reference/RecordBatchReader.html).
  Common inputs include arrow queries (`arrow_dplyr_query`),
  `FileSystemDataset`, arrow `Table`, or `data.frame`/`tibble`.

- FUN:

  `function`. A function to apply. The first parameter should accept a
  batch subset of the full dataset from `X`.

- ...:

  additional params to pass to `FUN`

- .schema:

  arrow `Schema` (optional). If not explicitly provided, the schema will
  be inferred from the first batch.

- .lazy:

  `logical`. If `TRUE` (default), returns lazy reader. If `FALSE`,
  materializes all batches immediately.

## Value

An arrow
[RecordBatchReader](https://arrow.apache.org/docs/r/reference/RecordBatchReader.html)
with transformed batches

## Examples

``` r
if (FALSE) { # \dontrun{
# Add row indices to batches
ds <- arrow::open_dataset("path/to/data")
.arrow_map_batches(ds, function(batch) {
    batch$row_id <- seq_len(nrow(batch))
    batch
})
} # }
```
