# data.table Utility: Add row index

Internal that adds a sequential row index column to a data.frame using
data.table's by-reference modification. This modifies the input object
in-place for memory efficiency with large in-memory datasets.

This is the in-memory counterpart to
[`.arrow_add_row_index()`](https://giotto-suite.github.io/GiottoDisk/reference/dot-arrow_add_row_index.md)
for use with materialized data.frames rather than lazy arrow queries.

## Usage

``` r
.dt_set_row_index(x, offset = 1L, col = "row_index")
```

## Arguments

- x:

  `data.frame` or `data.table`. Will be converted to data.table if not
  already.

- offset:

  `integer`. Value to add to the row index. Default is `1`, meaning the
  first row will be assigned index `1`. Use higher offsets when
  processing data chunks that continue from previous batches.

- col:

  `character`. Name of the column to create for the row index. Default
  is `"row_index"`.

## Value

The modified data.table with the new index column added. Note that the
input `x` is also modified by reference.

## Examples

``` r
if (FALSE) { # \dontrun{
# Add row_index starting from 1
df <- data.frame(x = 1:100, y = letters[1:100])
.dt_set_row_index(df)

# Add row_index starting from 1000 (for tile/chunk continuation)
.dt_set_row_index(df, offset = 1000)

# Custom column name
.dt_set_row_index(df, col = "id")
} # }
```
