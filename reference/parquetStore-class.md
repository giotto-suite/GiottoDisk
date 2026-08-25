# Parquet Store

S4 Class extending
[fileStore](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md)
for indexed storage of tabular data in Apache Parquet format.
`parquetStores` provide delayed/query-based access to data rather than
loading into memory. They are not intended for in-place updates or
edits. Create with
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md).

Stores represent a query interface to on-disk data and do not contain
the actual table data as state.

`[, j]` indexing is supported. `[i, ]` row indexing is not – use
filter/crop-based access instead.

[`nrow()`](https://rdrr.io/r/base/nrow.html) and
[`dim()`](https://rdrr.io/r/base/dim.html) calls return as `numeric`
instead of `integer` to support row counts exceeding R's int32 limit
(~2.1 billion).

## Slots

- `path`:

  character. Local file path or directory, or remote URI (s3://, gs://,
  az://). For remote paths, authentication is handled via environment
  variables (AWS_ACCESS_KEY_ID, etc.) or credential files. Can point to
  a single file or a directory with optional hive-style partitioning.

- `fields`:

  character. Cached column names from the parquet dataset.

- `read_fun`:

  function. Preset to
  [`arrow::open_dataset()`](https://arrow.apache.org/docs/r/reference/open_dataset.html)
  for standard parquet access. Can be customized for edge cases (see
  [fileStore](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md)).

- `window`:

  numeric(4). (`parquetGeomStore`-inheriting) Spatial window as xmin,
  xmax, ymin, ymax.

- `geomtype`:

  `character` (`parquetGeomStore`-inheriting) Type of geometry contained
  (i.e. polygons/points)

- `tiles`:

  tileIterator. (parquetGeomTileStore only) {tilework} object defining
  which tile(s) in a tile plan this store is responsible for.

- `datatype`:

  named `list` of `character` (optional). Datatype to read a column as.
  These are applied to the arrow schema at dataset opening.

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
[`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md),
[`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)

## Examples

``` r
# creation, writing, and existence tests
ps <- parquetStore()
storeExists(ps)
ps <- storeWrite(ps, mtcars)
storeExists(ps)

# object shape
nrow(ps) # numeric
ncol(ps) # integer
dim(ps) # numeric

# different outputs
storeRead(ps)
storeRead(ps, output = "tibble")

# filtering ops
ps_1 <- subset(ps, subset = cyl > gear)
ps_1
storeRead(ps_1, output = "tibble")
head(ps_1, 2) |> storeRead(output = 'tibble')
# sample size is approximate
rowSample(ps_1, 4) |> storeRead(output = "tibble")
```
