# Create a Parquet Edge Store handle

Construct a
[parquetEdgeStore](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md)
handle around an existing store root, or a yet-to-be-written one (the
directory and child parquets are materialized later by
[`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)).
The node sidecar handle is auto-derived from the store root — caller
does not pass it.

Typically not called directly — produced by
[`sourceWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite.md)
dispatching on an edge input marker
([`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md),
[`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md),
[`nnSearchInput()`](https://giotto-suite.github.io/GiottoDisk/reference/nnSearchInput.md)).

## Usage

``` r
parquetEdgeStore(
  path = .dump_tempfile(),
  type = c("kNN", "sNN", "spatial"),
  directed = FALSE,
  n_edges = NA_real_,
  ...
)
```

## Arguments

- path:

  character. Path to the store root directory. May not exist yet;
  [`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  creates it. Edges + nodes parquets are written into `<path>/edges/`
  and `<path>/nodes/` respectively.

- type:

  character. Network type ("kNN" / "sNN" / "spatial").

- directed:

  logical. Storage convention. Default FALSE.

- n_edges:

  numeric. Edge count. Auto-counted lazily if NA.

- ...:

  additional slots passed to `new()`.

## Value

[parquetEdgeStore](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md)
object.

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
[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md),
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
