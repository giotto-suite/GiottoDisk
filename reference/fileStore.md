# Create a File Store

Create a `fileStore` to represent a non-standard file on disk. A
`read_fun` function must be provided where the first param of the
function can read the `path` to load in the data in a useful format.

## Usage

``` r
fileStore(path = .dump_tempfile(), read_fun, ...)
```

## Arguments

- path:

  `character`. Disk path to file

- read_fun:

  `function`. Used to read in the data based on the path.

## See also

Other store constructors:
[`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md),
[`cellbinGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput.md),
[`csvWideInput()`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput.md),
[`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md),
[`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md),
[`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md),
[`nnSearchInput()`](https://giotto-suite.github.io/GiottoDisk/reference/nnSearchInput.md),
[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md),
[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md),
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
