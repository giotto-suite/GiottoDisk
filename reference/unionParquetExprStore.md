# Construct a unionParquetExprStore

Validates that all substores share an identical `feat_ids` vector,
concatenates `cell_ids`, and returns the union handle. Substore files on
disk are untouched.

Equivalent to `cbind2()` over `parquetExprStore` objects — see also the
`cbind2` methods.

## Usage

``` r
unionParquetExprStore(stores)

# S4 method for class 'parquetExprStore,parquetExprStore'
cbind2(x, y, ...)

# S4 method for class 'unionParquetExprStore,parquetExprStore'
cbind2(x, y, ...)

# S4 method for class 'parquetExprStore,unionParquetExprStore'
cbind2(x, y, ...)

# S4 method for class 'unionParquetExprStore,unionParquetExprStore'
cbind2(x, y, ...)
```

## Arguments

- stores:

  list of `parquetExprStore` objects

## Value

A
[unionParquetExprStore](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
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
[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md),
[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md)
