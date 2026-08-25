# rbind stores

Combine stores by rbinding. Methods provided for parquetStore which
implement this behavior using union partquet dataset access. Note that
rbind for parquetStore only works on pristine stores with no lazy ops
attached. Please write out the lazy steps first via
[`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
if rbinding stores with lazy ops.

## Usage

``` r
# S4 method for class 'parquetStore,parquetStore'
rbind2(x, y, src = NULL, ...)

# S4 method for class 'parquetStore,unionParquetStore'
rbind2(x, y, src = NULL, ...)

# S4 method for class 'unionParquetStore,parquetStore'
rbind2(x, y, src = NULL, ...)

# S4 method for class 'unionParquetStore,unionParquetStore'
rbind2(x, y, src = NULL, ...)

# S4 method for class 'parquetGeomStore,parquetGeomStore'
rbind2(x, y, ...)

# S4 method for class 'unionParquetGeomStore,parquetGeomStore'
rbind2(x, y, ...)

# S4 method for class 'parquetGeomStore,unionParquetGeomStore'
rbind2(x, y, ...)

# S4 method for class 'unionParquetGeomStore,unionParquetGeomStore'
rbind2(x, y, ...)
```

## Arguments

- x, y:

  store objects to rbind

- src:

  (optional) `gsource` object to write to if rbinding the same store to
  itself. In these cases, a second copy is written to disk to make sure
  that all data is repesented on disk. Providing `src` writes to the
  `gsource`. Otherwise, it is written to the default dump location (see
  [artifact_dump](https://giotto-suite.github.io/GiottoDisk/reference/artifact_dump.md))

- ...:

  additional params to pass
