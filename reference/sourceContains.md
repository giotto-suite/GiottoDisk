# Test if a Store is Managed by a Source

Returns `TRUE` if the store is already registered in the artifact vault
of a
[gsource](https://giotto-suite.github.io/GiottoDisk/reference/gsource.md).

- `fileStore`: checks whether the store's uid is present in the
  manifest.

- `SpatRaster`: checks whether all source file paths are inside the
  vault directory. In-memory rasters always return `FALSE`.

- Union stores: `TRUE` only if all substores are contained.

## Usage

``` r
# S4 method for class 'gDirSource,fileStore'
sourceContains(src, store, ...)

# S4 method for class 'gDirSource,SpatRaster'
sourceContains(src, store, ...)

# S4 method for class 'gDirSource,unionParquetStore'
sourceContains(src, store, ...)

# S4 method for class 'gDirSource,ANY'
sourceContains(src, store, ...)
```

## Arguments

- src:

  `gsource` object

- store:

  object to test (`fileStore`, `SpatRaster`, or union store)

- ...:

  additional params (unused)

## Value

`logical(1)`
