# Store UIDs

Return the UID(s) for a store. For union stores, one UID is returned per
substore. For `overlapPointDisk`, returns a named list with `poly` and
`feat` UID vectors capturing the provenance of the overlap computation.

## Usage

``` r
# S4 method for class 'fileStore'
storeUID(x, ...)

# S4 method for class 'unionParquetStore'
storeUID(x, ...)

# S4 method for class 'overlapPointDisk'
storeUID(x, ...)
```

## Arguments

- x:

  `store` or `overlapPointDisk` object

- ...:

  unused

## Value

`character` vector of UIDs, or named `list` for `overlapPointDisk`
