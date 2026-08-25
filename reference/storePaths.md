# Store Paths

Return the artifact-level path(s) for a store. For composite stores
(e.g. union stores), one path is returned per substore.

## Usage

``` r
# S4 method for class 'fileStore'
storePaths(x, ...)

# S4 method for class 'unionParquetStore'
storePaths(x, ...)
```

## Arguments

- x:

  `store` object

- ...:

  additional params to pass

## Value

`character` vector of paths
