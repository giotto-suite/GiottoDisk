# Special Store Columns

Documents columns that are special to particular store implementations.
These columns are enforced on disk and available when accessing the data
as a query. These columns may be automatically consumed depending on the
materialization format.

## Usage

``` r
# S4 method for class 'ANY'
specialCols(store)

# S4 method for class 'parquetBase'
specialCols(store)

# S4 method for class 'parquetGeomBase'
specialCols(store)

# S4 method for class 'unionParquetStore'
specialCols(store)
```

## Arguments

- store:

  `dataStore`-inheriting class

## Value

`character`
