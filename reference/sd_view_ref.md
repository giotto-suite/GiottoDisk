# Get the registered view name for a GiottoDisk SedonaDB dataframe

Returns the DataFusion view name as a double-quoted SQL identifier,
suitable for embedding directly in `sedonadb::sd_sql()` queries. Only
valid on `sedonadb_dataframe` objects returned by
`storeRead(output = "sedona")`.

## Usage

``` r
sd_view_ref(sdf)
```

## Arguments

- sdf:

  A `sedonadb_dataframe` from `storeRead(output = "sedona")`.

## Value

A length-1 character: the double-quoted view name, e.g. `'"gd_abc123"'`.
