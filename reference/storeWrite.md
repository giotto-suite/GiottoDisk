# Write to a `dataStore`

Direct data.table dispatch: the in-memory caller (e.g. GiottoClass's
`createNetwork()` machinery) hands a pre-canonicalized edge data.table
with character `from` / `to` columns plus optional `weight` / `distance`
/ `shared`. Trusts the caller for canonical form — no swap, no dedup.

igraph dispatch: extracts edges via
[`igraph::as_data_frame()`](https://r.igraph.org/reference/graph_from_data_frame.html)
and the vertex universe via `V(g)$name`, then delegates to the
data.table path. Directedness defaults to `igraph::is_directed(data)` if
not explicitly passed. Used by GiottoClass's network setters when a
giotto object has a gsource backend.

Write data to a `dataStore` inheriting class (usually a subclass of
[fileStore](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md)).
This documentation covers the basic expected API. For more functional
examples, see documentation for specific store implementations.

- [storeWrite-parquetStore](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetStore.md)

- [storeWrite-parquetGeomStore](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomStore.md)

- [storeWrite-parquetGeomTileStore](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomTileStore.md)

## Usage

``` r
# S4 method for class 'parquetEdgeStore,edgeInput'
storeWrite(
  store,
  data,
  type = c("kNN", "sNN", "spatial"),
  directed = FALSE,
  node_universe = NULL,
  ...
)

# S4 method for class 'parquetEdgeStore,data.table'
storeWrite(
  store,
  data,
  type = c("kNN", "sNN", "spatial"),
  directed = FALSE,
  node_universe = NULL,
  node_meta = NULL,
  ...
)

# S4 method for class 'parquetEdgeStore,igraph'
storeWrite(
  store,
  data,
  type = c("kNN", "sNN", "spatial"),
  directed = NULL,
  node_universe = NULL,
  node_meta = NULL,
  ...
)

# S4 method for class 'parquetExprStore,exprInput'
storeWrite(store, data, ...)

# S4 method for class 'parquetExprStore,parquetExprStore'
storeWrite(store, data, ...)

# S4 method for class 'parquetExprStore,unionParquetExprStore'
storeWrite(store, data, ...)

# S4 method for class 'parquetExprStore,memoryMatrix'
storeWrite(store, data, ...)

# S4 method for class 'fileStore,fileStore'
storeWrite(store, data, ...)

# S4 method for class 'h5ArrayStore,memoryMatrix'
storeWrite(store, data, ...)

# S4 method for class 'tileDBArrayStore,memoryMatrix'
storeWrite(store, data, ...)

# S4 method for class 'bpcMatrixStore,memoryMatrix'
storeWrite(store, data, ...)

# S4 method for class 'bpcMatrixStore,ANY'
storeWrite(store, data, ...)

# S4 method for class 'bpcMatrixStore,mtxInput'
storeWrite(store, data, ...)

# S4 method for class 'bpcMatrixStore,tenxH5Input'
storeWrite(store, data, ...)
```

## Arguments

- store:

  `dataStore` inheriting class

- data:

  data to write

- ...:

  additional params to pass

## See also

Other storeWrite methods:
[`storeWrite-parquetGeomStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomStore.md),
[`storeWrite-parquetGeomTileStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomTileStore.md),
[`storeWrite-parquetStore`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetStore.md)
