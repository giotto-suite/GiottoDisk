# Wrap an RANN / dbscan kNN result as an edge input

Expects a list with `nn.idx` (n x k integer) and `nn.dists` (n x k
float). Each input row expands to k edges at write time.

## Usage

``` r
nnSearchInput(nn, cell_ids, node_meta = NULL)
```

## Arguments

- nn:

  list with `nn.idx`, `nn.dists`.

- cell_ids:

  character. Vertex labels matching rows of `nn.idx`.

- node_meta:

  data.frame or NULL. Optional per-node sidecar.

## Value

nnSearchInput object.

## See also

Other store constructors:
[`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md),
[`cellbinGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput.md),
[`csvWideInput()`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput.md),
[`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md),
[`fileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore.md),
[`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md),
[`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md),
[`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md),
[`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md),
[`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md),
[`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md),
[`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md),
[`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
