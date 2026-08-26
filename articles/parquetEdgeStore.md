# Disk-backed networks with parquetEdgeStore

``` r

library(GiottoDisk)
```

`parquetEdgeStore` is the disk-backed counterpart to in-memory networks
(`nnNetObj@network`, `spatialNetworkObj@network`). It stores edges in
long-format parquet plus a small node-ID sidecar, supports lazy filter
pushdown via `[`, and materializes back to `igraph` / `data.table` /
`arrow` on demand. It’s the storage layer behind
`createNetwork(..., output = "parquet")` and the auto-write path on
network setters for backed `giotto` objects.

## A small synthetic example

We use a 20-vertex undirected random graph throughout this article so
every chunk is runnable end-to-end. Substitute your own network / giotto
pipeline where the comments say so.

``` r

set.seed(42)
ig <- igraph::sample_gnm(20, 50, directed = FALSE)
igraph::V(ig)$name <- paste0("cell_", sprintf("%02d", seq_len(20)))
igraph::E(ig)$weight <- runif(50)
ig
#> IGRAPH a74513a UNW- 20 50 -- Erdos-Renyi (gnm) graph
#> + attr: name (g/c), type (g/c), loops (g/l), m (g/n), name (v/c),
#> | weight (e/n)
#> + edges from a74513a (vertex names):
#>  [1] cell_01--cell_02 cell_03--cell_04 cell_01--cell_05 cell_03--cell_05
#>  [5] cell_02--cell_06 cell_04--cell_06 cell_06--cell_07 cell_02--cell_08
#>  [9] cell_04--cell_08 cell_07--cell_08 cell_02--cell_09 cell_03--cell_09
#> [13] cell_08--cell_09 cell_03--cell_10 cell_04--cell_10 cell_05--cell_10
#> [17] cell_05--cell_11 cell_08--cell_11 cell_01--cell_12 cell_02--cell_12
#> [21] cell_10--cell_12 cell_11--cell_12 cell_01--cell_13 cell_11--cell_13
#> [25] cell_02--cell_14 cell_06--cell_14 cell_07--cell_14 cell_11--cell_14
#> + ... omitted several edges
```

## Direct opt-in — `createNetwork(..., output = "parquet")`

The explicit way: set `output = "parquet"` on the `networkParam` and
`createNetwork()` returns a `parquetEdgeStore` directly. (We simulate
the call here with a `storeWrite` of the in-mem igraph, which is what
`.finalize_network` does under the hood.)

``` r

# In real code:
#   p <- networkParam("sNN", k = 30, output = "parquet")
#   edge_store <- createNetwork(pca, param = p)
edge_store <- storeWrite(
    storeCreate(type = "parquetEdgeStore"),
    ig, type = "sNN"
)
edge_store
#> <parquetEdgeStore> type=sNN directed=FALSE
#>   n_cells: 20  n_edges: 50
#>   path:    /tmp/Rtmp50Y0WK/gdisk_dump/file365559b92a10
```

If a `gDirSource` is supplied, the artifact lands in the project vault
instead of a tempdir. With the auto-write described below, you usually
don’t have to set this flag explicitly.

## Auto-write on backed gobjects

When a `giotto` object has a `gsource` backend attached, the network
setters (`setNearestNetwork()` / `setSpatialNetwork()`) automatically
write any in-memory igraph through to a `parquetEdgeStore` before
storing it. The user-facing call site is unchanged from the in-memory
case.

``` r

gdir <- file.path(tempdir(), "vignette_autowrite")
unlink(gdir, recursive = TRUE)

mat <- matrix(rpois(20 * 50, 2), nrow = 50, ncol = 20,
              dimnames = list(paste0("g_", 1:50),
                              paste0("cell_", sprintf("%02d", 1:20))))
g <- GiottoClass::createGiottoObject(expression = mat, backend = gdir)
#> Setting up Giotto project directory at: 
#> /tmp/Rtmp50Y0WK/vignette_autowrite
#> checking default envname 'giotto_env'
#> a system default python environment was found
#> Using python path:
#>  "/usr/bin/python3"
#> Warning: Some of Giotto's expected python module(s) were not found:
#> pandas, igraph, leidenalg, community, networkx, sklearn
#> (This is fine if python-based functions are not needed)
#> 
#> ** Python path used: "/usr/bin/python3"

nn <- methods::new("nnNetObj",
    network    = ig,           # in-memory igraph
    nn_type    = "sNN",
    name       = "sNN.demo",
    spat_unit  = "cell",
    feat_type  = "rna",
    provenance = "cell"
)
options("giotto.check_valid" = FALSE)
g <- GiottoClass::setNearestNetwork(g, nn, verbose = FALSE)

nn_back <- GiottoClass::getNearestNetwork(g, output = "nnNetObj",
    spat_unit = "cell", feat_type = "rna",
    nn_type   = "sNN", name = "sNN.demo")
class(nn_back@network)
#> [1] "parquetEdgeStore"
#> attr(,"package")
#> [1] "GiottoDisk"
nn_back@network
#> <parquetEdgeStore> type=sNN directed=FALSE
#>   n_cells: 20  n_edges: 50
#>   path:    /tmp/Rtmp50Y0WK/vignette_autowrite/artifacts/13909_vHKg5r52/data
```

The auto-write fires when both:

1.  `gobject@source` is non-null (a backend is attached), AND
2.  The incoming `@network` slot is NOT already a `dataStore` (so we
    never double-write an already-disk-backed network).

For gobjects without a backend, the setter is a no-op on storage shape:
in-memory networks stay in memory.

## Bulk promotion — `sourceWrite(src, gobject)`

If you have an existing in-memory `giotto` (e.g. from an analysis that
predates the backed-pipeline switch), the whole thing can be promoted in
one shot:

``` r

# An in-memory gobject built the normal way
g_inmem <- GiottoClass::createGiottoObject(expression = mat)
#> python already initialized in this session
#>  active environment : '/usr/bin/python3'
#>  python version : 3.12
g_inmem
#> An object of class giotto 
#> >Active spat_unit:  cell 
#> >Active feat_type:  rna 
#> dimensions    : 50, 20 (features, cells)
#> [SUBCELLULAR INFO]
#> [AGGREGATE INFO]
#> expression -----------------------
#>   [cell][rna] raw
#> spatial locations ----------------
#>   [cell] raw
#> 
#> 
#> Use objHistory() to see steps and params used
class(g_inmem@source)   # NULL — no backend
#> [1] "NULL"

# Promote to a backed project
gdir_bulk <- file.path(tempdir(), "vignette_bulk_promote")
unlink(gdir_bulk, recursive = TRUE)
g_backed <- sourceWrite(gDirSource(path = gdir_bulk), g_inmem)
#> Setting up Giotto project directory at: 
#> /tmp/Rtmp50Y0WK/vignette_bulk_promote
class(g_backed@source)
#> [1] "gDirSource"
#> attr(,"package")
#> [1] "GiottoDisk"
class(GiottoClass::getExpression(g_backed)[])
#> [1] "parquetExprStore"
#> attr(,"package")
#> [1] "GiottoDisk"
```

Internally this builds a fresh `giotto` with the gsource attached as the
backend, then re-attaches every data subobject via
\[[`GiottoClass::setGiotto()`](https://giotto-suite.github.io/GiottoClass/reference/setGiotto.html)\].
Each setter’s existing backend-aware write path handles the conversion:
expression matrices land as `parquetExprStore`, polygons as
`parquetGeomStore`, points as `parquetGeomStore`, NN and spatial
networks as `parquetEdgeStore`.

The operation is idempotent at the subobject level — every backend-
aware setter checks `inherits(x@<slot>, "dataStore")` and skips re-
writing if already on disk. Calling `sourceWrite(src, g_backed)` again
with a fully-backed gobject doesn’t grow the manifest or re-write
artifacts. Useful if some subobject snuck through in-memory on a project
that should be fully backed — re-running this call catches the stray.

Subobjects whose setters don’t yet have backend-aware write logic
(e.g. some image variants) stay in-memory; this is a known limitation
that doesn’t affect correctness.

## Reading back

Three output modes.

``` r

# Lazy Arrow query (zero RAM, streaming downstream)
ds <- storeRead(edge_store, output = "arrow")
head(dplyr::collect(ds))
#>    from_id to_id    weight
#>      <int> <int>     <num>
#> 1:       1     2 0.3467482
#> 2:       1     5 0.7846928
#> 3:       1    12 0.2405447
#> 4:       1    13 0.4793986
#> 5:       1    15 0.5816040
#> 6:       1    20 0.3330720
```

``` r

# Collected data.table with character node IDs (joined from sidecar)
head(storeRead(edge_store, output = "tibble"))
#>    from_id   to_id    weight
#>     <char>  <char>     <num>
#> 1: cell_01 cell_02 0.3467482
#> 2: cell_01 cell_05 0.7846928
#> 3: cell_01 cell_12 0.2405447
#> 4: cell_01 cell_13 0.4793986
#> 5: cell_01 cell_15 0.5816040
#> 6: cell_01 cell_20 0.3330720
```

``` r

# In-memory igraph (vertex names preserved via V(g)$name)
g_back <- storeRead(edge_store, output = "igraph")
g_back
#> IGRAPH 6f01812 UNW- 19 50 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from 6f01812 (vertex names):
#>  [1] cell_01--cell_02 cell_01--cell_05 cell_01--cell_12 cell_01--cell_13
#>  [5] cell_01--cell_15 cell_01--cell_20 cell_02--cell_06 cell_02--cell_08
#>  [9] cell_02--cell_09 cell_02--cell_12 cell_02--cell_14 cell_02--cell_15
#> [13] cell_03--cell_04 cell_03--cell_05 cell_03--cell_09 cell_03--cell_10
#> [17] cell_03--cell_17 cell_03--cell_19 cell_03--cell_20 cell_04--cell_06
#> [21] cell_04--cell_08 cell_04--cell_10 cell_04--cell_17 cell_04--cell_20
#> [25] cell_05--cell_10 cell_05--cell_11 cell_05--cell_18 cell_06--cell_07
#> [29] cell_06--cell_14 cell_07--cell_08 cell_07--cell_14 cell_07--cell_15
#> + ... omitted several edges
head(igraph::V(g_back)$name)
#> [1] "cell_01" "cell_02" "cell_03" "cell_04" "cell_05" "cell_06"
```

The `igraph` materialization ranks node ints down to a dense `1..V`
range over the subsetted edges before construction — so the in-memory
igraph only ever sees small int IDs, regardless of how large the source
universe is. `V(g)$name` carries the original character barcodes.

## Subsetting

`[` follows igraph’s adjacency semantics, with one-arg shorthand for the
common induced-subgraph case.

### Induced subgraph

``` r

v_set <- c("cell_01", "cell_02", "cell_05", "cell_10")
sub <- edge_store[v_set]
storeRead(sub, output = "tibble")
#>    from_id   to_id    weight
#>     <char>  <char>     <num>
#> 1: cell_01 cell_02 0.3467482
#> 2: cell_01 cell_05 0.7846928
#> 3: cell_05 cell_10 0.2712866
```

### Adjacency slice

``` r

storeRead(edge_store[c("cell_01"), c("cell_05", "cell_10")],
          output = "tibble")
#>    from_id   to_id    weight
#>     <char>  <char>     <num>
#> 1: cell_01 cell_05 0.7846928
```

### Endpoint filter

``` r

# Edges incident to cell_01 (undirected: matches either side)
storeRead(edge_store[, "cell_01"], output = "tibble")
#>    from_id   to_id    weight
#>     <char>  <char>     <num>
#> 1: cell_01 cell_02 0.3467482
#> 2: cell_01 cell_05 0.7846928
#> 3: cell_01 cell_12 0.2405447
#> 4: cell_01 cell_13 0.4793986
#> 5: cell_01 cell_15 0.5816040
#> 6: cell_01 cell_20 0.3330720
```

### Named-args (same as positional)

``` r

identical(
    storeRead(edge_store["cell_01", "cell_05"], output = "tibble"),
    storeRead(edge_store[from = "cell_01", to = "cell_05"],
              output = "tibble")
)
#> [1] TRUE
```

### Complement

``` r

# All edges NOT incident to cell_01
storeRead(edge_store[, "cell_01", negate = TRUE], output = "tibble")[1:5]
#>    from_id   to_id    weight
#>     <char>  <char>     <num>
#> 1: cell_02 cell_06 0.7487954
#> 2: cell_02 cell_08 0.2610880
#> 3: cell_02 cell_09 0.9828172
#> 4: cell_02 cell_12 0.0429888
#> 5: cell_02 cell_14 0.7193558
```

All `[` calls are lazy — they append a filter to the store’s `@ops` and
return a new handle. No I/O happens until
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md).

### Undirected vs directed semantics

For `@directed = FALSE` stores (sNN, spatial, undirected Delaunay),
edges are stored once in canonical form (`from_id <= to_id`). The `[`
methods make this transparent — asymmetric filters automatically OR both
orientations:

``` r

# Same set of edges regardless of orientation:
nrow(storeRead(edge_store[c("cell_01"), c("cell_05")], output = "tibble"))
#> [1] 1
nrow(storeRead(edge_store[c("cell_05"), c("cell_01")], output = "tibble"))
#> [1] 1
```

For `@directed = TRUE` stores (kNN), slices stay strict — the original
orientation is preserved.

``` r

ig_dir <- igraph::make_graph(c("a","b","b","c","c","a","d","a"),
                             directed = TRUE)
store_dir <- storeWrite(
    storeCreate(type = "parquetEdgeStore"),
    ig_dir, type = "kNN"
)
storeRead(store_dir["a", "b"], output = "tibble")
#>    from_id  to_id
#>     <char> <char>
#> 1:       a      b
storeRead(store_dir["b", "a"], output = "tibble")  # empty: no b->a
#> Empty data.table (0 rows and 2 cols): from_id,to_id
```

## Node-ID sidecar

The sidecar is a small parquet at `<store_root>/nodes/` that maps
character node IDs (cell barcodes, point IDs) to dense integer IDs used
inside the edge file. It’s auto-generated at write time.

``` r

storeRead(edge_store@nodes, output = "tibble") |> head()
#>    node_id int_id
#>     <char>  <int>
#> 1: cell_01      1
#> 2: cell_02      2
#> 3: cell_03      3
#> 4: cell_04      4
#> 5: cell_05      5
#> 6: cell_06      6
```

You don’t normally interact with it — `[` and `storeRead` handle the
char ↔︎ int translation transparently. It exists so that the in-memory
class doesn’t carry a billion-element character vector for the vertex
universe; instead the universe lives on disk and is joined when needed.

### Writer-supplied node universe

When constructing the store directly via
[`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md),
an optional `node_universe` argument preserves isolated vertices:

``` r

# Edges only mention 3 vertices; full universe has 5
small_dt <- data.table::data.table(
    from = c("v1", "v1"), to = c("v2", "v3")
)
s <- storeWrite(
    storeCreate(type = "parquetEdgeStore"),
    small_dt,
    type = "spatial",
    directed = FALSE,
    node_universe = c("v1", "v2", "v3", "v4", "v5")
)
s@n_cells   # 5 — including isolated v4, v5
#> [1] 5
storeRead(s@nodes, output = "tibble")
#>    node_id int_id
#>     <char>  <int>
#> 1:      v1      1
#> 2:      v2      2
#> 3:      v3      3
#> 4:      v4      4
#> 5:      v5      5
```

Without `node_universe`, the sidecar derives from `unique(c(from, to))`
on the edge table — vertices with no edges would not appear.

## Tile-and-loop iteration

For per-tile spatial work (Moran’s I, neighborhood enrichment), the
right shape is to load one tile’s worth of edges into an in-memory
graph, then run the inner algorithm in a loop within that tile:

``` r

# pseudocode — uses a separately built tile plan
plan <- tilework::quadtreePlan(geom_store, threshold = 50000)

results <- lapply(seq_along(plan), function(i) {
    v_set  <- storeRead(geom_store[plan[i]], output = "tibble")$cell_id
    g_tile <- storeRead(edge_store[v_set], output = "igraph")
    # ... inner per-cell loop on g_tile ...
})
```

One materialization per tile, amortized over many inner ops. The rank
step inside `storeRead(output = "igraph")` is cheap at typical chunk
sizes (~30 ms at 100k edges) and invisible at the algorithm level.

## Direct construction (escape hatches)

When the network is already in memory and you want to disk-back it
without going through `createNetwork()`:

``` r

# From an in-memory igraph (no marker class needed)
store_from_igraph <- storeWrite(
    storeCreate(type = "parquetEdgeStore"),
    ig, type = "sNN"
)
store_from_igraph
#> <parquetEdgeStore> type=sNN directed=FALSE
#>   n_cells: 20  n_edges: 50
#>   path:    /tmp/Rtmp50Y0WK/gdisk_dump/file36551589ff22
```

``` r

# From an arbitrary edge data.table with non-standard col names
my_dt <- data.table::data.table(
    src = c("a", "a", "b"),
    dst = c("b", "c", "c"),
    w   = c(0.5, 0.7, 0.9)
)
inp <- edgeDTInput(my_dt, from_col = "src", to_col = "dst", weight_col = "w")
storeWrite(storeCreate(type = "parquetEdgeStore"), inp,
           type = "kNN", directed = TRUE)
#> <parquetEdgeStore> type=kNN directed=TRUE
#>   n_cells: 3  n_edges: 3
#>   path:    /tmp/Rtmp50Y0WK/gdisk_dump/file365535e24619
```

These input markers are escape hatches — the primary path for production
work is `createNetwork(..., output = "parquet")` (or the backed-gobject
auto-write). They exist for migration from existing pipelines,
round-trip tests, and cases where the network was built outside
`createNetwork`.

## Snapshot integration

`parquetEdgeStore` participates in `gDirSource` snapshot lifecycle the
same as any other disk-backed store. With the auto-write on setters, the
typical flow is:

``` r

g <- createGiottoObject(expression = mat, backend = "~/projects/my_atlas")

# build + slot network — auto-writes through to vault
nn <- createNetwork(g, param = networkParam("sNN", k = 30))
g <- setNearestNetwork(g, nn)

# snapshot ties the artifact to the giotto save
snapshotSave(g@source, g)

# later session: artifact is recovered from the manifest, no rebuild
g <- snapshotLoad("~/projects/my_atlas")
```

The store directory contains `edges/` and `nodes/` subdirectories;
future hive-partitioned writes (e.g. `edges/source_id=<uid>/...`) drop
in without changing the read API.

## When to use `parquetEdgeStore` vs in-memory igraph

| Scenario | Recommendation |
|----|----|
| \< 10k cells, single session | in-memory igraph |
| 10k–500k cells, fits in RAM but heavy | either; parquet saves RAM during multi-step pipelines |
| 500k–10M cells | parquet — in-memory becomes a real bottleneck |
| 10M+ cells | parquet, tiled construction |
| Pipelines that round-trip to disk between R sessions | parquet — survives `snapshotSave` |

For Leiden clustering and other whole-graph compute, you can still
`storeRead(edge_store, output = "igraph")` to get an in-memory graph
when needed — the disk-backed store is the persistent canonical form,
not a replacement for in-memory graph compute.

## See also

- [`createNetwork()`](https://drieslab.github.io/GiottoClass/) on
  GiottoClass for the in-memory and parquet-output dispatch
- `parquetExprStore` — disk-backed expression matrices, same pattern
- `parquetGeomStore` — disk-backed spatial points / polygons
