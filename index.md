# GiottoDisk

On-disk processing framework for large-scale spatial genomics data.

## Overview

{GiottoDisk} provides tools for on-disk data processing and ETL. It
leverages {arrow} and {dplyr} for their powerful out-of-memory and lazy
operations, and extends them with escape hatches for chunked in-memory
operations when lazy methods are not available (e.g., geometry
construction, custom transformations).

This is powered through a unified storage abstraction layer for file
schema, disk reading and writing, enabling clean and maintainable ETL
from heterogeneous file formats characteristic of the spatial omics
field to lazy queries and/or file-backed representations.

Additionally, {GiottoDisk} provides a project management system
(`gDirSource`) for tracking on-disk artifacts, managing snapshots, and
pruning unused data.

## Installation

GiottoDisk currently builds against **development branches** of the
Giotto suite, not the release branches. `Remotes:` in `DESCRIPTION`
records them, so installing with a tool that honours it pulls the right
ones automatically:

``` r

remotes::install_github("drieslab/GiottoDisk")   # honours Remotes
# or
pak::pkg_install("drieslab/GiottoDisk")
```

To install the upstream packages explicitly:

``` r

remotes::install_github("drieslab/tilework")
remotes::install_github("giotto-suite/GiottoUtils@dev")
remotes::install_github("giotto-suite/GiottoClass@gsource")
remotes::install_github("giotto-suite/Giotto@gsource")
```

`GiottoVisuals` and `GiottoData` are not GiottoDisk dependencies — any
branch that satisfies Giotto works. Note that Giotto’s general
integration branch is `suite_dev`, not `dev`, and it does **not** carry
what GiottoDisk needs.

Why each pin exists, and what has to land upstream before it can be
dropped, is in *Upstream branch pins* in `AGENTS.md`.

## Quick Start

### Basic Usage

``` r

library(GiottoDisk)

# Create a store from spatial data
store <- parquetGeomStore(path = "path/to/storage")

# Write spatial data (SpatVector, data.frame with coordinates)
storeWrite(store, spatial_data, type = "polygon")

# Read back - returns arrow query (lazy)
data <- storeRead(store, extent = my_roi, output = "query")

# Materialize only what you need
polygons <- storeRead(store, extent = my_roi, output = "terra")
```

### Project Management

``` r

# Set up a managed project directory
src <- gDirSource("path/to/project")

# Write individual data as tracked artifacts
sourceWrite(src, my_matrix)        # matrices
sourceWrite(src, my_spatvector)    # spatial vectors
sourceWrite(src, my_dataframe)     # tables

# Save a snapshot (full giotto object state)
snapshotSave(src, gobject, name = "my_analysis")

# Load a snapshot (most recent if name omitted)
gobject <- snapshotLoad(src)

# Delete a snapshot
snapshotDelete(src, "my_analysis")

# Prune unprotected artifacts
sourcePrune(src)
```

### Store Types

#### `parquetStore` - Tabular Data

``` r

# For general tabular data with row indexing
store <- parquetStore(path = "data.parquet")
storeWrite(store, my_dataframe)

# Read as arrow query or tibble
query <- storeRead(store, output = "query")
df <- storeRead(store, fields = c("gene", "counts"), output = "tibble")
```

#### `parquetGeomStore` - Spatial Data

``` r

# For spatial features with extent tracking
store <- parquetGeomStore(path = "spatial_data")
storeWrite(store, polygons)  # SpatVector

# Spatial filtering happens in arrow (no full load)
roi <- ext(c(100, 200, 300, 400))
subset <- storeRead(store, extent = roi, output = "sf")
```

#### `parquetGeomTileStore` - Large Spatial Datasets

``` r

# For datasets too large for memory
store <- parquetGeomTileStore(path = "tiled_storage")

# Write with automatic tiling
storeWrite(store,
    data = large_dataset,
    n_tiles = 100,
    type = "polygon"
)

# Read specific tiles or spatial regions
tile_data <- storeRead(store, tile = 5, output = "terra")
region_data <- storeRead(store, extent = roi, output = "sf")
```

#### Matrix Stores

``` r

# HDF5-backed matrices (via HDF5Array)
store <- h5ArrayStore(path = "matrix.h5")

# BPCells on-disk matrices
store <- bpcMatrixStore(path = "matrix_dir")

# TileDB-backed matrices
store <- tileDBArrayStore(path = "matrix.tiledb")
```

#### `fileStore` - Custom Formats

``` r

# For non-standard file formats
store <- fileStore(
    path = "data.rds",
    read_fun = readRDS
)

data <- storeRead(store)
```

## Architecture

### Storage Hierarchy

    dataStore (VIRTUAL)
    └── fileStore - Disk-backed storage
        ├── queryableStore - Supports lazy queries
        │   └── parquetStore - Apache Parquet format
        │       └── parquetGeomStore - Spatial parquet with extent
        │           └── parquetGeomTileStore - Tiled spatial parquet
        ├── h5ArrayStore - HDF5-backed arrays
        ├── tileDBArrayStore - TileDB-backed arrays
        └── bpcMatrixStore - BPCells on-disk matrices

### Source Hierarchy

    gsource (VIRTUAL)
    └── gDirSource - Directory-based project management
        ├── Artifact tracking via JSON manifest
        ├── Snapshot save/load/delete
        └── Pruning with dependency protection

### Key Methods

| Method                           | Description                              |
|----------------------------------|------------------------------------------|
| `storeRead(store, ...)`          | Query data with optional filtering       |
| `storeWrite(store, data, ...)`   | Write data to store                      |
| `sourceWrite(src, data, ...)`    | Write data as a tracked artifact         |
| `snapshotSave(src, x, ...)`      | Save a giotto object snapshot            |
| `snapshotLoad(src, ...)`         | Load a snapshot                          |
| `snapshotDelete(src, name, ...)` | Delete a snapshot                        |
| `sourcePrune(src)`               | Remove unprotected artifacts             |
| `resolveSource(path)`            | Detect and regenerate a source from path |

## Advanced Features

### Lazy Operations

Filtering, sampling, and deduplication are recorded lazily in the store
and applied together at read time — no data is touched until
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
or a terminal step like
[`as.vector()`](https://rdrr.io/r/base/vector.html).

``` r

# ops compose — only one parquet scan happens at storeRead()
x |>
    subset(cell_type == "neuron") |>
    unique() |>
    head(100) |>
    storeRead(output = "tibble")
```

Column selection with `[, j]` is also lazy:

``` r

# select columns first, then deduplicate
unique(x[, c("feat_ID", "cell_type")]) |> storeRead(output = "tibble")
```

For extracting a single column, use
[`as.vector()`](https://rdrr.io/r/base/vector.html) — it requires
exactly one column selected and returns a named list, matching base R’s
[`as.vector()`](https://rdrr.io/r/base/vector.html) on a data.frame:

``` r

# distinct gene names — returns list(feat_ID = c(...))
unique(x[, "feat_ID"]) |> as.vector()
```

This keeps pipelines interchangeable between in-memory data.frames and
parquet stores. The distinct step is safe even at billions of rows since
the unique values of an ID column are always small.

### Callback for Custom Queries

``` r

# Apply custom operations to arrow queries before output
df <- storeRead(store,
    output = "tibble",
    callback = function(query) {
        query |>
            dplyr::filter(gene == "EPCAM") |>
            dplyr::mutate(log_count = log1p(count))
    }
)
```

### Hive Partitioning

``` r

# Tiled stores use hive partitioning automatically during storeWrite
# On disk: storage/tile_index=001/, storage/tile_index=002/, etc.
# Arrow pushdown automatically prunes partitions during reads
```

### Artifact Protection

``` r

# Snapshots protect their referenced artifacts from pruning
snapshotSave(src, gobject, name = "checkpoint_1")

# Only artifacts not referenced by any snapshot are pruned
sourcePrune(src)
```

## Options

### Default Store Types

When writing data via
[`sourceWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite.md),
the store format is chosen by global options. These can be set to change
the default backend for each data type:

| Option | Default | Description |
|----|----|----|
| `giotto.gdsrc_matrix_format` | `"h5"` | Store type for matrices. Alternatives: `"tiledb"`, `"bpcells"` |
| `giotto.gdsrc_spatvector_format` | `"parquetGeom"` | Store type for spatial vectors |
| `giotto.gdsrc_dataframe_format` | `"parquet"` | Store type for tabular data |

``` r

# Use BPCells for matrix storage instead of HDF5
options(giotto.gdsrc_matrix_format = "bpcells")

# All subsequent sourceWrite() calls for matrices will use BPCells
sourceWrite(src, my_matrix)
```

### Other Options

| Option | Default | Description |
|----|----|----|
| `giottodisk.use_locking` | `TRUE` | Use {filelock} for concurrent manifest access |
| `giottodisk.uid_include_pid` | `TRUE` | Include process ID in artifact UIDs |
| `giottodisk.uid_include_node` | `FALSE` | Include node name in artifact UIDs |

## Performance Tips

1.  **Filter early**: Use `extent` and `fields` parameters to reduce
    data loading
2.  **Leverage tiles**: For \>1M features, use `parquetGeomTileStore`
3.  **Lazy evaluation**: Keep data as arrow queries until final
    `collect()`
4.  **Use callbacks**: Apply custom filtering within
    [`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
    to stay lazy

## Documentation

In-depth docs live in `vignettes/articles/`:

- [`design.Rmd`](https://giotto-suite.github.io/GiottoDisk/vignettes/articles/design.Rmd)
  — architectural rationale: class hierarchy, lazy op system, transform
  back-projection, spatial predicates (`spat_relate`) and the
  multi-engine narrow path, manifest design, extent tracking, scale
  targets.
- [`gsource.Rmd`](https://giotto-suite.github.io/GiottoDisk/vignettes/articles/gsource.Rmd)
  — `gDirSource` walkthrough (the directory-backed `gsource` backend,
  currently the only one shipped): source verbs (`sourceWrite` /
  `sourceContains` / `sourceAdopt` / `sourcePrune`), snapshot lifecycle,
  deployment patterns (local, HPC/NFS, shared analyses).
- [`roadmap.Rmd`](https://giotto-suite.github.io/GiottoDisk/vignettes/articles/roadmap.Rmd)
  — public direction.
- [`parquetEdgeStore.Rmd`](https://giotto-suite.github.io/GiottoDisk/vignettes/articles/parquetEdgeStore.Rmd)
  — edge-store (graph) specifics.

For contributors and LLM-assisted work,
[`AGENTS.md`](https://giotto-suite.github.io/GiottoDisk/AGENTS.md) is
the navigation entry point — package structure, class hierarchy,
constraints and invariants, op-type dispatch, output formats,
conventions. It links out to the vignettes for architectural deep dives.

## See Also

- [Giotto Suite](https://github.com/drieslab/Giotto) - Spatial genomics
  analysis
- [tilework](https://github.com/drieslab/tilework) - Tiling framework
- [Apache Arrow](https://arrow.apache.org/) - Columnar data format
- [terra](https://github.com/rspatial/terra) - Spatial data operations
