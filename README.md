# GiottoDisk

Memory-efficient storage and processing framework for large-scale spatial
genomics data.

## Overview

{GiottoDisk} provides tools for on-disk data processing and ETL. It
leverages {arrow} and {dplyr} for their powerful out-of-memory and lazy
operations, and extends them with escape hatches for chunked in-memory 
operations when lazy methods are not available (e.g., geometry construction,
custom transformations).

This is powered through a unified storage abstraction layer for file schema,
disk reading and writing, enabling clean and maintainable ETL from
heterogeneous file formats characteristic of the spatial omics field to lazy
queries and/or file-backed representations.



## Installation
```r
# Install from GitHub
remotes::install_github("drieslab/GiottoDisk")
```

## Quick Start

### Basic Usage
```r
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

### Store Types

#### `parquetStore` - Tabular Data
```r
# For general tabular data with row indexing
store <- parquetStore(path = "data.parquet")
storeWrite(store, my_dataframe)

# Read as arrow query or tibble
query <- storeRead(store, output = "query")
df <- storeRead(store, fields = c("gene", "counts"), output = "tibble")
```

#### `parquetGeomStore` - Spatial Data
```r
# For spatial features with extent tracking
store <- parquetGeomStore(path = "spatial_data")
storeWrite(store, polygons)  # SpatVector

# Spatial filtering happens in arrow (no full load)
roi <- ext(c(100, 200, 300, 400))
subset <- storeRead(store, extent = roi, output = "sf")
```

#### `parquetGeomTileStore` - Large Spatial Datasets
```r
# For datasets too large for memory
store <- parquetGeomTileStore(path = "tiled_storage")

# Write with automatic tiling
storeWrite(store, 
    data = large_fileStore,  # Could be CSV, HDF5, etc.
    n_tiles = 100,
    type = "polygon"
)

# Read specific tiles or spatial regions
tile_data <- storeRead(store, tile = 5, output = "terra")
region_data <- storeRead(store, extent = roi, output = "sf")
```

#### `fileStore` - Custom Formats
```r
# For non-standard file formats
h5_store <- fileStore(
    path = "data.h5ad",
    read_fun = function(path) {
        anndata <- reticulate::py$scanpy$read_h5ad(path)
        # Return data.frame-like object
        as.data.frame(anndata$obs)
    }
)

data <- storeRead(h5_store)
```


## Architecture

### Storage Hierarchy
```
dataStore (VIRTUAL)
├── memoryStore - In-memory storage
└── fileStore - Disk-backed storage
    └── parquetStore - Apache Parquet format
        └── parquetGeomStore - Spatial parquet with extent
            └── parquetGeomTileStore - Tiled spatial parquet
```

### Key Methods

- `storeRead(store, ...)` - Query data with optional filtering
- `storeWrite(store, data, ...)` - Write data to store
- `nrow(store)` - Get row count without loading data
- `colnames(store)` - Get column names
- `ext(store)` - Get spatial extent (for geometric stores)
- `crop(store, extent)` - Crop to spatial region

## Advanced Features

### Batch Processing with Arrow
```r
# Efficient batch processing without full materialization
store <- parquetStore("large_table")

# Transform in batches using arrow
transformed <- storeWrite(
    new_store,
    data = storeRead(store),  # Arrow query
    callback = function(batch) {
        # Applied to each batch separately
        batch$normalized <- batch$counts / sum(batch$counts)
        batch
    }
)
```

### Hive Partitioning
```r
# Tiled stores use hive partitioning automatically
store <- parquetGeomTileStore("storage")
# Creates: storage/tile_index=001/, storage/tile_index=002/, etc.
```

## Performance Tips

1. **Filter early**: Use `extent` and `fields` parameters to reduce data loading
2. **Leverage tiles**: For >1M features, use `parquetGeomTileStore`
3. **Lazy evaluation**: Keep data as arrow queries until final `collect()`
4. **Batch processing**: Use `callback` parameter for memory-efficient transformations

## Dependencies

- **arrow** - Columnar data processing
- **terra** - Spatial vector operations
- **sf** - Simple features (optional)
- **dplyr** - Data manipulation
- **data.table** - Efficient data.frame operations

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.


## See Also

- [Giotto Suite](https://github.com/drieslab/Giotto) - Spatial genomics analysis
- [Apache Arrow](https://arrow.apache.org/) - Columnar data format
- [terra](https://github.com/rspatial/terra) - Spatial data operations
