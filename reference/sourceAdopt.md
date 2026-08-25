# Adopt a Store into a Source

Move an existing store into the managed artifact vault of a
[gsource](https://giotto-suite.github.io/GiottoDisk/reference/gsource.md),
registering it in the manifest.

- `fileStore`: the store's `@uid` is preserved so existing `source_id`
  partition references remain valid. Store must already be written to
  disk.

- `SpatRaster`: in-memory rasters are written as COG; on-disk rasters
  are moved into the vault. A fresh uid is assigned. Returns the updated
  `SpatRaster` pointing to the vault path.

- Union stores: each substore is adopted independently.

## Usage

``` r
# S4 method for class 'gDirSource,fileStore'
sourceAdopt(src, store, meta = NULL, giottosave = NULL, depends = NULL, ...)

# S4 method for class 'gDirSource,parquetEdgeStore'
sourceAdopt(src, store, meta = NULL, giottosave = NULL, depends = NULL, ...)

# S4 method for class 'gDirSource,SpatRaster'
sourceAdopt(src, store, meta = NULL, giottosave = NULL, ...)

# S4 method for class 'gDirSource,unionParquetStore'
sourceAdopt(src, store, meta = NULL, giottosave = NULL, ...)

# S4 method for class 'gDirSource,ANY'
sourceAdopt(src, store, meta = NULL, giottosave = NULL, ...)
```

## Arguments

- src:

  `gsource` object

- store:

  object to adopt (`fileStore`, `SpatRaster`, or union store)

- meta:

  `list` (optional). Additional metadata to attach to the artifact
  entry.

- giottosave:

  `character` (optional). Giottosave name to tag the artifact with.

- ...:

  additional params (unused)

## Value

updated store/raster object pointing to the managed vault location
