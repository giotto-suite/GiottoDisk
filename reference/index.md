# Package index

## Sources

A source owns a project directory and the stores inside it.

- [`gsource`](https://giotto-suite.github.io/GiottoDisk/reference/gsource.md)
  : Giotto Sources

- [`gDirSource()`](https://giotto-suite.github.io/GiottoDisk/reference/gDirSource.md)
  : Create a Giotto Directory Source

- [`gDirSource-class`](https://giotto-suite.github.io/GiottoDisk/reference/gDirSource-class.md)
  : Giotto Directory Source

- [`sourceCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceCreate.md)
  : Create a Giotto Source

- [`sourceWrite`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite.md)
  : Write Data to a Source

- [`sourceWrite(`*`<gDirSource>`*`,`*`<memoryMatrix>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)
  [`sourceWrite(`*`<gDirSource>`*`,`*`<SpatVector>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)
  [`sourceWrite(`*`<gDirSource>`*`,`*`<data.frame>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)
  [`sourceWrite(`*`<gDirSource>`*`,`*`<ANY>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)
  [`sourceWrite(`*`<gDirSource>`*`,`*`<igraph>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)
  [`sourceWrite(`*`<gDirSource>`*`,`*`<giotto>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)
  :

  Write Data to a `gDirSource`

- [`sourceAdopt(`*`<gDirSource>`*`,`*`<fileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceAdopt.md)
  [`sourceAdopt(`*`<gDirSource>`*`,`*`<parquetEdgeStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceAdopt.md)
  [`sourceAdopt(`*`<gDirSource>`*`,`*`<SpatRaster>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceAdopt.md)
  [`sourceAdopt(`*`<gDirSource>`*`,`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceAdopt.md)
  [`sourceAdopt(`*`<gDirSource>`*`,`*`<ANY>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceAdopt.md)
  : Adopt a Store into a Source

- [`sourceContains(`*`<gDirSource>`*`,`*`<fileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceContains.md)
  [`sourceContains(`*`<gDirSource>`*`,`*`<SpatRaster>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceContains.md)
  [`sourceContains(`*`<gDirSource>`*`,`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceContains.md)
  [`sourceContains(`*`<gDirSource>`*`,`*`<ANY>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourceContains.md)
  : Test if a Store is Managed by a Source

- [`sourcePrune(`*`<gDirSource>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourcePrune.md)
  [`sourcePrune(`*`<giotto>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/sourcePrune.md)
  : Prune a project source

- [`resolveSource()`](https://giotto-suite.github.io/GiottoDisk/reference/resolveSource.md)
  : Source Detection and Regeneration

- [`as.list(`*`<gDirSource>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/coerce_gdsrc.md)
  [`as.data.frame(`*`<gDirSource>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/coerce_gdsrc.md)
  : Coerce gDirSource

- [`sd_view_ref()`](https://giotto-suite.github.io/GiottoDisk/reference/sd_view_ref.md)
  : Get the registered view name for a GiottoDisk SedonaDB dataframe

## Stores

The generic store interface, shared by every backend.

- [`dataStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md)
  [`store`](https://giotto-suite.github.io/GiottoDisk/reference/dataStore-class.md)
  : Data Storage

- [`storeCreate()`](https://giotto-suite.github.io/GiottoDisk/reference/storeCreate.md)
  : Create a Store

- [`storeRead(`*`<parquetEdgeStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<mtxInput>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<tenxH5Input>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<cellbinGefInput>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<binGefInput>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<csvWideInput>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<parquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<unionParquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<fileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<queryableStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<parquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<unionParquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<parquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<parquetGeomTileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<h5ArrayStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<tileDBArrayStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  [`storeRead(`*`<bpcMatrixStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
  :

  Read a `dataStore`

- [`storeWrite(`*`<parquetEdgeStore>`*`,`*`<edgeInput>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<parquetEdgeStore>`*`,`*`<data.table>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<parquetEdgeStore>`*`,`*`<igraph>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<parquetExprStore>`*`,`*`<exprInput>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<parquetExprStore>`*`,`*`<parquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<parquetExprStore>`*`,`*`<unionParquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<parquetExprStore>`*`,`*`<memoryMatrix>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<fileStore>`*`,`*`<fileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<h5ArrayStore>`*`,`*`<memoryMatrix>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<tileDBArrayStore>`*`,`*`<memoryMatrix>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<bpcMatrixStore>`*`,`*`<memoryMatrix>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<bpcMatrixStore>`*`,`*`<ANY>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<bpcMatrixStore>`*`,`*`<mtxInput>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  [`storeWrite(`*`<bpcMatrixStore>`*`,`*`<tenxH5Input>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  :

  Write to a `dataStore`

- [`storeWrite(`*`<parquetStore>`*`,`*`<data.frame>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetStore.md)
  [`storeWrite(`*`<parquetStore>`*`,`*`<ANY>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetStore.md)
  : Write to a Parquet Storage Spec

- [`storeWrite(`*`<parquetGeomStore>`*`,`*`<SpatVector>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomStore.md)
  [`storeWrite(`*`<parquetGeomStore>`*`,`*`<data.frame>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomStore.md)
  : Write to a Parquet Geometry Storage

- [`storeWrite(`*`<parquetGeomTileStore>`*`,`*`<queryableStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomTileStore.md)
  [`storeWrite(`*`<parquetGeomTileStore>`*`,`*`<parquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomTileStore.md)
  [`storeWrite(`*`<parquetGeomTileStore>`*`,`*`<parquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite-parquetGeomTileStore.md)
  : Write to a Parquet Geometry Tiled Storage

- [`storeExists(`*`<fileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeExists.md)
  [`storeExists(`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeExists.md)
  : Store Existence

- [`storePaths(`*`<fileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storePaths.md)
  [`storePaths(`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storePaths.md)
  : Store Paths

- [`storeUID(`*`<fileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeUID.md)
  [`storeUID(`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeUID.md)
  [`storeUID(`*`<overlapPointDisk>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeUID.md)
  : Store UIDs

- [`storeChunkInfo(`*`<parquetExprBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/storeChunkInfo.md)
  : Report how a store's streaming windows are chosen

- [`specialCols(`*`<ANY>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/specialCols.md)
  [`specialCols(`*`<parquetBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/specialCols.md)
  [`specialCols(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/specialCols.md)
  [`specialCols(`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/specialCols.md)
  : Special Store Columns

## Parquet stores

Parquet-backed expression, geometry and edge storage.

- [`parquetStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md)
  [`parquetGeomStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md)
  [`parquetGeomTileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore.md)
  : Create a Parquet Store
- [`parquetStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md)
  [`parquetGeomStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md)
  [`parquetGeomTileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetStore-class.md)
  : Parquet Store
- [`parquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore.md)
  : Create a Parquet Expression Matrix Store
- [`parquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md)
  [`parquetExprBase-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetExprStore-class.md)
  : Parquet Expression Matrix Store (streaming)
- [`parquetEdgeStore()`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore.md)
  : Create a Parquet Edge Store handle
- [`parquetEdgeStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/parquetEdgeStore-class.md)
  : Parquet Edge Store (streaming)
- [`unionParquetExprStore()`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
  [`cbind2(`*`<parquetExprStore>`*`,`*`<parquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
  [`cbind2(`*`<unionParquetExprStore>`*`,`*`<parquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
  [`cbind2(`*`<parquetExprStore>`*`,`*`<unionParquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
  [`cbind2(`*`<unionParquetExprStore>`*`,`*`<unionParquetExprStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore.md)
  : Construct a unionParquetExprStore
- [`unionParquetExprStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/unionParquetExprStore-class.md)
  : Virtual Union of Expression Stores

## Other backends

- [`fileStore()`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore.md)
  : Create a File Store
- [`fileStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md)
  [`queryableStore-class`](https://giotto-suite.github.io/GiottoDisk/reference/fileStore-class.md)
  : File Store
- [`h5ArrayStore()`](https://giotto-suite.github.io/GiottoDisk/reference/arrayStore.md)
  [`tileDBArrayStore()`](https://giotto-suite.github.io/GiottoDisk/reference/arrayStore.md)
  [`bpcMatrixStore()`](https://giotto-suite.github.io/GiottoDisk/reference/arrayStore.md)
  : Array Storage

## Inputs

Wrappers describing on-disk or in-memory data to be ingested.

- [`exprInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/exprInput-class.md)
  : Expression Matrix Input (virtual)
- [`mtxInput()`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput.md)
  : Create a 10x / Xenium MatrixMarket triple input
- [`mtxInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/mtxInput-class.md)
  : 10x / Xenium MatrixMarket Triple Input
- [`tenxH5Input()`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input.md)
  : Create a 10x cell_feature_matrix.h5 input
- [`tenxH5Input-class`](https://giotto-suite.github.io/GiottoDisk/reference/tenxH5Input-class.md)
  : 10x cell_feature_matrix.h5 Input
- [`csvWideInput()`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput.md)
  : Create a wide-format CSV input
- [`csvWideInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/csvWideInput-class.md)
  : Wide-format CSV Input
- [`binGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput.md)
  : Create a Stereo-seq bin GEF input
- [`binGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/binGefInput-class.md)
  : Stereo-seq bin GEF Input
- [`cellbinGefInput()`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput.md)
  : Create a Stereo-seq cellbin GEF input
- [`cellbinGefInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/cellbinGefInput-class.md)
  : Stereo-seq cellbin GEF Input
- [`edgeInput-class`](https://giotto-suite.github.io/GiottoDisk/reference/edgeInput-class.md)
  : Base class for edge-table inputs
- [`edgeDTInput()`](https://giotto-suite.github.io/GiottoDisk/reference/edgeDTInput.md)
  : Wrap an in-memory data.table as an edge input
- [`igraphInput()`](https://giotto-suite.github.io/GiottoDisk/reference/igraphInput.md)
  : Wrap an in-memory igraph as an edge input
- [`nnSearchInput()`](https://giotto-suite.github.io/GiottoDisk/reference/nnSearchInput.md)
  : Wrap an RANN / dbscan kNN result as an edge input

## Disk-backed importers

Read a vendor run straight to disk, without loading it into memory.

- [`importXeniumDisk()`](https://giotto-suite.github.io/GiottoDisk/reference/importXeniumDisk.md)
  : Import a 10X Xenium assay (disk-backed)
- [`importCosMxDisk()`](https://giotto-suite.github.io/GiottoDisk/reference/importCosMxDisk.md)
  : Import a NanoString CosMx assay (disk-backed)
- [`importStereoSeqDisk()`](https://giotto-suite.github.io/GiottoDisk/reference/importStereoSeqDisk.md)
  : Import a BGI Stereo-seq assay (disk-backed)

## Geometry

- [`createGiottoPoints(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/createGiottoPoints.md)
  : Create a giottoPoints from a Parquet Geometry Store

- [`createGiottoPolygon(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/createGiottoPolygon.md)
  : Create a giottoPolygon from a Parquet Geometry Store

- [`centroids(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/centroids.md)
  : Get Polygon Centroids from a Parquet Geometry Store

- [`expanse(`*`<parquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/expanse.md)
  [`expanse(`*`<parquetGeomTileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/expanse.md)
  : Get the area of individual polygons

- [`ext(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/ext.md)
  : Spatial Extent

- [`crop(`*`<parquetGeomBase>`*`,`*`<ANY>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/crop.md)
  :

  Crop a `parquetGeomStore`

- [`window(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/window.md)
  [`` `window<-`( ``*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/window.md)
  :

  Window a `parquetGeomStore`

- [`getBoundedData(`*`<queryableStore>`*`,`*`<SpatExtent>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/getBoundedData.md)
  [`getBoundedData(`*`<parquetGeomStore>`*`,`*`<SpatExtent>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/getBoundedData.md)
  : Get Bounded Data

- [`getTile(`*`<queryableStore>`*`,`*`<spatialTilePlan>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/getTile.md)
  [`getTile(`*`<queryableStore>`*`,`*`<freeTilePlan>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/getTile.md)
  : Get Tile Data

- [`rasterize(`*`<parquetGeomStore>`*`,`*`<SpatRaster>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rasterize.md)
  : Rasterize Parquet Point Data

- [`calculateOverlap(`*`<parquetGeomStore>`*`,`*`<parquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/calculateOverlap.md)
  [`calculateOverlap(`*`<parquetGeomStore>`*`,`*`<parquetGeomTileStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/calculateOverlap.md)
  : Calculate Overlap

- [`overlapToMatrix(`*`<overlapPointDisk>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/overlapToMatrix.md)
  [`overlapToMatrix(`*`<parquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/overlapToMatrix.md)
  : Aggregate Overlap Results to Sparse Matrix

## Spatial transforms

- [`affine(`*`<parquetGeomBase>`*`,`*`<affine2d>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/affine.md)
  : Lazily record an affine transform on a parquetGeomBase store
- [`flip(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/flip.md)
  : Flip a parquetGeomBase store
- [`spin(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/spin.md)
  : Rotate a parquetGeomBase store
- [`shear(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/shear.md)
  : Shear a parquetGeomBase store
- [`rescale(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rescale.md)
  : Scale a parquetGeomBase store
- [`spatShift(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/spatShift.md)
  : Translate a parquetGeomBase store
- [`t(`*`<parquetGeomBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/t.md)
  : Transpose a parquetGeomBase store

## Subset, combine and inspect

- [`subset(`*`<parquetBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/subset.md)
  : Subset a parquet store
- [`rowSample(`*`<parquetBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rowSample.md)
  : Subsample rows of a parquetStore
- [`rbind2(`*`<parquetStore>`*`,`*`<parquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  [`rbind2(`*`<parquetStore>`*`,`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  [`rbind2(`*`<unionParquetStore>`*`,`*`<parquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  [`rbind2(`*`<unionParquetStore>`*`,`*`<unionParquetStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  [`rbind2(`*`<parquetGeomStore>`*`,`*`<parquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  [`rbind2(`*`<unionParquetGeomStore>`*`,`*`<parquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  [`rbind2(`*`<parquetGeomStore>`*`,`*`<unionParquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  [`rbind2(`*`<unionParquetGeomStore>`*`,`*`<unionParquetGeomStore>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/rbind.md)
  : rbind stores
- [`head(`*`<parquetBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/parquet-headtail.md)
  [`tail(`*`<parquetBase>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/parquet-headtail.md)
  : Head and tail
- [`plot(`*`<parquetGeomStore>`*`,`*`<missing>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/plot.md)
  [`plot(`*`<sedonadb_dataframe>`*`,`*`<missing>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/plot.md)
  : Visualize a Store

## Streaming analysis

Parameter objects driving out-of-memory computation.

- [`analyzeData(`*`<parquetExprBase>`*`,`*`<featStatsParam>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/analyzeData-featStatsParam.md)
  : Streaming per-feature statistics
- [`analyzeData(`*`<parquetExprBase>`*`,`*`<scranMarkersParam>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/analyzeData-scranMarkersParam.md)
  : Streaming pairwise marker detection
- [`gramEigenPcaParam()`](https://giotto-suite.github.io/GiottoDisk/reference/gramEigenPcaParam.md)
  : Gram-eigen streaming PCA parameter

## Snapshots and artifacts

- [`snapshotSave(`*`<gDirSource>`*`,`*`<giotto>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/snapshotSave.md)
  : Write Giotto Snapshot
- [`snapshotLoad(`*`<character>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/snapshotLoad.md)
  [`snapshotLoad(`*`<gDirSource>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/snapshotLoad.md)
  : Load Giotto Snapshot
- [`snapshotDelete(`*`<character>`*`,`*`<character>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/snapshotDelete.md)
  [`snapshotDelete(`*`<gDirSource>`*`,`*`<character>`*`)`](https://giotto-suite.github.io/GiottoDisk/reference/snapshotDelete.md)
  : Delete a Giotto Snapshot
- [`setArtifactDumpDir()`](https://giotto-suite.github.io/GiottoDisk/reference/artifact_dump.md)
  [`getArtifactDumpDir()`](https://giotto-suite.github.io/GiottoDisk/reference/artifact_dump.md)
  : GiottoDisk artifact dump management
- [`artifact_uid`](https://giotto-suite.github.io/GiottoDisk/reference/artifact_uid.md)
  [`uid`](https://giotto-suite.github.io/GiottoDisk/reference/artifact_uid.md)
  : Artifact Unique Identifier
