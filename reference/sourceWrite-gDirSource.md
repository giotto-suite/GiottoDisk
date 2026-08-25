# Write Data to a `gDirSource`

Write data to a
[gDirSource](https://giotto-suite.github.io/GiottoDisk/reference/gDirSource.md)
managed project directory. Defaults for different data types are
settable via global options:

- *sparse matrices:* `giotto.gdsrc_sparsematrix_format` (default =
  "parquetExpr")

- *dense matrices:* `"giotto.gdsrc_densematrix_format"` (default = "h5")

- *dataframes:* `"giotto.gdsrc_dataframe_format"` (default = "parquet")

- *spatvector:* `"giotto.gdsrc_spatvector_format"` (default =
  "parquetGeom")

giotto dispatch: promotes an in-memory `giotto` object to a disk-backed
one. Builds a fresh `giotto` with `src` attached as the backend, then
re-attaches every data subobject via
[`setGiotto()`](https://giotto-suite.github.io/GiottoClass/reference/setGiotto.html).
Each backend-aware setter (setExpression / setPolygonInfo /
setFeatureInfo / setNearestNetwork / setSpatialNetwork) routes its
in-memory subobject to the corresponding disk-backed store
(parquetExprStore / parquetGeomStore / parquetEdgeStore).

Subobjects whose setters don't (yet) have backend-aware write logic stay
in-memory but are still attached to the new gobject.

## Usage

``` r
# S4 method for class 'gDirSource,memoryMatrix'
sourceWrite(src, data, meta = NULL, store_type = NULL, ...)

# S4 method for class 'gDirSource,SpatVector'
sourceWrite(
  src,
  data,
  meta = NULL,
  store_type = getOption("giotto.gdsrc_spatvector_format", "parquetGeom"),
  ...
)

# S4 method for class 'gDirSource,data.frame'
sourceWrite(
  src,
  data,
  meta = NULL,
  store_type = getOption("giotto.gdsrc_dataframe_format", "parquet"),
  ...
)

# S4 method for class 'gDirSource,ANY'
sourceWrite(src, data, meta = NULL, ...)

# S4 method for class 'gDirSource,igraph'
sourceWrite(src, data, meta = NULL, ...)

# S4 method for class 'gDirSource,giotto'
sourceWrite(src, data, ...)
```

## Arguments

- src:

  `gsource` object. Used to provide default save locations and save
  formats for a managed Giotto backend.

- data:

  data to write

- meta:

  `list`. Additional metadata to attach to this Giotto backend

- store_type:

  `character`. Store type to write as

- ...:

  additional params passed to `storeWrite` method managed artifact.

## Value

written `store` object

`giotto` with `@source = src` and backable subobjects disk-backed.

## See also

Other sourceWrite methods:
[`sourceWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite.md)

## Examples

``` r
testdir <- file.path(tempdir(), "testdirsource")
gsrc <- sourceCreate(path = testdir)
# arrays ----------------------------------------- #
m <- matrix(1:9, nrow = 3)
m_written <- sourceWrite(gsrc, m)
storeRead(m_written)
# spatvector ----------------------------------------- #
sv <- terra::vect(system.file("ex/lux.shp", package="terra"))
sv_written <- sourceWrite(gsrc, sv)
storeRead(sv_written) # default is output = "query"
storeRead(sv_written, output = "sf")
storeRead(sv_written, output = "terra")
storeRead(sv_written, output = "tibble")
# tables ----------------------------------------- #
store <- sourceWrite(gsrc, iris)
force(store)
storeRead(store) # default output is an arrow query
storeRead(store, output = "tibble") # pull into memory as tibble
```
