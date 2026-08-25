# Get Bounded Data

Get the region of the data defined by a bounding `SpatExtent`.
[`crop()`](https://giotto-suite.github.io/GiottoDisk/reference/crop.md)
operations on `fileStore` objects return `dataStore`, but this function
directly produces a usable interface.

## Usage

``` r
# S4 method for class 'queryableStore,SpatExtent'
getBoundedData(
  x,
  bound,
  sdimx = "x",
  sdimy = "y",
  inclusive = TRUE,
  group_col = "id",
  envelope = FALSE,
  output = "query",
  ...
)

# S4 method for class 'parquetGeomStore,SpatExtent'
getBoundedData(
  x,
  bound,
  sdimx = "x_index",
  sdimy = "y_index",
  inclusive = TRUE,
  ...
)
```

## Arguments

- x:

  `queryableStore` inheriting object

- bound:

  `SpatExtent`

- sdimx:

  `character`. Which column contains spatial x values to filter on

- sdimy:

  `character`. Which column contains spatial y values to filter on.

- inclusive:

  `logical` (length 1 or 4. Default = TRUE). Whether bounds should be
  inclusive (inclusive is \>=/\<= vs \>/\<). Ordering of inputs is
  bottom, left, top, right. If a single logical is provided, it will be
  replicated and used for all 4 bounds. Default is fully inclusive.

- group_col:

  `character`. Only used when `envelope = TRUE`. Which column contains
  grouping IDs that are used to calculate the spatial envelope.

- envelope:

  `logical` (default = FALSE). Whether to perform bound filtering based
  on the envelope centroid instead of individual x and y values so that
  grouped rows are not separated by the filtering.

  This is generally not needed for `parquetGeomStore` inheriting objects
  since each row is its own geometry (point or polygon, etc) and xy
  filtering is performed on the point or the centroid.

- output:

  `character`. Output format passed to
  [`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md).

- ...:

  additional params to pass to
  [`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)

## Value

output as determined by `output` param – see
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)

## See also

Other tile methods:
[`getTile`](https://giotto-suite.github.io/GiottoDisk/reference/getTile.md)
