# Get Tile Data

Get the region of data defined by a tilework tile\* object. This
function determines the tile(s) requested based on the tilePlan then
calls
[`getBoundedData()`](https://giotto-suite.github.io/GiottoDisk/reference/getBoundedData.md)
for the actual retrieval step.

## Usage

``` r
# S4 method for class 'queryableStore,spatialTilePlan'
getTile(x, tiles, i = NULL, j, contiguous = FALSE, pad = NULL, ...)

# S4 method for class 'queryableStore,freeTilePlan'
getTile(x, tiles, i = NULL, contiguous = FALSE, pad = NULL, ...)
```

## Arguments

- x:

  `fileStore` object (currently `parquetStore` inheriting only)

- tiles:

  tilework `tile*` object

- i:

  **ANY `tile*` except tileIterator** tile vector index or row index if
  `j` is also provided

- j:

  **ANY `tile*` except tileIterator** tile col index

- contiguous:

  `logical` (default = FALSE). Whether to retrieve tiles using
  inclusivity rules (see
  [getBoundedData](https://giotto-suite.github.io/GiottoDisk/reference/getBoundedData.md))
  that prevent double counts between neighboring tiles for values that
  fall on the boundary line.

  When used, `get_params$inclusive` settings will be ignored.

  This is mainly relevant only for padding = 0 situations with
  `tilePlan` inheriting `tile` inputs and operations that will process
  all tiles in the `tilePlan`. Inclusivity rules are calculated based on
  the entire set of tiles defined in the `tilePlan`.

- pad:

  `numeric` (optional) additional padding to apply before tile
  retrieval. Useful for temporarily increasing padding without affecting
  `tile*` object.

- ...:

  additional named params passed through to
  [`getBoundedData()`](https://giotto-suite.github.io/GiottoDisk/reference/getBoundedData.md)
  (e.g. `sdimx`, `sdimy`, `inclusive`, `envelope`, `output`).

## Value

a lazy arrow/dplyr query

## Additional params:

- `advance` **`tileIterator` only** `logical` (default = TRUE). Whether
  to advance the iterator.

## See also

[tilework::tilePlan](https://drieslab.github.io/tilework/reference/tilePlan.html)

Other tile methods:
[`getBoundedData`](https://giotto-suite.github.io/GiottoDisk/reference/getBoundedData.md)
