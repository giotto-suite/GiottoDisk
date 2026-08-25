# Write Data to a Source

Write data to a
[gsource](https://giotto-suite.github.io/GiottoDisk/reference/gsource.md)
that manages data formats and file organization. Docs for specific
`gsource` implementations have further information.

- [sourceWrite-gDirSource](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)

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

## See also

Other sourceWrite methods:
[`sourceWrite-gDirSource`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite-gDirSource.md)
