# Delete a Giotto Snapshot

Delete a Giotto snapshot. By extension, also removes the protection of
the associated giottosave tag for artifact pruning.

## Usage

``` r
# S4 method for class 'character,character'
snapshotDelete(src, name, ...)

# S4 method for class 'gDirSource,character'
snapshotDelete(src, name, ...)
```

## Arguments

- src:

  `gsource` object or `character` filepath to project dir if
  `gDirSource` controlled.

- name:

  `character` (optional) name of specific snapshot to delete. if NULL,
  the most recent is selected

- ...:

  additional params to pass (none implemented)

## Value

`TRUE` invisibly if succeeds
