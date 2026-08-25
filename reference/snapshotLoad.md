# Load Giotto Snapshot

Load a gobject snapshot. If no name is provided, the most recent
snapshot is automatically loaded.

## Usage

``` r
# S4 method for class 'character'
snapshotLoad(src, ...)

# S4 method for class 'gDirSource'
snapshotLoad(src, name = NULL, load_params = list(), verbose = NULL, ...)
```

## Arguments

- src:

  `gsource` object or `character` filepath to project dir if
  `gDirSource` controlled.

- ...:

  additional params to pass (none implemented)

- name:

  `character` (optional) name of specific snapshot to load. if NULL, the
  most recent is selected

- load_params:

  `list` additional parameters (if any) for loading or reading giotto
  object

- verbose:

  verbosity

## Value

`giotto` object (possibly with additional downstream load steps needed)
