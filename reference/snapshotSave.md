# Write Giotto Snapshot

Register a gobject snapshot into the source

## Usage

``` r
# S4 method for class 'gDirSource,giotto'
snapshotSave(
  src,
  x,
  name = format(Sys.time(), "%Y%m%d_giottosave"),
  method = c("rds", "qs"),
  method_params = list(),
  overwrite = FALSE,
  export_image = TRUE,
  verbose = NULL,
  ...
)
```

## Arguments

- src:

  `gsource` object. Used to provide default save locations and save
  formats for a managed Giotto backend.

- x:

  `giotto` object to write

- method:

  `character` (default = "rds").

- method_params:

  `list` (optional) list of additional named params for save method.

- overwrite:

  `logical` (default = FALSE). Whether to overwrite if a snapshot of the
  same name exists.

- export_image:

  `logical` (default = TRUE) Whether to make a copy of external images
  in the project directory

- verbose:

  verbosity

- ...:

  additional params to pass (none implemented)

## Value

the modified gobject (invisible). Mutated by the internal adoption pass
— captured by `snapshotSave(gDirSource, giottoMulti)` so the multi-level
`.rds` sees post-adoption file handles in each child.
