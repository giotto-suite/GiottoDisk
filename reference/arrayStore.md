# Array Storage

Array Storage

## Usage

``` r
h5ArrayStore(path = .dump_tempfile(), name = HDF5Array::getHDF5DumpName(), ...)

tileDBArrayStore(
  path = file.path(tempdir(), .make_uid()),
  name = TileDBArray::getTileDBAttr(),
  ...
)

bpcMatrixStore(path = file.path(tempdir(), .make_uid()), ...)
```

## Arguments

- path:

  storage directory

- name:

  `character` Remote naming used by some storage types, for example
  HDF5Array and TileDB.
