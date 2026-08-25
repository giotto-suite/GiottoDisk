# GiottoDisk artifact dump management

Utilities to control the location of automatically created GiottoDisk
artifacts. Some GiottoDisk processes perform operation that materialize
to disk. The dump is used as the default location artifacts will be
written to. When setting objects into the Giotto object, these objects
may be reclaimed from the dump and placed under control of the `gsource`
managing the object.

Passing a `gDirSource` object points the dump directly at the source's
vault (`artifacts/`). Artifacts are then written to their final vault
location immediately, making adoption a pure registration step with no
file movement.

## Usage

``` r
setArtifactDumpDir(x, verbose = NULL)

getArtifactDumpDir()
```

## Arguments

- x:

  `character` path to directory, or a `gDirSource` object. When a
  `gDirSource` is supplied, the dump is set to the vault directory
  (`<x@path>/artifacts/`). When missing, resets to the default session
  temp location.

- verbose:

  verbosity

## Value

artifact dump directory path (invisibly)
