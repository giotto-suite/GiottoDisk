# Source Detection and Regeneration

Detect if an input path is to a known
[gsource](https://giotto-suite.github.io/GiottoDisk/reference/gsource.md)
type. If so, regenerate the source object. Otherwise, return FALSE.

## Usage

``` r
resolveSource(path)
```

## Arguments

- path:

  `character` filepath to managed directory

## Value

`gsource` inheriting object if recognized source type, `FALSE` if not.
