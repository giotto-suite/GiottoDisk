# Prune a project source

On-disk processing can generate many intermediate artifacts. Manually
call this function in order to drop artifacts that are not tagged with a
giottosave that depends on them (the giottosave .rds file must also
exist).

## Usage

``` r
# S4 method for class 'gDirSource'
sourcePrune(src, ...)

# S4 method for class 'giotto'
sourcePrune(src, ...)
```

## Arguments

- src:

  `gsource` object. Used to provide default save locations and save
  formats for a managed Giotto backend.
