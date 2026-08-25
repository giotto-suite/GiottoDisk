# Report how a store's streaming windows are chosen

Streaming passes read a store in windows sized to stay within a fraction
of free RAM. The window is **derived per read** rather than stored, from
the store's shape and its cached `@stats` marginals, so it adapts to the
machine doing the reading. This reports what those windows come out to
and why.

Two read shapes are shown because the package has two, and they differ
by roughly 4x in memory per stored value:

- `sparse matrix`:

  chunks land in a `dgCMatrix` (4-byte index + 8-byte value). Used by
  the PCA band loops and the
  [`storeWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/storeWrite.md)
  bake.

- `triplet frame`:

  chunks are collected as `row_id` / `col_id` / `value` / `source_id`
  plus arrow and data.table overhead, measured at 47-54 bytes per stored
  value. Used by R-side statistic accumulation when the op chain cannot
  be lowered to Acero.

## Usage

``` r
# S4 method for class 'parquetExprBase'
storeChunkInfo(
  x,
  ram_frac = c(0.1, 0.15, 0.2, 0.25, 0.35, 0.5),
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  a `parquetExprStore` or `unionParquetExprStore`.

- ram_frac:

  numeric. Fractions to tabulate. Defaults to a spread around the
  configured value.

- verbose:

  logical. Print the report (default `TRUE`).

- ...:

  unused.

## Value

A `data.table` of one row per (pass, ram_frac), invisibly.

## Steering the window

- `giottodisk.chunk_ram_frac`:

  fraction of free RAM to budget per chunk (default 0.25). Lower it on a
  busy machine; raise it on a dedicated one. Scales every pass without
  knowing any store's shape.

- `giottodisk.chunk_size`:

  pins an absolute window in rows, overriding the derivation entirely.
  Also the value used when the derivation cannot run.

## Examples

``` r
if (FALSE) { # \dontrun{
storeChunkInfo(pe)
options(giottodisk.chunk_ram_frac = 0.10)   # halve every streaming window
} # }
```
