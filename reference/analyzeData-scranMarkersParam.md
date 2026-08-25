# Streaming pairwise marker detection

[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
method for a
[Giotto::scranMarkersParam](https://rdrr.io/pkg/Giotto/man/analyze_param.html)
on a disk-backed expression store.

The per-group statistic pass is a grouped aggregate pushed into Acero,
so the store is scanned once and the pairwise comparisons that follow
never touch it. Results are verified elementwise against `findMarkers`.

`test_type` must be `"t"`: it is the only test that reduces to per-group
moments. `"wilcox"` needs the per-cell values to rank, and `"binom"` is
not wired up.

For those, fall back to the in-memory path. That means materializing the
whole store, **so it works only if the whole matrix fits in memory**:

    m <- storeRead(x, output = "dgcmatrix", max_rows = Inf, max_cols = Inf)
    analyzeData(m, markersParam(method = "scran", test_type = "wilcox"),
                groups = groups)

The cap arguments are required:
[`storeRead()`](https://giotto-suite.github.io/GiottoDisk/reference/storeRead.md)
refuses a full-size materialization unless one axis is small, which is
what forces that fit to be considered rather than discovered. Narrowing
with `[` is not an alternative here — a rank test needs every gene and
every cell.

## Usage

``` r
# S4 method for class 'parquetExprBase,scranMarkersParam'
analyzeData(x, param, ..., groups = NULL)
```

## Arguments

- x:

  a `parquetExprBase` store.

- param:

  a
  [Giotto::scranMarkersParam](https://rdrr.io/pkg/Giotto/man/analyze_param.html).

- ...:

  additional arguments (none used).

- groups:

  vector of cluster assignments, one per cell of the store's current
  view. `NA` excludes a cell from every group.

## Value

A `SimpleList` of `DataFrame`s, one per group, as `findMarkers` returns.

## See also

[`Giotto::markersParam()`](https://rdrr.io/pkg/Giotto/man/analyze_param.html),
[`Giotto::findScranMarkers()`](https://rdrr.io/pkg/Giotto/man/findScranMarkers.html)
