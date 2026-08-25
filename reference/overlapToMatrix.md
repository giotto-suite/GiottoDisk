# Aggregate Overlap Results to Sparse Matrix

Aggregates the output of
[`calculateOverlap()`](https://giotto-suite.github.io/GiottoDisk/reference/calculateOverlap.md)
into a feature x cell sparse count matrix, written as a Matrix Market
(.mtx) directory. The directory layout is 10x-compatible (`matrix.mtx`,
`barcodes.tsv`, `features.tsv`) and can be loaded directly by
`BPCells::import_matrix_market()`,
[`Matrix::readMM()`](https://rdrr.io/pkg/Matrix/man/externalFormats.html),
or Python's `scipy.io.mmread`.

To obtain a BPCells on-disk matrix pass the returned path to
`BPCells::import_matrix_market()` followed by
`BPCells::write_matrix_dir()`.

## Usage

``` r
# S4 method for class 'overlapPointDisk'
overlapToMatrix(
  x,
  name = "raw",
  sort = TRUE,
  count_col = NULL,
  store_type = getOption("giotto.gdsrc_sparsematrix_format", "parquetExpr"),
  path = .dump_tempfile(),
  output = c("store", "exprObj"),
  ...
)

# S4 method for class 'parquetStore'
overlapToMatrix(
  x,
  path = .dump_tempfile(),
  feat_id_col = "feat_ID",
  poly_id_col = "poly_ID",
  count_col = NULL,
  store_type = getOption("giotto.gdsrc_sparsematrix_format", "parquetExpr"),
  sort = TRUE,
  all_feat_ids = NULL,
  all_cell_ids = NULL,
  verbose = NULL,
  ...
)
```

## Arguments

- x:

  `overlapPointDisk` output from
  [`calculateOverlap()`](https://giotto-suite.github.io/GiottoDisk/reference/calculateOverlap.md)

- count_col:

  `character` (optional) column to sum instead of counting rows. Useful
  when feature detections carry a `count` field.

- path:

  `character` output directory for the Matrix Market files

- ...:

  additional params to pass

- feat_id_col:

  `character` feature ID column name (default `"feat_ID"`)

- poly_id_col:

  `character` polygon ID column name (default `"poly_ID"`)

## Value

`character` path to the output directory (invisibly)
