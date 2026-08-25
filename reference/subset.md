# Subset a parquet store

Filter rows of a parquet store by a logical expression. The operation is
recorded lazily in the store's `@ops` slot and applied at read time.

Expressions are evaluated using non-standard evaluation. Local variables
referenced in the expression are automatically inlined at capture time,
making the recorded operation self-contained and safe to pass across
sessions or to parallel workers.

## Usage

``` r
# S4 method for class 'parquetBase'
subset(x, subset, select, negate = FALSE, quote = TRUE, ...)
```

## Arguments

- x:

  `parquetBase`-inheriting store object

- subset:

  logical `expression` to filter rows. Column names in the store are
  referenced directly. Local variables are inlined automatically – no
  `!!` injection needed.

- select:

  `expression`, indicating columns to keep. `-` can be used to drop
  columns.

- negate:

  `logical`. If `TRUE`, the filter expression is negated, keeping
  rows/cols that do NOT satisfy the condition.

- quote:

  `logical`. If `TRUE` (default), `subset` is captured via NSE. If
  `FALSE`, `subset` is treated as a pre-built R `call` object, allowing
  programmatic construction of filter expressions.

- ...:

  additional arguments (ignored)

## Value

the store with the filter step appended to `@ops`

## Examples

``` r
# standard NSE
subset(store, gene == "EPCAM")

# local variable -- inlined automatically
my_genes <- c("EPCAM", "CDH1")
subset(store, gene %in% my_genes)

# negation
subset(store, gene == "EPCAM", negate = TRUE)

# pre-built expression
expr <- quote(gene == "EPCAM")
subset(store, expr, quote = FALSE)
```
