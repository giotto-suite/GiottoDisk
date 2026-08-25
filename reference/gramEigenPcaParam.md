# Gram-eigen streaming PCA parameter

Streaming PCA that accumulates the P × P Gram matrix `AᵀA` in two disk
passes and eigendecomposes it in memory. Forming `AᵀA` squares the
condition number; when the predicted rel-err on d_k exceeds
`fallback_relerr`, delegates to the Halko path.

## Usage

``` r
gramEigenPcaParam(
  ncp = 50L,
  center = TRUE,
  scale = TRUE,
  feats_to_use = NULL,
  set_seed = TRUE,
  seed_number = 1234L,
  n_oversamples = 10L,
  n_power_iter = 2L,
  fallback_relerr = 0.01,
  ...
)
```

## Arguments

- ncp:

  number of components. Default `50`.

- center:

  logical. Center columns. Default `TRUE`.

- scale:

  logical. Per-gene z-score (implicit — σ rescaling is absorbed into the
  projection matrix; sparsity preserved, no extra passes). Default
  `TRUE`.

- feats_to_use:

  character vector of feature IDs. `NULL` = use all.

- set_seed, seed_number:

  Halko-fallback seed control.

- n_oversamples, n_power_iter:

  Halko-fallback knobs. Defaults `10`, `2`.

- fallback_relerr:

  predicted-relerr threshold above which the method delegates to Halko.
  Default `0.01`.

- ...:

  reserved.

## Value

A `gramEigenPcaParam` object.

## Examples

``` r
# p <- gramEigenPcaParam(ncp = 30, feats_to_use = hvg_ids)
# res <- reduceData(pes, p)
```
