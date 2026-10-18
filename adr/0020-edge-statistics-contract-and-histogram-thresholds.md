# 0020. Edge statistics cross backends as sufficient statistics; order statistics come from streamed histograms

- **Status:** Accepted
- **Date:** 2026-10-16
- **Supersedes:** —
- **Superseded by:** —

## Context

Spatial gene detection (binSpect, Moran's I, Geary's C) had only in-memory
implementations in Giotto. binSpect joined every gene's binarized values onto
the edge list (O(genes x edges) rows) or built a dense cells x cells matrix,
and it tested a 2x2 edge table with Fisher's exact test, treating edges as
independent. Moran's I and Geary's C looped over genes and shuffled for
significance.

Every one of these statistics, and its exact randomisation moments, is a
function of a few per-feature sums over cells and one sum over edges,
`sum_e w x_i x_j`. A null conditioned on a cell grouping needs the same sums
per group plus small group-to-group edge tables.

binSpect also needs a per-feature threshold (a 1-D 2-means split or a rank
cut), which depends on each feature's whole value distribution. Accumulators
over cell windows cannot produce an order statistic.

## Decision

Backends produce one contract, `analyzeData(x, edgeStatsParam(), network =,
groups =)`, defined in GiottoClass; the statistics are computed from it once,
in Giotto. GiottoDisk implements the contract from its accumulators (one pass
grouped by degree class, a weighted `by_cell` pass for neighbour-group sums)
and one arrow plan for the edge sum (`.pe_edge_products`).

Per-feature thresholds are solved on a per-feature value histogram, which is
additive over cells (`.pe_value_hist`: the `max` accumulator, a binning
`multiply`, `nnz` grouped `by_value`). The solver, Giotto's `histThreshold()`,
takes exact or binned histograms, so in-memory and streamed backends make the
same call. The threshold is applied in bin units (the same `multiply`, then a
per-row `compare`), so it agrees with the binning entry for entry.

## Consequences

- A new backend gets Moran's I, Geary's C and binSpect, with and without the
  covariate null, by implementing `edgeStatsParam` and
  `thresholdParam("binarize")`. A new statistic built from the same sums needs
  no backend code.
- Streamed binarization is exact when each feature's distinct positive values
  fall in distinct bins (raw counts up to `n_bins`); on continuous values,
  values sharing a bin act as ties and a split or rank cut can move by one bin.
- The view's pairs are read once into R (about 16 B a pair) and every
  edge-only quantity is derived from them. `.pe_edge_products` joins one cell
  window at a time against the triplets of that window's neighbours, so its
  peak follows the window budget plus the neighbours' entries. Windows are row
  ranges; when the stored cell order is not local (nearly every pair within a
  window), the view is first rewritten out of core in breadth-first network
  order into a temporary store, so the neighbours stay a thin rim. Measured
  source orders differ: CosMx (FOV order) is local, Xenium partly, Visium HD
  barcode order not. The rewrite costs about three passes and temporary disk
  the size of the view's values; a reader that writes cells in spatial order
  avoids it. With GiottoKernels installed, the edge sums run in its
  `edge_products_stream()` kernel instead: every cell of a block of features
  in memory with random access, so neither windows nor rewrite are needed and
  memory is bounded by the feature block.
- The covariate join-count null assumes an unweighted network.

## Alternatives considered

- **Port binSpect's edge join to the store** — O(genes x edges) intermediate
  rows and still the edge-independence Fisher test.
- **Exact quantiles by sorting each feature** — needs a global order along the
  feature axis, which a cell-major store cannot window.
- **Permutation testing on the store** — one pass per permutation; the closed
  form gives the exact moments of the same null at the cost of one pass.

## References

- Cliff, A. D. and Ord, J. K. (1981). Spatial Processes: Models and
  Applications. Pion.
- GiottoClass `R/analyze-edgeStats.R`; Giotto `R/spatautocor-stats.R`,
  `histThreshold()` in `R/normalize.R`.
