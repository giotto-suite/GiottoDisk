# 0009. Detection thresholds gate counts, not magnitudes

- **Status:** Accepted
- **Date:** 2026-08-07
- **Supersedes:** —
- **Superseded by:** —

## Context

The per-axis statistic verbs take a threshold — `detection_threshold` on
`featStatsParam` / `cellStatsParam`, `expression_threshold` on the COV HVF
methods, and Giotto's filter verbs. What it is allowed to affect had drifted
apart across GiottoDisk's implementations.

`.stream_feat_qc_stats` and `.stream_cell_qc_stats` applied
`dplyr::filter(value > thr)` *before* the aggregate, so the surviving rows fed
both the count and the sums. `total_expr` and `mean_expr` therefore shrank as
the threshold rose. `.stream_gene_stats` did the opposite: it summed
unconditionally and thresholded only the count.

Giotto is unambiguous, on the code path HVF actually uses
(`.calc_expr_general_stats`):

```r
nr_cells   = rowSums_flex(expr_values > expression_threshold)   # thresholded
total_expr = rowSums_flex(expr_values)                          # not
mean_expr  = rowMeans_flex(expr_values)                         # not
sd         = .rowSds_flex(expr_values)                          # not
```

`cov := sd / mean_expr` — the actual HVF ranking quantity — is built from two
unthresholded statistics. `featStatsParam` matches. Two independent
corroborations: the BPCells fast path just above is gated on
`expression_threshold == 0` *only* because `matrix_stats`' `nonzero` cannot take
an arbitrary threshold — if the threshold touched mean or variance, substituting
`rs["mean", ]` would be invalid at any threshold, not just non-zero ones. And
the one threshold-conditioned magnitude Giotto does compute, `mean_expr_det`,
lives in its own `_det`-suffixed column via a separate `.mean_expr_det_test()`
call rather than being folded into `total_expr`.

Surveying every use across Giotto and GiottoClass, the threshold feeds a count,
a fraction, or a keep/drop decision. The single exception is `mean_expr_det`,
and even there nothing is clipped: the BPCells identity
`rowSums(x) - rowSums(min_scalar(x, t)) + count * t` reconstructs the *full*
value of every passing entry.

## Decision

A threshold is a **detection predicate**. It selects which entries are counted
as detected; it never transforms a value that participates.

So `sum` and `sumsq` are unconditional, and only `nnz` and `sum_det` see the
threshold. A statistic conditioned on detection is legitimate — it is a mean
over a different population, not a modified total — but it goes in its own
column, following `mean_expr_det`.

`.stream_expr_accum()` implements this with separate accumulators rather than a
pre-filter, so one pass yields both the unconditional totals and the detected
subset:

```
sum      sum(value)                      every stored entry
sumsq    sum(value^2)                    every stored entry
nnz      count(value > thr)              detection count
sum_det  sum(value where value > thr)    detected entries only
```

`inclusive = TRUE` switches the comparison to `>=`, which is what *filtering*
means by a threshold — Giotto's `expression_threshold = 1` is "expressed if the
count is at least 1", unlike the statistic verbs' strict `>`.

## Consequences

The QC verbs' `total_expr` and `mean_expr` no longer shrink with the threshold.
This was latent at the default: `thr = 0` removes nothing, since stored entries
are non-zero and stay positive under both norm ops. It appeared at any non-zero
threshold, and on a store holding explicit zeros or negatives.

`mean_expr_det` was already correct and stays so — filtered sum over filtered
count is exactly the detected mean. The fix was to stop *reusing* that filtered
sum for the totals, not to remove the filter.

The cost is one extra accumulator: `sum_det` measured at ~21 ms on 19.2M
nonzeros, against a ~130 ms pass. Cheaper than the pre-filter shape it replaced,
which measured slower (0.135 s vs 0.124 s) as well as wrong.

The constraint on future statistics: if one is conditioned on detection, give it
its own column and its own accumulator. Do not filter the stream and let the
conditioning leak into everything computed from it.

## Alternatives considered

- **Filter once, before the aggregate** — what QC did. One fewer accumulator and
  measurably slower, and it makes the threshold silently change quantities that
  are defined without reference to it.
- **Two passes**, one filtered and one not. Correct, and twice the scan.
- **Match Giotto's `featStatsParam` shape exactly** by calling something like
  `.mean_expr_det_test()` separately rather than carrying `sum_det` in the same
  pass. Rejected: same arithmetic, extra traversal, and the accumulator's
  `stats` selector already makes it opt-in per caller.

## References

- `R/methods-analyzeData.R` — `.stream_expr_accum()`, `.pe_featstats()`
- `Giotto/R/variable_genes.R` — `.calc_expr_general_stats()`,
  `featStatsParam` / `cellStatsParam` methods
- `Giotto/R/auxiliary_giotto.R` — `.mean_expr_det_test()`
- ADR 0004 (roles and invariants)
