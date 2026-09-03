# 0013. No second stat accumulator engine: the win is per-window overhead, and it is small

- **Status:** Accepted
- **Date:** 2026-09-03
- **Supersedes:** —
- **Superseded by:** —

## Context

Once `storeRead(output = "duckdb")` gave `parquetExprStore` a native DuckDB scan
(adr/0012), a DuckDB *stat accumulator* looked like the obvious next step: the
cell-window seam hands its consumer a store rather than a query, so which engine
reads each window is the consumer's to choose, and a third consumer alongside
`.pe_accum_acero_windowed()` would have cost no new carrier.

A prototype benchmark supported it. On 5,000 x 50,000 at 3% density, 12 groups,
20 windows, four shapes were compared, every arm elementwise identical:

| shape | vs Acero |
|---|---|
| Acero, shipped accumulator | 1.00x |
| `storeRead(output = "duckdb")` per window | ~0.96x |
| one hoisted `storeRead`, per-window dplyr predicate | ~0.97x |
| render the query once, `dbBind` the window bounds per pass | **~5.5x** |

`roadmap.Rmd` carried that 5.5x and a design to match.

## Decision

Do not build it. Both engines are the same speed on this work, and the shape
that looked 5.5x faster was measuring a configuration real window sizing never
produces.

Timing the shipped Acero accumulator across window counts on one store fits
exactly:

```
total = 0.026 s of work + 37 ms x n_windows        (R^2 = 1.000)
```

| windows | nonzeros/window | total | overhead share |
|---:|---:|---:|---:|
| 1 | 4.0M | 0.060 s | 62% |
| 5 | 800k | 0.216 s | 85% |
| 20 | 200k | 0.766 s | **96%** |
| 40 | 100k | 1.500 s | 98% |

The whole scan-and-aggregate over 4M nonzeros costs **26 ms**; everything else
is a fixed per-window cost of 37 ms. The prototype pinned 20 windows onto a
store the budget reads in **one**, so its pass was 96% fixed overhead — and the
5.5x is the fraction of that overhead a prepared statement removes. It says
nothing about how fast either engine processes data, which is what the two
parity arms already showed.

The overhead share does not improve with job size, because it is not a function
of job size:

```
overhead share = 37 ms / (time to process ONE budget-sized window)
```

Both terms are independent of the total. A larger store means *more* windows,
not slower ones. At a ~1 GB budget a window holds ~60M nonzeros, which is ~0.4 s
of work at the throughput measured here, so the share sits near 9% and the
ceiling near **1.10x** — and that is the optimistic end, since the fixture was
page-resident where a real store is I/O-bound and the share falls further.

The absolute prize is bounded the same way. It is `n_windows x ~37 ms`; a
10M-cell atlas runs on the order of 540 windows, so the whole saving is about
**20 seconds** on a job measured in tens of minutes.

The ratio is therefore largest exactly where the stakes are lowest. There is no
configuration in which both the multiplier and the absolute saving are worth
having, because the eliminated term is fixed per window and windows are sized to
the budget.

## Consequences

- No `giottodisk.stat_accum_engine` option, no `.pe_accum_duckdb_windowed()`, no
  engine fallback path, no SQL aggregate builder, and no union `source_id`
  special case. `.pe_accum_acero_windowed()` remains the only windowed
  accumulator.
- `storeRead(output = "duckdb")` is unaffected. Nothing here argues against the
  carrier — adr/0012 stands, and a caller who wants a `tbl_dbi` still gets one.
  What is rejected is a *second accumulator* selected by engine.
- The optimization target moves to **window count**, and specifically to the
  byte accounting that sets it. `n_windows x 37 ms` is the only term either
  engine can still reduce, and windows are currently smaller than the budget
  allows because the sizing charges every path for what the most expensive path
  holds. Three reservations are wrong for the Acero accumulator:
  `bytes_per_nz = 48` (assumes a materialized triplet frame; it streams and
  retains only the aggregate), and a `(k + p) x 8` per-cell sketch it never
  allocates — `p = 10` is hardcoded in `.chunk_budget()` so even `k = 0L` is
  charged 80 B/cell, and `.exprbase_chunk_size()` passes `k = 50L` so the
  `storeWrite` bake is charged 480 B/cell.

  The two sketch terms are the important ones, because they scale with
  `n_cells` and so eat a growing share of a fixed budget: at 8 GB free, the
  accumulator's 80 B/cell reservation is 0.04 GB at 500k cells but 4.0 GB at
  50M, which exceeds the whole 25% budget and drops sizing onto the 5%
  emergency floor. That makes the bound stop tracking the model, not merely
  makes it slow. Unlike the per-window setup cost, which stays a fixed fraction
  of runtime, this error compounds with job size — which is why the
  window-count saving grows from ~1% on a 500k-cell store to ~33% at 50M cells
  on constrained RAM. See `vignettes/articles/roadmap.Rmd`.
- adr/0011 says a DuckDB stat accumulator "is expected" after the native scan.
  That expectation is withdrawn here. 0011's own decision — window rather than
  spill — is untouched and does not depend on it.
- Hoisting is worth naming separately, because it is the half that *does* have
  an Arrow analogue and it is still not worth much. Acero re-opens the Dataset
  once per window (measured: 12 windows, 12 `open_dataset` calls), and one
  `open_dataset` is 2.0 ms against 39 ms of per-window cost — a ~5% ceiling.
  The other half, a prepared parameterized plan, has no Acero equivalent at all:
  an `arrow_dplyr_query` cannot be bound, so every window rebuilds and
  re-optimizes. That asymmetry is structural rather than an oversight.
- **Revisit if** per-window overhead stops being fixed — an op whose setup cost
  scales with window payload would change the arithmetic — or if a verb appears
  whose windows are necessarily small for a reason other than the budget.
  Re-derive the fit before reopening; do not inherit the 37 ms.

## Alternatives considered

- **Build it anyway, behind an opt-in option.** The roadmap design: no default,
  unset resolves to `"arrow"` with a one-shot nudge. Rejected on the arithmetic
  above — a second accumulator, a fallback for non-gapless windows, hand-written
  aggregate SQL and a union case, for ~20 s on a 30-minute job. The indirection
  outlasts the benefit.
- **Hoist the Dataset handle on the Acero path.** The Arrow half of the same
  idea, and genuinely available. Rejected as not worth it *here*: ~5% ceiling,
  and caching a handle across calls means holding engine state and reasoning
  about staleness when files move, which the store model has deliberately kept
  out. Left open as an independent question rather than as an engine decision.
- **Re-measure on a larger store before deciding.** The throughput used for the
  projection (154M nonzeros/s) comes from a page-resident 4M-nonzero fixture.
  Rejected as unnecessary to *this* decision: I/O-bound reads make per-window
  work larger and the overhead share smaller, so the conclusion moves further in
  the direction already taken. Worth doing if the number is ever load-bearing
  for something else.

## References

- `R/methods-analyzeData.R` — `.pe_accum_acero_windowed()`, the accumulator that
  stays; `.pe_accum_chunk_size()`, which sets the window and so the 37 ms count.
- `R/utils-pestore-ops.R` — `.pe_windows()`, the seam a second consumer would
  have attached to.
- adr/0011 (cell windowing over spill) — the decision this leaves intact, and
  the source of the withdrawn expectation.
- adr/0012 (one predicate classifier, many carriers) — the native DuckDB scan,
  unaffected.
