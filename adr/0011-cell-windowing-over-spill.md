# 0011. Bounded expression passes window the cell axis rather than spill

- **Status:** Accepted
- **Date:** 2026-08-20
- **Supersedes:** —
- **Superseded by:** —

## Context

A grouped expression statistic joins a per-cell group assignment onto the values
before aggregating. The aggregate is `O(features × groups)`; the join in front of
it emits `O(nonzeros)`. Acero has no spill path, so the whole-store form of that
plan does not degrade — it fails. Measured on a 5,000 × 50,000 store at 3%
density, 12 groups: 187 MB engine peak for the whole-store pass. At atlas scale
(≈556M stored values) the same plan wanted roughly 21 GB against 8 GB of RAM.

Callers were working around it by batching the **feature** axis. That is the
expensive workaround: stores are written cell-major, so a feature-side predicate
prunes no row groups and every batch rescans in full. Cost is linear in batch
count, not in features per batch — 20/10/5 batches measured at 1.77/0.57/0.31 s
on the store above. One user run took 15 minutes for 76 batches.

The statistics involved (`sum`, `sumsq`, `nnz`, `sum_det`) are exactly additive
over cells. `.pe_pool_moments()` already relied on that to build one-vs-rest
groups without a second scan.

## Decision

Bounded passes window the **cell** axis; the statistic accumulators combine the
parts by addition. The window is derived per read from the store's shape against
a fraction of free RAM (`.recommend_chunk_size()`); a budget that covers the view
yields one window, which is the single plan that existed before. Windowing is not
a mode and no option switches it on.

Scope: windowing itself predates this decision — PCA and the `storeWrite()` bake
already read a chunk at a time, because they materialize one by nature. What is
decided here is that the **grouped statistic** windows too, rather than being
handed to an engine that would spill. It is the only pass whose alternative was
failure rather than a different chunk size.

Partials are folded as they arrive rather than collected and reduced at the end,
so retained state is `O(groups)` rather than `O(groups × windows)`.

## Consequences

- Engine peak on the measured store falls 187 MB → 44 MB across 20 windows,
  and the atlas-scale pass becomes possible rather than merely slower.
- Windows are exact, not approximate — the additivity is load-bearing. A
  statistic that is **not** additive over cells cannot use this path. Anything
  needing a global ordering along the feature axis (rank, median, Wilcoxon) has
  no streaming path as a result; `scranMarkersParam` refuses `test_type =
  "wilcox"` for this reason.
- Windows cost. Roughly 45–65 ms per window of plan setup and fold on the
  measured store: 0.27 s at one window, 1.60 s at twenty. This is an argument for
  taking the largest window the budget allows — which is what the sizing returns
  — not for a fixed window count.
- Reassociating a float sum is not bitwise invariant. Counts stay exact; float
  accumulators move by 2–3 ULP across window counts. Results are reproducible to
  tolerance, not bitwise, and tests must be written that way.
- The bound now depends on the sizing arithmetic rather than on an engine
  guarantee. `bytes_per_nz` is calibrated per read shape, and the `10000`-cell
  floor is a floor, not a guarantee — on very wide data it can exceed the budget.
- **Revisit if** a verb needs a non-decomposable statistic (the rank case), or
  if the op chain grows a step that lands on `@post_ops` in the normal pipeline —
  that makes the ungrouped accumulator windowed too, which puts a discrete
  feature-selection cut downstream of a reassociated float sum.

## Alternatives considered

- **Let Acero run the whole-store plan** — what existed. No spill, so it fails
  rather than degrading. Rejected by the failure itself.
- **DuckDB with spill enabled** — survives, but only by materializing the
  intermediate, writing it out and reading it back. The statistic decomposes, so
  nothing needs to exist at once; spill pays to store work that need not be
  produced, and forgoes the row-group pruning a cell window gets. Also:
  `@post_ops` has no SQL form, a spilling engine still needs the same memory
  limit the window is derived from, and `parquetExprStore` extends
  `queryableStore` rather than `parquetStore`, so `output = "duckdb"` is
  Arrow-backed today — the native scanner is not wired to this class and wiring
  it needs the axis predicates rendered as SQL.
- **Batch the feature axis** — what callers were doing. Prunes nothing, cost
  linear in batch count. Strictly worse than a cell window for the same bound.
- **Window, but reduce partials at the end** — simpler, and one fewer fold per
  window. Retains `O(groups × windows)`, so tightening the window to save memory
  would cost memory. Rejected as backwards.
- **A gene-major second copy of the store** — would make feature narrowing prune
  and give rank statistics a path. Doubles storage and adds a write path; a
  storage decision rather than an execution one, and not needed for any statistic
  that decomposes. Left open.

## References

- `vignettes/chunking.Rmd` (`vignette("chunking")`) — the long-form argument,
  with the measurements and the worked sizing.
- `R/utils-pestore-ops.R` — `.pe_windows()`, `.pe_chunk_ranges()`,
  `.pe_window_store()`: the seam every windowed pass attaches to.
- `R/methods-analyzeData.R` — `.pe_accum_raw()` and the two windowed
  accumulators; `.pe_fold_partial()` carries the reassociation note.
- `R/chunk-sizing.R` — `.recommend_chunk_size()` and `storeChunkInfo()`, which
  carries the two options.
- adr/0003 (payloads keyed by on-disk id — why a window needs no payload
  slicing), adr/0008 (axis predicate shapes — why a cell range prunes and a
  feature predicate does not).
