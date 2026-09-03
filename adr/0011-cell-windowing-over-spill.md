# 0011. Cell windowing over spill for a large grouped join

- **Status:** Accepted
- **Date:** 2026-08-20
- **Supersedes:** —
- **Superseded by:** —

## Context

A grouped expression statistic needs a join: the per-cell group assignment has to
meet the values before the aggregate can key on it. The aggregate itself is small
(`O(features × groups)`), but the join emits one row per stored value, so on a
large store its output is large — `O(nonzeros)`.

Large joins exceed memory easily, and this one does. Measured on a 5,000 × 50,000
store at 3% density with 12 groups: 187 MB engine peak. At atlas scale (≈556M
stored values) the same plan wanted roughly 21 GB against 8 GB of RAM. Acero has
no spill path, so it failed outright rather than slowing down.

That leaves two ways to get a large join under a memory ceiling:

1. **Brute force it** — run the join whole and let the engine spill the excess to
   disk. DuckDB has the machinery for this.
2. **Do less work** — split the input, join and aggregate each part, combine the
   parts. Available only when the statistic can be combined from parts.

Here it can. `sum`, `sumsq`, `nnz` and `sum_det` are exactly additive over cells,
and `.pe_pool_moments()` already relied on that to build one-vs-rest groups
without a second scan.

Callers had meanwhile started batching the **feature** axis to get under the
ceiling — option 2, along the wrong axis. Stores are written cell-major, so a
feature-side predicate prunes no row groups and every batch rescans in full: cost
linear in batch count rather than in features per batch (20/10/5 batches at
1.77/0.57/0.31 s on the store above; one real run took 15 minutes for 76
batches).

## Decision

Take option 2 on the **cell** axis: partition the cells, aggregate each part, add
the parts. The window is derived per read from the store's shape against a
fraction of free RAM (`.recommend_chunk_size()`); a budget that covers the view
yields one part, which is the single plan that existed before. Windowing is not a
mode and no option switches it on.

Chosen over spilling because it does less work. Spilling writes the
`O(nonzeros)` join output to disk and reads it back, to produce a result that
never needed it to exist all at once; windowing never builds it. A contiguous
cell range also prunes row groups, which a whole-store scan cannot. Where both
fit the comparison is fewer operations against more, and under a tight ceiling
the spilling engine was not observed to stay bounded at all (see *Alternatives*)
— so the availability of one does not by itself reopen this.

This decides *whether to window*, not *which engine reads a window*. Those are
separate, and the second is settled elsewhere: adr/0012 makes DuckDB a native
carrier for the expression scan, and a DuckDB stat accumulator is expected to
follow. Both window. Nothing here argues against either.

Partials are folded as they arrive rather than collected and reduced at the end,
so retained state is `O(groups)` rather than `O(groups × windows)`.

Scope: windowing itself predates this decision — PCA and the `storeWrite()` bake
already read a chunk at a time, because they materialize one by nature. What is
decided here is the grouped statistic's join, which had no bound at all.

## Consequences

- Engine peak on the measured store falls 187 MB → 44 MB across 20 windows,
  and the atlas-scale pass becomes possible rather than merely slower.
- Windows are exact, not approximate — the additivity is load-bearing. A
  statistic that is **not** additive over cells cannot use this path. Anything
  needing a global ordering along the feature axis (rank, median, Wilcoxon) has
  no streaming path as a result; `scranMarkersParam` refuses `test_type =
  "wilcox"` for this reason.
- Windows cost. Roughly 45–65 ms per window of plan setup and fold on the
  measured store: 0.27 s at one window, 1.60 s at twenty. So where the whole-store
  plan *fits*, it is the faster one (0.20 s), which is exactly why a budget that
  covers the view yields a single window. Decomposition wins on work done, not on
  constant factors, and the argument is for the largest window the budget allows —
  not for a fixed window count.
- The comparison against a spilling engine **was** measured, but narrowly: one
  store shape, one DuckDB version, and a failure whose mechanism is unexplained
  (see *Alternatives*). What is measured more robustly is that windowing never
  builds the intermediate at all (187 MB → 44 MB engine peak at 20 windows) and
  that a cell range prunes row groups. Someone reversing this decision should
  re-measure rather than inherit either number.
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

- **Let Acero run the whole join** — what existed. The join output does not fit
  and Acero cannot spill, so the pass fails rather than degrading. Rejected by
  the failure.
- **Spill it, with DuckDB** — the brute-force option, and it did not behave as
  expected. Measured on this shape the whole join fit at a 320 MB limit and raised
  `Out of Memory Error: could not allocate block` at 256 MB; windowed, the same
  plan ran at 128 MB. Not a claim that DuckDB cannot spill — `temp_directory` is
  set on these connections — but at the caps tried the whole join either fit or
  died, so the window, not the engine, is what bounded the pass. Rejected on work
  done besides, which is why a window wins even where both fit: spilling writes
  the `O(nonzeros)` join output and reads it back where windowing never builds it,
  and a whole-store scan forgoes the row-group pruning a contiguous cell range
  gets. Keep spill in view for an aggregate that genuinely *cannot* be split — a
  rank, a median — and note two limits on the SQL route: `@post_ops` has no SQL
  form, and a spilling engine still needs the same memory limit the window is
  derived from. **This rejects *spilling*, not DuckDB.** A native DuckDB scan
  landed separately (adr/0012), a DuckDB stat accumulator is expected after it,
  and both window exactly as Acero does. Read this bullet as being about the
  execution strategy, not about the engine.
- **Batch the feature axis** — what callers were doing. The same split-and-combine
  idea on the axis that cannot prune: every batch rescans in full, so cost is
  linear in batch count. Strictly worse than a cell window for the same bound.
- **Window, but reduce partials at the end** — simpler, and one fewer fold per
  window. Retains `O(groups × windows)`, so tightening the window to save memory
  would cost memory. Rejected as backwards.
- **A gene-major second copy of the store** — would make feature narrowing prune
  and give rank statistics a path. Doubles storage and adds a write path; a
  storage decision rather than an execution one, and not needed for any statistic
  that decomposes. Left open.

## References

- `vignettes/expression_windows.Rmd` (`vignette("expression_windows")`) — the
  long-form argument, with the measurements and the worked sizing.
- `R/utils-pestore-ops.R` — `.pe_windows()`, `.pe_chunk_ranges()`,
  `.pe_window_store()`: the seam every windowed pass attaches to.
- `R/methods-analyzeData.R` — `.pe_accum_raw()` and the two windowed
  accumulators; `.pe_fold_partial()` carries the reassociation note.
- `R/chunk-sizing.R` — `.recommend_chunk_size()` and `storeChunkInfo()`, which
  carries the two options.
- adr/0003 (payloads keyed by on-disk id — why a window needs no payload
  slicing), adr/0008 (axis predicate shapes — why a cell range prunes and a
  feature predicate does not).
