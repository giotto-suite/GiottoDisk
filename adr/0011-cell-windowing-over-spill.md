# 0011. A decomposable statistic is computed by decomposition, not by a join over the whole store

- **Status:** Accepted
- **Date:** 2026-08-20
- **Supersedes:** —
- **Superseded by:** —

## Context

The grouped expression statistics (`sum`, `sumsq`, `nnz`, `sum_det`) are exactly
additive over cells: the value for a set of cells is the sum of the values for
any partition of it. `.pe_pool_moments()` already relied on this to build
one-vs-rest groups without a second scan.

The plan they were being computed with does not use that. Joining a per-cell
group assignment onto the values and then aggregating produces `O(nonzeros)`
intermediate rows in order to return `O(features × groups)` — the intermediate
exists because of the plan, not because of the problem. Nothing needs to see all
those rows at once, and nothing downstream reads them.

Acero's lack of a spill path is what made that visible rather than merely
wasteful: the whole-store form does not degrade, it fails. Measured on a
5,000 × 50,000 store at 3% density with 12 groups, 187 MB engine peak; at atlas
scale (≈556M stored values) it wanted roughly 21 GB against 8 GB of RAM. The
failure is the symptom. The mismatch between an additive statistic and a
materializing plan is the cause.

Callers had started batching the **feature** axis to get under the ceiling. That
is a worse plan again: stores are written cell-major, so a feature-side predicate
prunes no row groups and every batch rescans in full — cost linear in batch count
rather than in features per batch (20/10/5 batches at 1.77/0.57/0.31 s on the
store above; one real run took 15 minutes for 76 batches).

## Decision

Compute these statistics by decomposition: partition the **cell** axis, aggregate
each part, add the parts. The window is derived per read from the store's shape
against a fraction of free RAM (`.recommend_chunk_size()`); a budget that covers
the view yields one part, which is the single plan that existed before. Windowing
is not a mode and no option switches it on.

This is chosen as the *correct shape* for the computation, not as a way around an
engine limitation. Partitioning does strictly less work than materializing: no
`O(nonzeros)` frame is built, nothing is written out and read back, and a
contiguous cell range prunes row groups a whole-store scan cannot. A spilling
engine would make the materializing plan survivable; it would not make it right,
which is why the arrival of one is not a reason to revisit this.

Partials are folded as they arrive rather than collected and reduced at the end,
so retained state is `O(groups)` rather than `O(groups × windows)`.

Scope: windowing itself predates this decision — PCA and the `storeWrite()` bake
already read a chunk at a time, because they materialize one by nature. What is
decided here is that the grouped statistic is *computed differently*, not merely
read in smaller pieces.

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
- The claim that this beats spilling is **reasoned, not measured**: no benchmark
  against a spilling engine was run. What is measured is that it never builds the
  intermediate at all (187 MB → 44 MB engine peak at 20 windows) and that a cell
  range prunes row groups. Someone reversing this decision should measure that
  comparison rather than inherit the assumption.
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

- **Let Acero run the whole-store plan** — what existed. Rejected because it is
  the wrong plan for an additive statistic, not because Acero lacks spill; the
  missing spill only decided whether the waste showed up as a crash or as
  latency.
- **DuckDB with spill enabled** — the same wrong plan, made survivable. Spilling
  materializes the `O(nonzeros)` intermediate, writes it out and reads it back, to
  produce a result that never required it to exist; it also forgoes the row-group
  pruning a contiguous cell range gets. Spill is the right tool for an aggregate
  that genuinely cannot be decomposed — a rank, a median — and these are not
  those. Three further practical limits: `@post_ops` has no SQL form, a spilling engine still needs the same memory
  limit the window is derived from, and `parquetExprStore` extends
  `queryableStore` rather than `parquetStore`, so `output = "duckdb"` is
  Arrow-backed today — the native scanner is not wired to this class and wiring
  it needs the axis predicates rendered as SQL.
- **Batch the feature axis** — what callers were doing. Also a decomposition, but
  along the axis that cannot prune: every batch rescans in full, so cost is linear
  in batch count. Strictly worse than a cell window for the same bound.
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
