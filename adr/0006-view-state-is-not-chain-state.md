# 0006. View state is not chain state; window-dependent ops bake at push time

- **Status:** Accepted
- **Date:** 2026-08-07
- **Supersedes:** —
- **Superseded by:** —

## Context

A `parquetExprStore` carries two independent kinds of pending state:

- **view state** — `@cell_idx` / `@gene_idx`, the window onto the file. Set by
  `[`, which rewrites nothing on disk.
- **chain state** — `@ops` / `@post_ops`, the sequence of steps to apply to
  values.

Some ops are pure functions of a value and care about neither: `log` is
`log1p(v)/log(base)` whatever window is open. Others are defined by a
*population statistic* over the window that was current when the op was
requested. `multiply(axis = "cell")` as produced by `libraryNormParam` is the
case that matters — its factors are `scalefactor / colSums`, and which columns
are summed depends on which features are in view.

That makes the two kinds of state interact, and the interaction has to be
decided rather than left to fall out. Concretely: after
`processData(pe, libraryNormParam(...))`, what should `pe[hvg_idx, ]` mean? The
factors were computed over all features. Do they still apply?

The question is not hypothetical. HVG selection after normalization is the
standard workflow, so the subset-after-normalize path is the common one, not an
edge case.

## Decision

**View state and chain state are independent.** `[` narrows the window and
never touches the chain; the chain never consults the window at read time.

**An op whose semantics depend on the window must bake that dependence into
frozen data at push time.** `libraryNormParam` runs its aggregate when the verb
is called and freezes the result into the record as a payload. From then on the
record is a pure function of a value and an axis id — no different from `log` in
what it needs at read time.

An op that *cannot* bake — one that would have to re-derive its statistic
against whatever window is open at read time — does not belong on the chain.
The escape is materialization: `storeWrite` bakes the chain into on-disk values
and returns a store with empty chains, after which the statistic can be taken
over the new population by running the verb again.

The corollary, which is the user-visible part: **re-running the producer is how
you ask for a statistic over a new population.** Subsetting is not.

## Consequences

Normalize-then-subset keeps full-library scaling. A cell's normalized value is
the same regardless of which features are in view, which is what makes values
comparable across analyses — and it is the behaviour every other toolbox has,
since normalization corrects for sequencing depth, a property of the cell rather
than of the gene set someone is looking at. Verified against dense references
in both orders; the two orders are different operations and both are reachable.

Subsetting therefore never invalidates a payload, which is what lets payloads be
keyed by on-disk id (0003) and lets `[` leave the chain alone entirely.

Producers pay an eager cost. `processData(pe, libraryNormParam(...))` is a full
pass over the store, not a queued recipe — it reads like the latter at the call
site and nothing signals otherwise. The payload it freezes is O(axis) and rides
on the store thereafter, through `saveRDS` and every subset.

Unions must reject ops-dirty substores, and do. A per-substore chain would carry
factors tuned to that substore's population, which is not the union's; there is
no sound way to compose them. The canonical order is `cbind` raw, then
`processData` on the union.

The constraint on future ops: if a new op type's meaning depends on the window,
its producer must do the work up front and store the answer. If that is not
possible, the op does not go on the chain — and if a consumer needs the
statistic recomputed for a narrowed view, the answer is a write and a re-run,
not a smarter chain.

Revisit if an op appears whose population statistic is cheap enough to recompute
per read *and* whose users expect it to track the window. Nothing in the current
inventory is like that, and the comparability argument above suggests
normalization never will be.

## Alternatives considered

- **Re-evaluate window-dependent ops at read time.** Rejected on semantics
  before performance: the same cell would normalize differently depending on
  which genes are in view, so no two analyses would be comparable. It also turns
  every read into an aggregate pass over the window.
- **Invalidate the op on subset — error, or silently drop it.** Rejected: the
  subset-after-normalize path is the standard workflow, so this errors on the
  common case.
- **Slice the payload to the surviving window on `[`.** This was built, for the
  cell axis, when the payload was a `(source_id, orig_row_id, scalef)` table. It
  is compatible with freezing — slicing a frozen payload does not recompute it —
  but it is unnecessary once payloads are keyed by on-disk id, and it was the
  machinery 0003 removed.
- **Forbid window-dependent ops entirely, requiring a `storeWrite` between
  normalize and anything else.** Rejected: library normalization is the central
  use case, and this would make the ordinary pipeline write an intermediate
  store.

## References

- ADR 0003 (payloads keyed by on-disk id), ADR 0005 (`@ops` as prefix)
- `R/class-parquetExprStore.R` — `unionParquetExprStore()` constructor, which
  enforces this contract for substores
- `R/stream-normalize.R` — `processData(·, libraryNormParam)`, the baking
  producer
- `tests/testthat/test-stream-normalize.R` — gene subset preserves the factor
  payload; cell subset likewise
