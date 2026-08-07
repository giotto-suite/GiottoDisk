# 0008. Arrow axis predicates never add a range "for good measure"

- **Status:** Accepted
- **Date:** 2026-06-12
- **Supersedes:** —
- **Superseded by:** —

## Context

`[` on a `parquetExprStore` records the kept positions in `@cell_idx` /
`@gene_idx`; `storeRead` turns those into Arrow predicates on `row_id` /
`col_id`. How the predicate is *shaped* decides whether Parquet row groups can
be skipped.

The naive form, `x %in% <large integer vector>`, never prunes and pays a hash
probe per returned row. A range, `x >= lo & x <= hi`, is checked against each
row group's min/max statistics, so untouched groups are never opened.

But a range only prunes on an axis whose layout it matches. The file is sorted
`(row_id, col_id)`, so a row group holds a contiguous span of cells and, within
each, the full spread of features. A `row_id` range therefore prunes; a `col_id`
range covers nearly every group and prunes nothing, while still costing a
comparison per row.

Measured on Atera:

- gene axis, 2k of 18k genes, HVG-ranked: stacking a range onto `is_in(kept)`
  cost **+6.4% for zero pruning**
- cell axis, contiguous chunk: a range was **3.9x faster** than the `%in%` form

## Decision

`.pe_axis_pred()` classifies the kept set and `.pe_axis_pred_exprs()` emits one
of three shapes:

1. **gapless** (unique count equals range span) — range alone. It matches the
   kept set exactly, so no membership test is needed at all.
2. **gaps, few dropped** — range **and** `!(x %in% dropped)`. Here the range is
   required for *correctness*, not speed: without it, rows outside `[lo, hi]`
   pass the negated test.
3. **gaps, few kept** — `x %in% kept` alone. Already correct and complete, and
   the range is deliberately **not** emitted.

All three are exact. The predicate admits precisely the in-view entries, so
nothing downstream has to re-filter.

Bounds come from `min`/`max`, never first/last: `idx` is not guaranteed sorted
(`feats_to_use` may be HVG-rank ordered). Gap detection runs on unique values so
duplicates cannot make `n == span` accidentally true, and `dropped` is
materialised only when case 2 wins.

## Consequences

Case 3 is the one that looks wrong and is not. A range there is free correctness
and cheap-looking, so the natural instinct is to add it for symmetry — that is
the +6.4% above, paid on every scan for nothing. This ADR exists mostly to
answer that instinct.

Cell subsets prune and gene subsets do not, which is a property of the sort
order rather than of this code. It is why the PCA path materialises a transient
HVG-narrowed store past a feature-ratio threshold
(`giottodisk.pca_bake_max_ratio`) instead of relying on the predicate.

Because the predicates are exact, consumers can treat what comes back as the
view. The `!is.na(match(...))` and inner-join guards downstream are backstops
against a store whose on-disk ids fall outside its declared index, not filters —
measured as dropping zero rows across no-subset, scattered-gene, gapped-cell,
and both-axis cases.

Revisit if the on-disk sort order changes. A feature-major or hybrid layout
would flip which axis prunes, and case 3's reasoning with it.

## Alternatives considered

- **Always `x %in% kept`.** Simplest and always correct. Rejected on the cell
  axis measurement — 3.9x — where a contiguous chunk is the common access
  pattern for every streaming reader.
- **Always emit the range alongside the membership test.** Symmetric and easy to
  reason about. Rejected: +6.4% on the gene axis for zero pruning, and the gene
  axis is where large subsets actually occur.
- **Sort the parquet by `(col_id, row_id)` instead**, making gene subsets the
  prunable axis. Rejected: cell-chunked streaming is the dominant access pattern
  and the CSC layout of the materialised `dgCMatrix` matches cell-major, so a
  chunk read needs no reordering.
- **Pre-bake every subset into a narrowed store.** What the PCA path does past a
  ratio threshold, and too expensive as a general rule — it rewrites the file
  for what is otherwise a free view.

## References

- `R/methods-parquetExprStore.R` — `.pe_axis_pred()`, `.pe_axis_pred_exprs()`
- ADR 0006 (view state is not chain state) — why a subset stays a view
- Measurements: Atera FFPE, 2k of 18k HVG-ranked genes; contiguous cell chunk
