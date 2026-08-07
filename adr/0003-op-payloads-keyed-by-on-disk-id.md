# 0003. Op payloads are keyed by on-disk id, so `[` needs no subset slicing

- **Status:** Accepted
- **Date:** 2026-07-22
- **Supersedes:** —
- **Superseded by:** —

## Context

A `parquetExprStore` view narrows by `[`, which only sets `@cell_idx` /
`@gene_idx` against the same file — the Parquet payload is untouched. Ops queued
before the subset carry payloads that are themselves indexed vectors: a
`multiply`'s factor vectors, one entry per cell or per feature.

The question is what those vectors are indexed *by*. Two candidates, and they
behave differently under `[`:

- **View position** — invalidated by every subset, so each `[` must filter every
  op payload down to the surviving positions and renumber.
- **On-disk id** (`row_id` / `col_id`) — unchanged by a subset, because narrowing
  a view does not move ids in the file.

The first was built. When the normalization op carried a `(source_id,
orig_row_id, scalef)` data.table, `.pe_op_table_keys` was a registry naming which
payload tables were view-position-keyed, with a dispatcher that filtered them to
surviving keys; eight call sites in `[` fed it. By the time the phase-split
landed, the norm payload had become a plain vector indexed by on-disk id and
**every one of those eight sites had become a no-op**.

The `@stats` marginals (`row_nnz`, `col_nnz`) had already made the same choice
for the same reason, and their docs state it outright: keyed by on-disk id, not
by name and not by view position, "which is what makes them invariant under
`[`". Names break under feature renaming; view positions are invalidated by
every subset.

## Decision

Op payloads are keyed by on-disk id. `[` carries `@ops` and `@post_ops` through
untouched — there is no subset-slice step, deliberately. `.pe_op_table_keys` and
its dispatcher were removed in `fa3ee64`.

Where a consumer needs the on-disk id for a view position, it converts at the
point of use (`.pe_orig_col`), rather than the store maintaining rewritten
payloads.

## Consequences

- `[` stays cheap and total: no op payload is re-derived, so subsetting cost does
  not grow with chain length.
- The invariant is load-bearing and unenforced. A future op that wants a payload
  keyed by view position would silently read the wrong entries after a subset —
  nothing checks the keying. That is the trap this record exists to name.
- Two related invariants must hold alongside it:
  - `@stats` vector lengths are the *file's* dimensions, not the view's
    (`@n_cells` / `@n_genes` describe the view).
  - `storeWrite()` **renumbers** ids when the input is subset, so a written
    store's payloads and marginals are computed against the new file rather than
    inherited. On-disk-id keying is invariant under `[`, *not* across a write.
- Revisit only if an op genuinely needs a table keyed by view positions. The
  slicing step goes back where the removed registry was
  (`R/utils-pestore-ops.R`, *subset slice*) — but prefer redesigning the payload
  onto on-disk ids, which needs none.

## Alternatives considered

- **Keep the registry + dispatcher** (`.pe_op_table_keys`) — the mechanism was
  correct for view-position-keyed payloads, but after the payload redesign it
  filtered nothing at eight call sites. Rejected as machinery for a population
  of zero.
- **Key payloads by identifier name** (cell barcode / feature id) — survives
  subsetting, but breaks under feature renaming, and costs a string join per op
  application on the hot path.
- **Re-derive payloads at read time from the current view** — no stale-key risk,
  but it makes every read recompute what the op already knows, and the scale
  factors are not always recoverable from the narrowed data.

## References

- `fa3ee64` refactor: phase-split @ops (arrow) + @post_ops (R) execution
  (removal of `.pe_op_table_keys`)
- `18825e2` feat: @ops on union store + (source_id, row_id) composite key +
  subset-slice (the introduction it replaced)
- `R/utils-pestore-ops.R` — *subset slice*
- `R/class-parquetExprStore.R` — `@stats` slot docs, on-disk-id keying rationale
- `R/methods-parquetExprStore.R` — `.pe_orig_col`, `.pestore_write_baked`
- [0002](0002-post-ops-phase-split.md) — the phase split this landed with
