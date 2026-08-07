# 0001. No positional row indexing on parquet stores

- **Status:** Accepted
- **Date:** 2026-04-17 (earliest written record; the decision predates the docs commit)
- **Supersedes:** —
- **Superseded by:** —

## Context

`parquetStore` and friends stand in for in-memory tables in Giotto workflows,
and the in-memory idiom is `x[i, ]` — positional row access. Supporting it would
make stores drop-in wherever a `data.table` is used today.

Two things stand against it:

- **Parquet is not row-addressable.** Reaching row *i* means scanning row groups
  and counting, or maintaining a separate offset structure. There is no seek.
- **The index does not fit the scale target.** Point/transcript stores are
  designed for 2 trillion+ rows (design.Rmd, *Scale Targets*); a positional index
  over 2 billion rows already costs roughly 16GB on disk, and that is three
  orders of magnitude below the target.

The access patterns that actually motivated the request turned out to be
value-based (filter to these ids) or statistical (give me a subsample), neither
of which needs positions.

## Decision

Positional row indexing is not supported. Row access is value-based via
`subset()` / `[i, on = ]`, or statistical via `rowSample()`.

`row_index` exists as a special column, but it is intrinsic row *ordering within
a source* — a stable join key, not an addressable offset. `(source_id,
row_index)` is the composite identity used across rbind/union stores and by the
`spat_relate` narrow path.

## Consequences

- Stores are not fully substitutable for `data.frame`/`data.table`; any consumer
  written against `x[i, ]` must be rewritten before it can take a store. This is
  a recurring cost at every new integration seam, paid deliberately.
- Ops stay expressible as a single Arrow scan, since nothing in the chain needs
  to resolve a position.
- `nrow()` must query `COUNT(*)` rather than read an index — see the separate
  invariant that it returns a double (AGENTS.md).
- Revisit only if parquet gains addressable seeks, or if a bounded-size store
  class appears where an index is affordable and the class can be kept distinct
  from the transcript-scale ones.

## Alternatives considered

- **Materialize a positional index at write time** — ~16GB at 2B rows, and it
  must be rebuilt on every write. Rejected on size against the scale target.
- **Count row groups on the fly to resolve position** — no extra storage, but
  turns `x[i, ]` into a scan, so it silently makes an O(1)-looking call O(n).
  Rejected as a performance trap.
- **Support positional indexing only on small stores** — a method whose
  availability depends on row count is worse than its plain absence; consumers
  cannot write against it.

## References

- `vignettes/articles/design.Rmd` — *No positional row indexing*, *Scale Targets*
- `AGENTS.md` — *Key Design Decisions*
- `R/methods-accessors.R`, `R/methods-ops.R` (`subset`, `rowSample`)
