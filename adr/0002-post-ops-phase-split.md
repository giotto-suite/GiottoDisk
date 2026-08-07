# 0002. Two-phase `@ops` / `@post_ops` split on `parquetExprStore`

- **Status:** Accepted; revision pending (see *Consequences*)
- **Date:** 2026-07-22
- **Supersedes:** —
- **Superseded by:** —

## Context

`parquetExprStore` holds long-format expression as `(row_id, col_id, value,
source_id)` and records workflow steps (library normalization, log, HVF, PCA
prep) as a lazy op chain. Initially every op was Arrow-lazy, applied inside the
composed query.

Measured under giotto-view's streaming consumers, all-Arrow execution cost
**2-3x per chunk** on the normalization arithmetic specifically: Acero
`left_join` against an R-side scale-factor vector plus a `mutate`, versus a
positional R vector index over the collected triplet frame. The structural ops
(filter / distinct / head / tail / sample) showed the opposite profile — there
Arrow's scan optimization is the whole point.

A second problem was schema churn: an Arrow-side normalization could not mutate
`value` in place, so it wrote a `v_norm` sidecar column, and 19+ call sites
branched on `if ("v_norm" %in% names(df))`.

`parquetGeomBase` already had a two-slot shape (`@ops` Arrow-lazy, `@post_ops`
R-side affine transforms on WKB) and was cited as precedent.

## Decision

Split the chain into two slots with a monotonic phase rule:

- `@ops` — Arrow-lazy, composed into the query before collect.
- `@post_ops` — R-side, applied to the materialized `data.table`, for every
  output mode and for streaming consumers alike.

Once `@post_ops` is non-empty, subsequent pushes go to `@post_ops` regardless of
the op's natural phase (`.pe_push_op`). Ops are pure-data records; the phase is
determined by which slot the record lives in, not by a field on it. The escape
valve is materialization: `storeWrite` bakes `@post_ops` into on-disk values and
resets the chain.

The normalization ops (`multiply`, `log`) were moved onto `@post_ops`.

## Consequences

Both original wins hold: the 2-3x on the norm hot path is recovered, and `value`
mutates in place so the schema stays 4-column regardless of op state and the
`v_norm` branches are gone.

The cost is that the slot now carries **two different criteria**. On
`parquetGeomBase`, `@post_ops` means *cannot be lowered* — WKB affine transforms
have no Acero equivalent and performance never enters into it. On
`parquetExprStore` it means *faster in R*: `multiply` and `log` both lower to
Acero fine. There was no way to express "R-side by choice", so the choice got
encoded in the slot — the only place available. Two things follow:

1. **`storeRead` cannot distinguish a lowerable post-op from an R-only one.** A
   consumer wanting pushdown has to route around it via `output = "query"` and
   reimplement the ops, which is why the norm math exists twice with a
   keep-in-sync warning on both copies (`.stream_gene_stats_arrow`).
2. **The monotonic rule is wrong for these ops.** It rejects a lazy op queued
   after anything on `@post_ops`, which is correct for a genuinely R-side op
   (materialization has happened) but not for a norm. This will bite as soon as
   a real lazy pestore op exists — a pushdown `filter` after
   `processData(libraryNormParam)` would error for no reason.

The intended fix restores *expressibility* as the slot's criterion while keeping
the performance win, by making R-side execution a **chain edit** rather than a
permanent placement: demote-with-cascade plus a token op whose executor is the
specialized R implementation. `.pe_demote_ops` is in place; the token swap and
moving the two records back onto `@ops` are not. The pieces do not separate —
moving `multiply`/`log` to `@ops` without the demotion mechanism makes every
materializing read pay the join.

Explicit non-goal for that fix: `storeRead` must not infer a strategy from
output mode. Whether R or Acero wins depends on payload shape and cardinality,
not on the output class.

## Alternatives considered

- **Keep everything on `@ops` (all-Arrow).** Rejected on the measured 2-3x
  per-chunk cost on norm arithmetic, plus the `v_norm` sidecar and its 19+
  branch sites.
- **A `phase` field on the op record instead of two slots.** Would have kept
  both criteria expressible from the start. Not taken: it puts a mutable
  execution hint on what is otherwise a pure-data record, and every consumer
  folding the chain then has to honour it. The chain-edit approach reaches the
  same expressiveness by rewriting the chain, which `show()` can print.
- **Move only `multiply`, leave `log` on `@ops`.** Rejected — the two run
  adjacently on the hot path, so a phase boundary between them forces a collect
  in the middle of the norm.

## References

- `fa3ee64` refactor: phase-split @ops (arrow) + @post_ops (R) execution
- `e120c82` docs: record why the norm ops moved to @post_ops
- `R/utils-pestore-ops.R` header comment (`KNOWN MISCLASSIFICATION`)
- `vignettes/articles/roadmap.Rmd` — *Restore the `@ops` / `@post_ops` contract*
- `vignettes/articles/design.Rmd` — *Lazy Operation System*, *Spatial Transforms
  are Post-Ops*
