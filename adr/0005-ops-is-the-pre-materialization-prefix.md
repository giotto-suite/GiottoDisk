# 0005. `@ops` is the pre-materialization prefix, not the set of lowerable ops

- **Status:** Accepted
- **Date:** 2026-08-07
- **Supersedes:** 0002
- **Superseded by:** —

## Context

0002 split the chain into `@ops` (Arrow-lazy) and `@post_ops` (R-side) and moved
the normalization records onto `@post_ops`, on a measured 2-3x per-chunk win for
a positional R vector index over an Acero `left_join`. It recorded its own cost
honestly: the slot ended up carrying two criteria — *cannot be lowered* on
`parquetGeomBase`, *faster in R* here — and named the intended fix as
demote-with-cascade plus a token op, with `.pe_demote_ops` as the first half.

Two things changed since.

**The payload changed shape.** Under 0003 a `multiply`'s factors became a
per-substore vector keyed by on-disk id. The R executor no longer rebuilds a
lookup from a `(source_id, orig_row_id, scalef)` table on every call — it
indexes the vector directly — and the Arrow executor builds its joinable table
from the same vectors. The 2-3x that justified the placement was measured
against the table payload.

**The consequences 0002 predicted arrived.** Consumers wanting pushdown routed
around the slot: `.stream_gene_stats_arrow` took `output = "query"` and
re-lowered the records itself, so the norm maths existed twice with a
keep-in-sync warning on both copies. When that workaround was removed in favour
of honouring the chain, the stats verbs slowed **10x** — measured 0.154s against
1.66s on 19.2M nonzeros — for no reason other than records being parked in the
slot that means "we have left Acero".

The token op in 0002's intended fix turned out to be unnecessary. Its purpose
was to name a specialized R implementation, but the strategy choice (positional
index versus keyed join) is derivable from the payload, so the record already
carries everything the executor needs.

## Decision

`@ops` and `@post_ops` are **one ordered sequence split at the point where
execution leaves Acero**. `@ops` is the prefix; `@post_ops` is everything from
the first step that cannot run in Acero onward. Lowerability is a *consequence*
of a record's position, not the definition of the slot — a record that lowers
perfectly well sits on `@post_ops` legitimately when something ahead of it
forced materialization.

Producers push `phase = "lazy"`. A consumer that wants a record run R-side
expresses it as a chain edit — `.pe_demote_ops` moves that record and, of
necessity, everything after it — rather than by permanent placement or by
bypassing the chain. No token op.

The monotonic rule in `.pe_push_op` stands unchanged and is **not** a defect. It
was described in 0002 as "wrong for these ops"; it is not. `@ops` runs entirely
before the collect, so an op appended after something on `@post_ops` would
execute *before* it. What was wrong was the placement that made the rule fire.

## Consequences

`storeRead` pushes the whole chain down by default, so every consumer gets
Acero without asking. The duplicated norm maths and its keep-in-sync warning are
gone; `.stream_gene_stats_arrow` no longer exists.

The PCA band loops keep their hot path by demoting at the first `multiply`
(`.pe_pca_demote_chain`) — measured at parity with the pre-refactor path, so the
2-3x from 0002 is preserved where it was actually earned. This is
`.pe_demote_ops`' first real caller; before this it was exercised only against
synthetic records.

`storeWrite` now takes its fast branch for a normalized store, folding the chain
into the same Arrow plan as the row/col remap instead of running the chunk-stream
bake.

Both of 0002's original wins survive: `value` still mutates in place (no `v_norm`
sidecar), and the R-side positional index is still what the band loops use.

The constraint this places on future code: a consumer must not decide where a
step runs by inspecting which slot it is in. If it wants different execution it
edits the chain, and the edit is visible in `show()`.

Revisit if a consumer appears whose preferred execution differs *per record*
rather than *per suffix* — demotion cascades by construction, so it cannot
express "run this one in R but the next one in Acero". Nothing needs that today,
and materialization being one-way suggests nothing can.

## Alternatives considered

- **Keep the placement, restore the workaround.** Rejected: it is the duplicated
  maths that 0002 already flagged as the cost, and the copies have no mechanism
  keeping them in step.
- **Token op, as 0002 intended.** Rejected as unnecessary rather than wrong. The
  strategy choice is derivable from the payload, so a distinct type would name
  something the executor can already see. It becomes load-bearing only if one
  record type ever gets two R-side implementations.
- **Infer the phase from output mode in `storeRead`** — demote for materializing
  reads, lower for `query`. Rejected in 0002 and still rejected: whether R or
  Acero wins depends on payload shape and cardinality, not output class.
- **Revert 0002 entirely and put everything back on `@ops`.** Rejected: the band
  loops genuinely are faster with the positional index, and 0002's measurement
  stands. Demotion keeps that win without encoding it in the slot.

## References

- `08e6e0a` refactor(pestore): op chain by position, shared stats accumulator,
  derived chunk sizing
- ADR 0002 (superseded), ADR 0003 (payload keying), ADR 0004 (roles and
  invariants)
- Benchmark: `featStats` 1.66s (`@post_ops`) vs 0.154s (`@ops`), 4000 x 40000
  at density 0.12; PCA at parity with `e120c82`
- `R/stream-pca.R` — `.pe_pca_demote_chain()`
