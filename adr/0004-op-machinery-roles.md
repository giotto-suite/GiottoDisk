# 0004. Op machinery roles, and the invariants that connect them

- **Status:** Accepted
- **Date:** 2026-08-07
- **Supersedes:** —
- **Superseded by:** —

## Context

The op chain has more moving parts than it has names. Over one working session
on `parquetExprStore`, four separate wrong turns traced to the same cause — a
term meaning two things at different tiers, not a logic error:

- **`@ops` read as "the lowerable ops"** rather than "the prefix that runs
  before materialization". Consequence: the monotonic rule in `.pe_push_op` was
  written up as a defect to be removed. It is not; it is sequencing. An op
  queued after materialization *must* run after it, whatever it is capable of.
- **`norm_libsize` naming an intent** (library normalization) rather than an
  operation (multiply by a per-cell factor). Because the mapping from producer
  to record type was 1:1, `libraryNormParam` could find its own earlier record
  by type and rewrite it — which is only sound if nothing was inserted in
  between. It builds a basis-cut helper (`.pe_chain_before`) purely to serve
  that rewrite; both the rewrite and the helper were removed.
- **"scale" meaning multiply at the operation tier and standardize at the
  workflow tier** — `scaleParam`, the `"scaled"` expression slot and
  `ScaledMatrix` all include centring, while `standardise_flex(scale=)` and
  Halko's `scale=` are the multiplier alone. A `scale` op record would have
  landed in the ambiguous middle.
- **`@post_ops` read as "ops that cannot be lowered"** rather than "the suffix
  that runs after we left Acero", which made a perfectly ordinary placement look
  like a misclassification.

None of these is exotic. They are what happens when a vocabulary is implied by
the code rather than stated.

## Decision

Name the roles, and state the invariants that connect them. The roles:

| role | what it is | instances |
|---|---|---|
| **producer** | a dispatched verb that appends a record | `processData(·, libraryNormParam)`, `processData(·, logNormParam)` |
| **record** | pure-data `list(type, ...params)`, no closures | `multiply`, `log`, `add` |
| **payload** | per-axis state on a record | `factors`, `terms` |
| **carrier** | what values live in while an op runs | lazy Arrow query; collected triplet `data.table` |
| **executor** | applies one record to one carrier | `.op_multiply`, `.op_transform_log` |
| **fold** | composes a whole chain over a carrier | `.pe_apply_ops`, `.pe_apply_post_ops_df` |
| **chain editor** | rewrites the chain as data, touching no values | `.pe_push_op`, `.pe_demote_ops`, `.pe_chain_none` |
| **consumer** | a dispatched verb that reads through the chain to produce a result | `reduceData(·, randomPcaParam)`, `analyzeData(·, featStatsParam)`, `filterData(·, filterParam)`, plus `storeRead()` / `storeWrite()` which are not param-dispatched |

Producers and consumers are named at the **verb** tier — `reduceData(pe,
randomPcaParam(...))`, not `.stream_random_svd()`. That is the unit that
matters when adding an op type: the question is which verbs will encounter it,
and every verb reaches the chain through `storeRead`, so the internal worker it
dispatches to is an implementation detail. The workers are where a chain edit
gets applied when one is wanted — `.pe_pca_demote_chain()` inside the PCA
path — but they are not separate roles.

The invariants, which are the load-bearing part:

1. **Executors are keyed by (record type, carrier), not by phase.** An op that
   is a pure elementwise expression on `value` can have one executor serving
   both carriers — `.op_transform_log` does, because `dplyr::mutate` is generic
   over an Arrow query and a `data.table`. An op carrying per-axis state cannot:
   Arrow has no way to index an R vector from inside a plan, so `multiply` needs
   a joinable table on one side and a positional lookup on the other. That split
   is irreducible and is why those two are separate functions.

2. **Producers append and never revisit.** A record does its work at the
   position it occupies. Reaching back to read or rewrite an earlier record
   assumes nothing was inserted in between, which the chain cannot promise.
   Re-running a producer therefore composes rather than replaces — normalizing
   an already-normalized store yields factors of ~1, which is the right answer
   arrived at by construction rather than by special-casing.

3. **Consumers must not infer capability from slot.** `@post_ops` may hold a
   record that lowers perfectly well; it is there because something ahead of it
   forced materialization. A consumer that wants different execution says so
   with a chain edit, and never by reading the slot and reimplementing its
   contents.

4. **Payloads are keyed by on-disk id, per source.** See 0003. This is what
   makes them invariant under `[`, and it is why the payload for a union is a
   named list of per-substore vectors rather than one stacked vector — `row_id`
   restarts per substore, and `[` can drop a substore entirely.

5. **Chain editors are the only sanctioned way to move where something runs.**
   Bypassing the chain — taking `output = "query"` and re-lowering the records
   by hand — is the anti-pattern this replaces. It duplicates the op's maths at
   the call site with no mechanism keeping the copies in step.

## Consequences

Adding an op type is now a checklist rather than a reading exercise: pick a
primitive name, decide which carriers can express it, write one executor per
(type, carrier) that cannot be shared, and have the producing verb append it.

The instance table above will drift as op types are added. That is expected and
acceptable — the invariants are the durable half, and the registry in
`R/utils-pestore-ops.R` is the authoritative list. This ADR is not the place to
learn what op types exist.

Naming primitives rather than intents costs some readability at the record
level: `multiply` says nothing about *why*. That is deliberate. Intent belongs
to the param class, and a type string that names intent is what let a producer
claim ownership of its own output in the first place.

## Alternatives considered

- **A glossary section in `AGENTS.md`.** Rejected: it would describe current
  behaviour, and these are decisions with rejected alternatives. Also the wrong
  reading distance — you want this when adding an op, not when orienting.
- **Header comment in `R/utils-pestore-ops.R` only.** That header carries the
  instance table and extension protocol, and should. But the invariants are
  choices with alternatives that were tried and discarded, which is what an ADR
  is for; keeping them next to the registry means they get edited as the
  registry changes, which is how the previous framing went stale.
- **A `phase` field on the record.** Considered and rejected in 0002 for
  different reasons; it also breaks invariant 1, since it makes an executor's
  applicability a property of the record rather than of the carrier.
- **Leave the vocabulary implicit.** What was in place. The four confusions in
  *Context* are the argument against.

## References

- `08e6e0a` refactor(pestore): op chain by position, shared stats accumulator,
  derived chunk sizing
- ADR 0002 (phase split), ADR 0003 (payload keying)
- `R/utils-pestore-ops.R` — header: chain model, op registry, extension protocol
- `GiottoClass::standardise_flex()` — the `scale` / `center` naming precedent,
  and the two backends that record an additive vector without densifying
