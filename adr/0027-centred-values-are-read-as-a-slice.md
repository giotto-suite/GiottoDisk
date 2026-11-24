# 0027. Centred values are recorded, and read only as a bounded slice

- **Status:** Accepted
- **Date:** 2026-11-24
- **Supersedes:** —
- **Superseded by:** —

## Context

Giotto's `zscoreScaleParam` and `defaultScaleParam` centre and scale the
expression matrix. On a `parquetExprStore` the scaling half is a sparse
`multiply`, but centring adds a per-axis offset to every entry, implicit zeros
included, so the values are dense. The `add` op had been registered as a
refused stub for this (roadmap, "An `add` op for centred display values").

The consumers of centred values are displays: heatmaps, the `"scaled"` slot.
Every analysis that needs centring already folds it into its own algebra (PCA
with `scale = TRUE`, the gram path's `n·μμᵀ`), and every sufficient-statistic
verb is correct only because an absent entry means 0.

## Decision

`zscoreScaleParam` records its statistics as frozen payloads, computed when the
verb runs (adr/0006): a `multiply` by 1/sd and, with `center = TRUE`, an `add`
of -mean/sd on `@post_ops`. The `add` runs only in a read that hands the fold
its slice, `storeRead(output = "dgcmatrix")`: the frame densifies into triplet
form over the slice and every op after it keeps working on triplets. Every
other consumer — streaming passes, tabular reads — is refused, and a slice
denser than `giottodisk.dense_max` entries is refused. A producer that needs
statistics over a centred chain derives them from the chain before its one
trailing `add`; a chain with more than one `add` must be baked first.

## Consequences

- A centred store cannot feed HVF, markers or other streaming verbs; they error
  rather than compute wrong means. Analyses take the uncentred store.
- `.pe_check_dgc_dims()` bounds one axis only, so `giottodisk.dense_max` is the
  bound on the dense frame.
- Scaling without centring stays sparse and streams. A feature or cell with
  nothing stored reads 0 there, where in memory it is NaN.
- Revisit if a streaming consumer of centred values appears: it would need the
  densification windowed, sized against density 1.0 (`.pe_window_cells()`).

## Alternatives considered

- **Wrap the matrix in `ScaledMatrix` at read time** (the earlier TODO). Keeps
  the slice sparse, but an op after the centring (a `log`, a `compare`) would
  need a matrix-tier executor; densifying into triplets keeps one executor.
- **Densify inside streaming windows.** Exact, but turns every pass over the
  store into a dense pass over all of it, for consumers that do not need it.
- **General algebra for any number of `add` ops.** Each further step adds cross
  terms; two in a row (`defaultScaleParam`) is the case that occurs.

## References

`R/stream-scale.R`, `.op_add_df()` and `.pe_read_slice()` in
`R/utils-pestore-ops.R`, `tests/testthat/test-stream-scale.R`, adr/0006.
