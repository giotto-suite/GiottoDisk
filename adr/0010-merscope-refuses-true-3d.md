# 0010. The MERSCOPE reader refuses true 3D rather than reducing it

- **Status:** Proposed
- **Date:** 2026-08-10
- **Supersedes:** —
- **Superseded by:** —

## Context

MERSCOPE writes `cell_boundaries.parquet` with one row per cell *per z-plane*.
Two very different exports share that shape. In the common case the same 2D
segmentation is replicated across every plane — identical `EntityID` sets and
identical geometry, differing only in `ZIndex`. In the other, the planes carry a
genuine 3D segmentation and differ by design.

They cannot be told apart from the schema, and reading the replicated case
naively multiplies every cell by the plane count: a 7-plane export of 11,955
cells presents as 83,685 rows, and every per-cell aggregate is wrong by 7x
without anything looking obviously broken.

Giotto's in-memory `createGiottoMerscopeObject()` already solves the detection
half, comparing the first 100 vertex coordinates between two planes. It also has
somewhere to *put* true 3D: `aggregateStacks()` combines per-plane spat_units
into one aggregate unit. GiottoDisk has no equivalent — aggregating expression
across per-plane spat_units is a store-level operation, and none is implemented.

So the detection could be ported. The disposal of the 3D case could not.

## Decision

Detect the architecture, then diverge from Giotto on what happens next.

Detection ports Giotto's v2 logic with a stricter test: compare `EntityID` sets
across planes, then byte-equality of the serialized WKB geometry for a sample of
cells. Identical bytes imply identical vertices and nothing has to be decoded.

Disposal splits three ways:

```
2d              single plane, or no ZIndex column   -> use it
replicated_2d   IDs and geometry identical          -> keep lowest plane, report the drop
3d              either test fails                   -> stop()
```

The third case is the decision. The reader **refuses** genuine 3D rather than
ingesting one arbitrary plane, and the error names the available z-indices and
the explicit `poly_z_indices` argument that overrides it.

## Consequences

- Users with true 3D data are blocked rather than silently served a partial
  object. `poly_z_indices = <one index>` is the documented escape hatch, and the
  error says so.
- Detection is not free: distinct `ZIndex` values, then two full `EntityID`
  columns, then geometry for up to `n_sample` cells, before any real work.
- **The 2D verdict rests on a sample.** `n_sample = 200L` under `set.seed(1L)`.
  An export replicated across the 200 sampled cells but differing elsewhere is
  misclassified as replicated-2D, and the extra planes are dropped without
  warning. Raising `n_sample` narrows that window at linear cost; it is a
  confidence/latency trade, not a proof.
- The `.merscope_stop_3d()` error message is load-bearing rather than
  decorative, since refusing is only defensible if the way forward is legible
  from the message.
- Revisit when a disk-backed `aggregateStacks()` exists: the refusal should then
  become a delegation, and this ADR is superseded rather than amended.

## Alternatives considered

- **Always use z0, no detection** — cheapest and correct for the common export,
  silently discards most of the segmentation on a real 3D one. The failure is
  invisible: a plausible object with the wrong cells in it.
- **Load every plane as-is** — loses no data, but multiplies cell counts by the
  plane count on the *common* case, pushing a correctness problem downstream to
  every aggregate.
- **A user argument with no detection** — makes the caller answer a question they
  would have to inspect the parquet to answer, and gives a wrong default to
  whoever does not.
- **Reduce 3D by picking one plane with a warning** — rejected as the same class
  of failure as a silent wrong answer. Warnings scroll past; the object survives
  and gets analysed.
- **Compare all cells rather than a sample** — exact, but reads the full geometry
  column purely to make a routing decision, on the largest column in the file.

## References

- `R/convenience-merscope.R`: `.merscope_zplane_architecture()` (detection and
  the `n_sample` constant), `.merscope_stop_3d()` (the refusal and its message).
- `tests/testthat/test-merscope.R` — fixtures cover `replicated`, `3d` and
  single-plane. The 3D branch is only reachable with synthetic fixtures: every
  real export seen so far is replicated 2D.
- Giotto's `createGiottoMerscopeObject()` for the in-memory detection this ports,
  and `aggregateStacks()` for the 3D path GiottoDisk lacks.
- Verified on `202601301420_FFPE-MsSkin-...-Cdkn2a-vs311-LM`, region R1: 7 planes,
  detected replicated-2D, 6 dropped, 11,955 cells ingested. Transcript totals
  after overlap matched the vendor `cell_by_gene.csv` exactly (1,646,708).
