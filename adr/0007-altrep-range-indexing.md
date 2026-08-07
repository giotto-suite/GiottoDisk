# 0007. HDF5 slice reads index with `lo:hi`, not `seq.int()`

- **Status:** Accepted
- **Date:** 2026-05-30
- **Supersedes:** —
- **Superseded by:** —

## Context

The 10x / Xenium HDF5 ingest reads CSC slices out of `matrix/data` and
`matrix/indices` one cell batch at a time, through `hdf5r`'s `[` method.

`hdf5r` inspects the index it is handed. Given a contiguous integer sequence it
issues a single **hyperslab** read; given an arbitrary integer vector it falls
back to **point selection**, one element at a time. It detects the former by
recognising R's ALTREP compact sequence — the object `:` produces, which stores
only its bounds.

`seq.int(lo, hi)` does not produce one. It materialises a real integer vector.
The two expressions are indistinguishable at the call site, return identical
values, and `seq.int()` is arguably the more explicit form — but only one of
them is served as a hyperslab.

Measured on Atera, a 25k-cell batch, byte-identical output either way:

| index form | time | RSS |
|---|---|---|
| `lo:hi` | 0.27 s | +0.55 GB |
| `seq.int(lo, hi)` | 1.20 s | +2.16 GB |

There are two such reads per batch. At 7 concurrent ingest workers this single
distinction accounted for roughly **15 GB of peak RSS and about half the ingest
wall-clock**.

## Decision

Index ALTREP-aware readers with `lo:hi`. Never `seq.int(lo, hi)`, and never a
materialised vector where a range would do.

This applies wherever the reader inspects the index to choose an access
strategy — `hdf5r` today, and the same trap exists for any backend that
special-cases contiguous selections.

The distinction is invisible in the code, so the call sites carry a comment
saying which form is required and why. Those comments are load-bearing, not
commentary.

## Consequences

Two call sites in `R/methods-fileInputs.R` are constrained: the `.tenxh5_*`
readers must keep the `:` form. A tidy-up pass that "clarifies" them to
`seq.int()` gives back half the ingest throughput with no test failing and no
output changing — which is precisely why this is written down rather than left
to the comments alone.

The constraint does not extend to ordinary R vector construction. `seq.int()` is
fine, and preferable, anywhere the result is used as a value rather than handed
to a reader as an index. `rep.int(seq.int(c_lo, c_hi), nnz_per_cell)` in the
same functions is correct as written.

Nothing enforces this mechanically. A linter rule matching `\[seq\.int\(` on
`H5D` objects would be possible but is not worth the machinery for two sites;
the risk is a well-meant edit, and the mitigation is that the reason is now
findable.

Revisit if `hdf5r` gains contiguity detection for materialised vectors, or if
the ingest path moves off `hdf5r` entirely.

## Alternatives considered

- **`seq.int(lo, hi)` for explicitness.** This is what the trap looks like. It
  reads better and costs 4x in both time and memory.
- **Pass an explicit hyperslab spec to `hdf5r`** rather than relying on ALTREP
  detection. Rejected: more API surface at every call site to express what `:`
  already expresses, and it would still be silently correct-but-slow if someone
  bypassed it.
- **Wrap the reads in a helper that normalises the index.** Reasonable, and
  worth doing if a third call site appears. For two sites it adds a layer
  between the reader and the code without removing the underlying trap.

## References

- `R/methods-fileInputs.R` — `.tenxh5_read_range()` (primary note),
  and the `.tenxh5_*` batch reader that back-references it
- Measurement: Atera FFPE, 25k-cell batch, 7 concurrent workers
