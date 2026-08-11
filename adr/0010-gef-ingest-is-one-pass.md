# 0010. GEF ingest is a single pass: duplicates defer, coordinates ride out

- **Status:** Accepted
- **Date:** 2026-08-11
- **Supersedes:** —
- **Superseded by:** —

## Context

A Stereo-seq `.gef` stores expression as one gene-major HDF5 compound dataset
(`cellBin/geneExp`, or `geneExp/<bin>/expression`). `cellbinGefInput` and
`binGefInput` stream it in gene chunks via rhdf5 hyperslabs, which is the whole
point — a `tissue.gef` bin1 table is ~14M records and a full `.gef` ~25M, and
the backend exists so none of that has to be resident.

Streaming imposes two constraints that in-memory readers never meet, because
they hold the entire table and can do whatever they like to it.

**Duplicate gene names.** Two gene rows with distinct `geneID`s can carry the
same `geneName`, and their records must sum into one matrix entry. The gene
table is ordered by `geneID`, so the duplicates of a name scatter arbitrarily:
on `C04687E314.tissue.gef`, 16 duplicated names with gaps up to 25,400 rows.
`.gef_safe_chunks()` only keeps *consecutive* runs of the same name together —
it cannot reach a partner 25,000 rows away without collapsing the chunk plan
into one enormous chunk. Left alone, each chunk aggregates independently and
the store receives two rows for one `(cell, gene)` pair: 753 such rows on that
file, inflating `nnz` and every marginal derived from it. The in-memory reader
never sees this, because `Matrix::sparseMatrix()` sums duplicate `(i, j)` pairs
for free.

**Bin coordinates.** For `type = "bin"` there is no cell table. A bin's `(x, y)`
exists only inside the expression records, and `bin_ID` is assigned by first
appearance while streaming. The iterator therefore builds the `(x, y) -> bin_ID`
map as a side effect of ingest. Giotto's in-memory `.stereoseq_spatlocs()`
recovers the same map by re-reading the whole `geneExp/<bin>/expression`
dataset — the exact read the backend exists to avoid. But `storeWrite()` returns
an S4 store, and the `exprInput` is copied on the way in, so there is no obvious
channel to carry a side value back out to the caller.

## Decision

Both are solved inside the GEF readers, without touching the shared
`storeWrite(parquetExprStore, exprInput)` consumer.

**Duplicates defer.** `.gef_dup_cols()` identifies every `col_id` backed by more
than one raw gene row. Records for those columns are held back as each chunk is
emitted and flushed once, aggregated, as a final batch at end of stream
(`.gef_flush_deferred()`). Memory is bounded by the duplicated genes alone, not
the matrix.

**Coordinates ride a reference cell.** `binGefInput()` carries
`params$coord_env`, an environment. The iterator writes the finished map into it
on a complete pass, and `.stereoseq_expression_disk()` reads it back after
`sourceWrite()` returns. An environment because copy semantics defeat a plain
slot; publication is gated on exhaustion because a partial map would silently
yield spatial locations for a subset of the store's columns.

## Consequences

- The GEF iterators may emit one batch after the last chunk. Anything driving
  `next_batch()` must loop to `NULL` rather than counting chunks.
- `.gef_safe_chunks()`'s adjacency handling is now redundant but harmless; it
  survives as the chunk-size planner. Do not "fix" it to chase non-adjacent
  duplicates — that is what the deferral is for, and widening chunks to span
  them destroys the streaming property.
- `coord_env` looks like dead state on the input object. It is the only path by
  which bin spatial locations avoid a second full read of the gef; deleting it
  silently reintroduces that read via the inherited in-memory closure.
- Neither mechanism is reachable from `mtxInput` / `tenxH5Input` /
  `csvWideInput`. `storeWrite` is unchanged, so no other input type pays for
  either.
- Revisit if the GEF format ever guarantees unique `geneName`s, or if
  `storeWrite` grows a sanctioned way to return per-input metadata — at which
  point `coord_env` should move onto it.

## Alternatives considered

- **Capture the coordinates in `storeWrite`** — the original plan. It returns
  the S4 store, so the value has nowhere to go without either an attribute on
  an S4 object (fragile) or a second return path in a consumer shared by every
  input type. Rejected: it puts GEF-specific state on the universal seam.
- **Widen chunks to contain every duplicate group** — correct, but a 25,400-row
  gap forces a chunk spanning most of the gene axis, which is the streaming
  property gone.
- **Dedup the written Parquet afterwards** — a full extra read/write of the
  store to fix ~753 rows out of 5.2M.
- **Re-read the gef for bin spatial locations**, as the in-memory reader does —
  correct and simple, and what the code did before. It makes `backend =` a
  memory no-op on the bin path, which is the path that needs it most.

## References

- `R/utils-parquetExprStore.R` — `.gef_dup_cols()`, `.gef_flush_deferred()`,
  `.gef_safe_chunks()`
- `R/methods-fileInputs.R` — `storeRead()` for `cellbinGefInput` / `binGefInput`
- `R/class-fileInputs.R` — `binGefInput()`, `params$coord_env`
- `R/convenience-stereoseq.R` — `.stereoseq_expression_disk()`,
  `.stereoseq_spatlocs_from_coords()`
- `tests/testthat/test-stereoseq-gef.R` — cross-chunk duplicate summing,
  first-appearance bin-ID parity, no publication from a partial pass
- Validation on `C04687E314`: both GEF flavours match the in-memory reader
  exactly on nnz, IDs, per-gene and per-cell sums, and spatial locations.
