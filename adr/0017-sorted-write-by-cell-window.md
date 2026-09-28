# 0017. A parquet expression write sorts per cell window and writes each window itself

- **Status:** Accepted
- **Date:** 2026-09-28
- **Supersedes:** —
- **Superseded by:** —

## Context

`storeWrite(parquetExprStore, <parquetExprStore | union>)` on a lazy chain ended
in `arrange(row_id, col_id)` and handed the query to `arrow::write_dataset()`.
Stores are meant to be cell-major (AGENTS.md), and the `@post_ops` bake path
already wrote that layout by walking cell windows in order.

The lazy path did not. `write_dataset()` writes batches as its threads finish
them and has no order-preserving option in arrow 23, so the sort was paid for and
then discarded. Measured on a 169,420-cell Atera store: 170,944 runs of cells
after a whole-store write, and a different layout on every run. It went
unnoticed because every reader so far either filters by cell range or is additive
over entries. The Gram kernel in GiottoKernels is neither — a cell's entries
must arrive together — so the layout became load-bearing.

The same path was also unbounded: Acero's sort holds the entire output, 17.9 GB
of process peak writing that 300M-value store in one piece.

## Decision

Walk `.pe_windows()` (the same cell windows every bounded pass uses), and per
window run the remap query, sort it, collect it, and write it as its own file,
in cell order. Window size comes from the shared budget model at 96 bytes per
stored value; `giottodisk.chunk_size` pins it. A view that fits the budget is
one window.

## Consequences

- The written layout is cell-major and deterministic, for single stores and
  unions, matching what the bake path already produced.
- Memory is bounded by the window. Whole store (300M values): 17.9 GB peak at
  one window, 10.8 GB at two, 9.3 GB at four, and 22-24 s against 27 s. HVF
  subset (17.6M values): one window, 3.2 s against 2.7 s.
- Each window rescans the source to apply its cell range. On a store of 9,171
  small files that cost ~1.3 s per window (16 windows of the HVF subset: 22.9 s),
  so the budget, not a fixed count, must set the window. Compacting fragmented
  stores lowers it for every reader.
- Sorted output compresses better (720 MB against 1.0 GB for the store above)
  and arrives as one file per window, so a rewrite also compacts.
- `bytes_per_nz = 96` is calibrated, not derived: 60-100 bytes per value of
  process peak, the upper end at more windows because Arrow's allocator keeps
  freed pages. Re-measure if the write path changes shape.
- **Revisit if** arrow's R `write_dataset()` gains an order-preserving write, or
  a reader appears that needs a stronger layout than "each file sorted, files
  covering disjoint cell ranges".

## Alternatives considered

- **Collect the sorted query whole and `write_parquet()` it** — sorted and as
  fast (2.5 s on the HVF subset), but holds the whole output: unbounded for a
  large write, which is the case that matters.
- **Stream the sorted query through a sequential parquet writer** — sorted, but
  still one whole-output sort upstream, and slower (4.5 s).
- **DuckDB `COPY (... ORDER BY row_id, col_id)`** — sorted and correct, and it
  spills: bounded at ~650 MB with a 200 MB limit. But 4.6 s / +3 GB at its default
  limit and 6.3-8.5 s bounded, against 2.5-3.2 s. The same finding as adr/0011:
  a window that never builds the whole sort beats spilling one.
- **Sort once, keep `write_dataset()`, fix the order later** — no such step
  exists short of another full sort.

## References

- `.pestore_write_windowed()`, `.pe_write_chunk_size()` in
  `R/methods-parquetExprStore.R`; the bake path `.pestore_write_baked()`.
- `tests/testthat/test-parquetExprStore-write-layout.R`.
- adr/0011 (cell windowing over spill).
