# 0017. GEF ingest writes cell-major: `cellExp` for cellbin, a stripe spill for bin

- **Status:** Accepted
- **Date:** 2026-09-28
- **Supersedes:** 0010 (duplicate deferral; the `coord_env` hand-off stands)
- **Superseded by:** —

## Context

Every expression store is meant to be cell-major (AGENTS.md, "window the CELL
axis"): part files cover ascending cell ranges, so a cell window lowers to a
`row_id` range that prunes row groups. The mtx, 10x h5, zarr and CSV inputs all
write that way. The GEF readers did not: they streamed the gene-major
`geneExp` in batches of about 500 genes and wrote one part per batch, each
spanning every cell. Cell windows then pruned almost nothing, so every windowed
pass rescanned the store once per window, and marker detection on a small
memory budget got slower and heavier instead of lighter.

A cellbin GEF also carries a cell-major copy of the matrix, `cellBin/cellExp`,
indexed by `cellBin/cell$offset` / `geneCount`. A bin GEF has none.

## Decision

- **Cellbin** reads `cellExp` in contiguous cell ranges sized by the memory
  budget. Files without it take the bin path's spill.
- **Bin** streams `geneExp/<bin>/expression` in gene chunks as before, names
  bins in first-appearance order (the in-memory reader's numbering), and spills
  each record to a horizontal stripe of the chip. The stripes are emitted in
  ascending y; within one, bins take store positions in (y, x) order. The
  published coordinates gain a `pos` column so spatial locations follow store
  order.
- Duplicate gene names are summed by each batch's `(row_id, col_id)`
  aggregate. A batch owns whole cells, so there is nothing left to defer.

## Consequences

- Stores are cell-major from every importer; windowed passes may rely on it.
  A new input that cannot read in cell order must reorder before writing.
- Bin ingest writes and reads the matrix once more, through a temporary spill
  under `tempdir()`, sized like any other batch. Time still fell, because the
  reads were batched per gene either way and the store now has a few parts
  instead of one per gene batch.
- A bin's name (`bin_<first appearance>`) and its store position differ.
  Anything that assumes `bin_k` sits at position k is wrong; use names.
- Stores written before this keep their gene-batched layout until re-imported.
- This spill does not reopen adr/0011, which chose windowing over spilling for
  a grouped join on a store that is already cell-major. Here the input is not,
  and reordering it once at ingest is what lets those windows prune.
- The `coord_env` hand-off from adr/0010 carries over unchanged.

## Alternatives considered

- **Detect the layout at read time and stop windowing** — tried as a stopgap
  (footer statistics via nanoparquet). It fixed the symptom in one pass and left
  every other windowed pass and the stores themselves as they were.
- **Number bins in grid order** — simpler, but bin names would then differ from
  the in-memory reader's for the same file.
- **Densify `wholeExp/<bin>` to assign positions up front** — exact, but bin1
  is a 23,520² grid; the stripe spill needs only the coordinate range.

## References

`R/methods-fileInputs.R` (`storeRead` for `cellbinGefInput` / `binGefInput`),
`R/utils-parquetExprStore.R` (`.gef_spill_iterator`, `.gef_cell_ranges`,
`.gef_batch_rows`), `tests/testthat/test-stereoseq-gef.R`.
