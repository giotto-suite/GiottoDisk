# 0019. The overlap is a queryableStore carrier, built and aggregated in SQL engines

- **Status:** Accepted
- **Date:** 2026-10-06
- **Supersedes:** —
- **Superseded by:** —

## Context

`calculateOverlap()` on disk stores returned an `overlapPointDisk` wrapping a
`parquetStore`. That class requires an on-disk `row_index`, which no overlap
consumer reads: `overlapToMatrix` aggregates `(feat_ID, poly_ID)`, and joins
back to the points go through `(pt_tile_index, pt_row_index)`. The writers
satisfied the contract in two ways, neither a key. terra wrote `seq_len()` per
tile file into flat files (on Atera, 487,563,530 rows carried 895,750 distinct
`row_index` values). DuckDB wrote `row_number() OVER ()`, which forced the whole
spatial join onto one thread: about 12 minutes on Atera, against about 4 with
the column gone. The op chain, `colnames` and `subset` that came with the class
were never used on an overlap either.

The DuckDB engine also read raw parquet globs, so filters pending on either
store (a cell filter, a feature filter) were ignored, where terra applied them.

Measured on Atera (624M points, 170,057 cells; 487.6M overlap rows; 300.5M
values in the matrix):

| step | engine | elapsed |
|---|---|---|
| overlap | terra, 10 mirai workers | ~12 min |
| overlap | duckdb, `row_number()` | ~12-13 min, one core |
| overlap | duckdb, no `row_number()` | 221 s |
| overlap | sedonadb 0.4.1 | 27.9 s |
| matrix | Arrow (`count`, `arrange`) | 64.4 s, 17.8 GB peak |
| matrix | DuckDB (`GROUP BY`, `ORDER BY`) | 6.6 s, 10.4 GB peak |

Arrow's cost was serial phases: the 300M-group aggregate merged on one thread,
`arrange()` sorted on one thread, and the parquet write ran on the calling
thread. DuckDB ran all three in parallel at about 13 cores. sedona spent 8x
less CPU than DuckDB on the same join. All outputs were checked identical.

## Decision

The overlap is a `queryableStore` over parquet files with a fixed schema:
`poly_ID`, `feat_ID`, any keep columns, `pt_tile_index`, `pt_row_index`. A row
is keyed by `(pt_tile_index, pt_row_index, poly_ID)`; there is no `row_index`.
`calculateOverlap(engine = "duckdb" | "sedona")` joins the two stores'
`storeRead()` scans, so pending filters apply. `overlapToMatrix()` to a
`parquetExprStore` runs one DuckDB query that inner-joins the overlap onto
integer lookup tables built from `@feat_ids` / `@spat_ids`, aggregates, sorts
and writes. Snapshots holding the old `parquetStore` overlap drop it on load.

## Consequences

- Overlaps no longer satisfy the `parquetStore` API: no `[, j]`, `subset()`,
  `rowSample()` or `output = "sedona"` on `@data`. Nothing used them.
- The ID universes are the subset state. Narrowing `@spat_ids` / `@feat_ids`
  narrows the matrix, with no filter on the files. A future `[` on
  `overlapPointDisk` should narrow those slots, not record ops.
- A breaking change: older projects lose their overlaps on load and must rerun
  `calculateOverlap()`.
- DuckDB is the matrix path when installed; without it the Arrow path still
  runs, slower and holding the whole COO.
- Revisit if overlap content ever becomes user-queryable, which would want the
  op chain back.

## Alternatives considered

- **Keep `parquetStore`, write `pt.row_index` as `row_index`.** Fast, but keeps
  a contract nothing reads and the unused API around a scratch result.
- **`fileStore`.** Works for the vault, but drops `storeRead(output =
  "query")`, which the Arrow matrix path and the readers' carriers already use.
- **Shard the overlap by tile and aggregate groups on parallel workers.** Built
  and measured: 36.8 s on 8 forked workers, 143 s serial. Exact, but slower
  than one DuckDB query, and it reorders cells spatially to avoid a global sort.
- **Convert old overlaps on load.** Overlaps are cheap to regenerate and are
  scratch results; conversion code would exist for one release.

## References

- Through the package on Atera with the stored cell and feature filters
  applied (169,420 cells, 487,541,346 overlap rows, 300,480,211 values):
  sedona overlap 38.6 s, duckdb overlap 239 s, `overlapToMatrix()` 17.5 s
  including the store's stats pass. The sedona and duckdb matrices are
  identical, and equal the unfiltered build on the cells both hold.
- `R/methods-aggregate.R`: `.calculate_overlap_sql()`, `.overlap_store()`,
  `.overlap_to_pestore_duckdb()`; `R/methods-snapshotLoad.R`:
  `.drop_legacy_overlaps()`.
- GiottoDisk#76: sedona reads by hive discovery; a per-tile UNION ALL
  overflowed DataFusion's planner stack past a few hundred tiles.
