# 0022. A bin goes to one polygon: the one holding its centre, else the one covering most of it

- **Status:** Proposed
- **Date:** 2026-11-04
- **Supersedes:** —
- **Superseded by:** —

## Context

Overlaps of bin points with cell polygons tested only the bin's point
("inside"). Stereo-seq stores a bin1 point at the corner of its DNB, and the
cellbin GEF gives a cell every DNB its mask covers, so testing the point loses
DNBs at every cell edge. On the Stereo-seq mouse-eyeball mini (C04687E314,
1,368 cells, 3.44M bin1 records), "inside" recovered 0.950 of the GEF's cell
counts (median cell 0.959). The assignment has to come out the same in memory
(terra) and in the three disk engines (terra tiles, DuckDB, SedonaDB), and each
bin has to be counted once.

## Decision

A bin is the square of side `bin_width` around its point plus `bin_offset`.
It goes to the one polygon sharing interior with the square that holds the
square's centre or, when none or several do, to the one covering most of the
square, ties to the lowest `poly_ID` (byte order). Edge contact alone is no
hit. GiottoClass's `extract_bin_overlap()` is the in-memory rule and the terra
engine's; the SQL engines run the same rule in one query. The bin points carry
their geometry and default rule (`@bin_width`, `@bin_offset`,
`@bin_overlap`); the Stereo-seq readers set `"area"`, everything else keeps
`"inside"`.

## Consequences

- On the mini, "area" recovers 0.983 of the GEF counts (median cell 0.989);
  memory, terra, DuckDB and SedonaDB give identical assignments (506,816 bins).
- The engines agree only on valid polygons: area and touch predicates on a
  self-intersecting outline are undefined, and each engine answers
  differently (9 of the mini's 1,368 cellbin outlines). Readers repair
  invalid outlines by default (GiottoClass option `giotto.make_valid`); with
  the repair turned off, "area" can differ between engines on those cells.
- The terra engine tiles over the bins for "area" (each bin in one tile,
  polygons fetched padded by their radius plus a bin), the reverse of the
  "inside" path, so no bin is decided twice.
- The SQL query adds window functions partitioned by bin, which ADR 0019 warns
  about for global windows; not yet measured at Atera scale.
- `binOverlapsTable()` reads bin -> polygon pairs off the overlap. That is a
  read of finished content, not the op chain ADR 0019's revisit clause means.
- Revisit if a platform's bins are not axis-aligned squares in object
  coordinates (Visium HD in pixel space can be rotated), or if exact GEF
  agreement is required: GDAL's "all touched" rule (raster) reached 0.990
  (median 0.996) by counting edge contact on one side, which the SQL engines
  cannot express.

## Alternatives considered

- **Count edge contact too (closed squares).** 1.032 of the GEF counts
  (median 1.046): every DNB next to a cell joins it.
- **Rasterise the polygons with `touches = TRUE`.** Closest to the GEF, but
  the raster's draw order decides shared DNBs and the SQL engines cannot
  reproduce it.
- **Move Stereo-seq points to the DNB centre and keep "inside".** Changes the
  coordinates every other Stereo-seq output uses, and still drops DNBs a cell
  outline crosses.
- **Settle bins two tiles both took after the fact.** Read every overlap row
  and all polygons into R, then rewrote the files; unbounded at bin1 scale.

## References

- GiottoClass `extract_bin_overlap()`, `calculateOverlap(giottoPolygon,
  giottoBinPoints)`; `R/methods-aggregate.R`: `.bin_area_sql()`,
  `.calculate_overlap_terra_bins()`.
- `tests/testthat/test-binpoints-disk.R`: memory and every engine agree,
  including a tie where `poly_ID` order differs from row order.
