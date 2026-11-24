# 0026. Sedona reads CRS-less stores in a planar engineering CRS; files stay `crs: null`

- **Status:** Proposed
- **Date:** 2026-11-23
- **Supersedes:** —
- **Superseded by:** —

## Context

Spatial omics coordinates are planar and have no CRS; `@params$crs` is `""`.
Geometry stores write GeoParquet with `crs: null` (unknown) and
`edges: planar`. An omitted `crs` would mean OGC:CRS84 under the spec.

SedonaDB requires both sides of a predicate to share a CRS, and decides itself
what `crs: null` means: OGC:CRS84 in 0.4, no CRS in 0.5. The suite followed it,
copying the column's CRS onto query geometries (`ST_Crs(geom)`). That breaks in
three ways:

- Geometry built in SQL has to be given a CRS. The area rule for bin points
  builds each bin's square and centre (`ST_MakeEnvelope`, `ST_Point`) and
  compares them with the polygons. Built with no CRS, they mismatch the
  polygons' OGC:CRS84 on 0.4 (`Mismatched CRS arguments: ogc:crs84 vs None`).
  Copying the polygons' per row (`ST_Crs(poly.geom)`) turns the join into a
  nested loop; it has to be a constant. A constant that follows sedona needs a
  per-session probe of what sedona does with `crs: null`.
- Following sedona means the answer changes when sedona does. 0.5 reports no
  CRS for these stores, and query geometries tagged with `ST_Crs(geom)` stopped
  matching anything.
- A lng/lat label claims geographic coordinates. Predicates stay planar today
  because of `edges: planar`, but an operator that grows CRS awareness
  (area, length, distance in metres on an ellipsoid) would change results
  without an error, and Giotto coordinates run far outside lng/lat ranges.

Tested on sedonadb 0.4.1.9000: an engineering CRS given as PROJJSON is accepted
by `ST_SetCrs` on file geometry over multi-file scans and by the 3-argument
`ST_Point(x, y, crs)`; `ST_SetCrs(ST_Point())` fails to plan over multi-file
scans.

## Decision

The sedona base view puts every geom store without a CRS in one fixed PROJJSON
engineering CRS, `.SEDONA_PLANAR_CRS`, and a store with a CRS in its own:
`ST_SetCrs(geom, …)` for stored geometry, and `ST_Point(x_index, y_index, …)`
for points. Query geometries copy it with `ST_Crs(geom)`; geometry built in a
join takes it as a constant (`.sedona_crs_sql()`). Files keep `crs: null`; the
planar CRS exists only in the engine adapter.

Whether a store has a CRS is decided from our side, never by asking sedona:
`@params$crs` when set, else the `crs` its files record (a store opened from a
path has no params). Asking sedona (`ST_Crs(geom)`) cannot tell an intended
CRS from its own reading of `crs: null`, which is OGC:CRS84 on 0.4. The planar
constant claims no unit (`unity`) and no axis orientation.

Points are rebuilt from their coordinates because sedona reads a WKB column
whose file statistics have min == max (a one-point tile) as plain binary, and
`ST_SetCrs` has no kernel for binary.

## Consequences

- No dependence on how a sedona version reads `crs: null`, and no probe.
- Files commit to nothing: the planar CRS exists only in the adapter, so
  changing it migrates no artifacts.
- A store with a real CRS is read in it, and geometry built against it takes
  it too.
- A store with a CRS compared with one without (or with a different CRS) is a
  sedona error, not a silent overlap. No such pairs exist yet.
- The sedona dataframe from `storeRead(output = "sedona")` carries this CRS, so
  user SQL that joins it with its own geometry must set the same CRS. It did
  already: 0.4 gave the column OGC:CRS84. Its `geom` is now the last column.
- A polygon store whose `geom` statistics have min == max (one row, or
  identical shapes) is read as plain binary. Geometry functions already failed
  on it; a plain sedona collect of one now fails too.
- DuckDB and terra enforce no CRS and are unaffected.

## Alternatives considered

- **Probe sedona's reading of `crs: null` once per session** and use that
  constant. Works, but keeps the lng/lat label on 0.4 and stays version
  dependent.
- **Probe per view** (`ST_Crs(geom)` on one row) and tag query geometries only
  when the store has a CRS. Fixes the 0.5 case, but not geometry built in a
  join on 0.4, and adds a query per view.
- **Store the sedona CRS in `@params`.** It describes the sedona version doing
  the reading, not the data, so a stored value goes stale on upgrade.
- **Write the engineering CRS into the files.** Spec-valid, and outside tools
  would see a local frame, but every artifact would commit to a choice made
  before the suite has CRS support.
- **`EPSG:3857`, or a made-up `GIOTTO:PLANAR`.** Both work in 0.4. Web Mercator
  is a real projection that a CRS-aware operator could act on; a made-up code
  works only while sedona compares CRSs as strings.

## References

- `.sedona_view_sql()`, `.sedona_crs_sql()`, `.pgeom_store_crs()`,
  `.SEDONA_PLANAR_CRS`
  (`R/methods-storeRead.R`); `.sedona_query_geom_sql()`
  (`R/methods-spatRelate.R`); `.bin_area_sql()` (`R/methods-aggregate.R`);
  `.geoparquet_metadata()` (`R/utils-spatial.R`).
- adr/0019 (the overlap join), adr/0022 (the bin area rule).
