# 0023. Sedona reads CRS-less stores in a planar engineering CRS; files stay `crs: null`

- **Status:** Accepted
- **Date:** 2026-10-28
- **Supersedes:** —
- **Superseded by:** —

## Context

Spatial omics coordinates are planar and have no CRS; `@params$crs` is `""`.
Geometry stores write GeoParquet with `crs: null` (unknown) and
`edges: planar`. An omitted `crs` would mean OGC:CRS84 under the spec.

SedonaDB requires both sides of a predicate to share a CRS, and decides itself
what `crs: null` means: OGC:CRS84 in 0.4, no CRS in 0.5. The suite followed it,
copying the column's CRS onto query geometries (`ST_Crs(geom)`). Two things
broke that:

- Points without WKB (adr/0022) are built in SQL, so their CRS has to be chosen.
  Copying the polygons' per row (`ST_Crs(poly.geom)`) turns the overlap join
  into a nested loop; it has to be a constant. A constant that follows sedona
  needs a per-session probe of what sedona does with `crs: null`.
- A lng/lat label claims geographic coordinates. Predicates stay planar today
  because of `edges: planar`, but an operator that grows CRS awareness
  (area, length, distance in metres on an ellipsoid) would change results
  without an error, and Giotto coordinates run far outside lng/lat ranges.

Tested on sedonadb 0.4.1.9000: an engineering CRS given as PROJJSON is accepted
by `ST_SetCrs` on file geometry over multi-file scans and by the 3-argument
`ST_Point(x, y, crs)`; the overlap of 60 polygons and 4,000 points matches
terra (690) and plans as a spatial join.

## Decision

The sedona base view puts every store without a CRS in one fixed PROJJSON
engineering CRS, `.SEDONA_PLANAR_CRS`: `ST_SetCrs(geom, …)` for stored geometry,
`ST_Point(x_index, y_index, …)` for points. Query geometries copy it with
`ST_Crs(geom)`. Files keep `crs: null`; the CRS exists only in the engine
adapter.

## Consequences

- No dependence on how a sedona version reads `crs: null`, and no probe.
- Files commit to nothing: when the suite gets a real CRS concept, this is one
  function to change, with no artifacts to migrate.
- A store with a real `@params$crs` is not overridden; points in such a store
  get no CRS in SQL. No such stores exist yet.
- The sedona dataframe from `storeRead(output = "sedona")` carries this CRS, so
  user SQL that joins it with its own geometry must set the same CRS. It did
  already: 0.4 gave the column OGC:CRS84.
- sedona reads a `geom` column whose file statistics have min == max (one row,
  or identical shapes) as plain binary, and every geometry function fails on
  it. That predates this decision; the override adds a plain collect of such a
  store to what fails.
- DuckDB and terra enforce no CRS and are unaffected.

## Alternatives considered

- **Probe sedona's reading of `crs: null` once per session** and use that
  constant. Works, but keeps the lng/lat label on 0.4 and stays version
  dependent.
- **Store the sedona CRS in `@params`.** It describes the sedona version doing
  the reading, not the data, so a stored value goes stale on upgrade.
- **Write the engineering CRS into the files.** Spec-valid, and outside tools
  would see a local frame, but every artifact would commit to a choice made
  before the suite has CRS support.
- **`EPSG:3857`, or a made-up `GIOTTO:PLANAR`.** Both work in 0.4. Web Mercator
  is a real projection that a CRS-aware operator could act on; a made-up code
  works only while sedona compares CRSs as strings.

## References

- `.pgeom_view_sql()`, `.SEDONA_PLANAR_CRS` (`R/methods-storeRead.R`);
  `.sedona_query_geom_sql()` (`R/methods-spatRelate.R`);
  `.geoparquet_metadata()` (`R/utils-spatial.R`).
- adr/0022.
