# 0022. A points store is its coordinates; WKB is an opt-in export column

- **Status:** Accepted
- **Date:** 2026-10-28
- **Supersedes:** —
- **Superseded by:** —

## Context

Every geometry store wrote `x_index` / `y_index` (the centroid) and a WKB `geom`
column. For points the two carry the same numbers, and the WKB is the expensive
copy:

- **Size.** 5M points, `id + x_index + y_index` with and without `geom`:
  187 → 107 MB (snappy), 145 → 83 MB (zstd), about 42% smaller.
- **terra.** `terra::vect()` on 1M points takes 0.90 s from WKB, 0.14 s from
  x/y.
- **Engines.** A predicate on a point has to build the point either way, by
  decoding WKB or by `ST_Point(x, y)`, and parquet min/max statistics on
  `x_index` / `y_index` prune row groups as well as a GeoParquet bbox does,
  since writes are Morton-sorted. Earlier benches found WKB bought spatial SQL
  engines little over raw numerics.
- `terra::centroids()` already read a polygon store as points from
  `x_index` / `y_index` (`use_xy_as_geom`), so most of the reading path existed.
- CRS is not yet a real concept for spatial omics data (`@params$crs` is `""`).

The cost is interoperability: a parquet file with no geometry column is not
GeoParquet, so GDAL, QGIS and sedona's GeoParquet reader do not see points in it.

## Decision

Points stores write no `geom` column and no GeoParquet metadata. Every reader
builds points from `x_index` / `y_index` and ignores a `geom` column if one is on
disk. `options(giottodisk.points_wkb = TRUE)` writes `geom` again, for tools
outside the suite; the suite still never reads it. Points stay a
`parquetGeomStore` with `@geomtype = "points"`, not a class of their own.

## Consequences

- Code that names `geom` must not assume it exists on a points store. The SQL
  engines get it back from the base view (`.pgeom_view_sql()`); R code that
  selects special columns intersects them with the schema first.
- Stores written before this still read, and so do stores written with the
  option: `geom` is ignored on both.
- On sedona the built point needs the same CRS as the polygons it meets;
  adr/0023 sets one for both.
- No Z yet. A `z_index` column would follow the same rule, but the Morton sort
  is 2D, so it needs its own decision.

## Alternatives considered

- **A separate `parquetCoordsStore` class.** The hierarchy already varies two
  ways (single file vs union, tiled vs not), so a new class needs tiled and
  union versions, or a new virtual base. A centroid view is a points store over
  polygon files, and with one class that is just `@geomtype`; a separate class
  would need a conversion. The slots (`@window`, `@crop`, `@post_ops`) are the
  same.
- **Keep WKB and only speed up writing it** (GiottoDisk#110). Fixes write time;
  leaves the size and the terra read cost.
- **GeoParquet 1.1 native point encoding** (a separate x/y struct column). Keeps
  the files GeoParquet, but writes the same coordinates a second time.

## References

- `.writes_wkb()` (`R/methods-storeWrite.R`), `.pgeom_xy_points()` and
  `.pgeom_view_sql()` (`R/methods-storeRead.R`).
- Size bench layout: `id + x_index + y_index [+ geom]`, 5M uniform points,
  200k ids, arrow `write_parquet`.
