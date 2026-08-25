# gDirSource: the directory-backed Giotto source

**`gsource`** is the abstract Giotto backend type — the contract a
Giotto object’s stores can be attached to for managed, on-disk
durability. The verbs `sourceWrite`, `sourceContains`, `sourceAdopt`,
`sourcePrune`, and the `snapshot*` family all dispatch on this type, so
new backends can plug in without changing call sites.

The only `gsource` implementation currently shipped is **`gDirSource`**
— a managed project directory on the local filesystem. It tracks
artifacts in a JSON manifest, provides a vault layout for storing files,
and lets you snapshot a Giotto object to durable storage with the
artifacts it depends on intact. Future backends (e.g. a SpatialData
adapter, a cloud-object-store backend) would implement the same verbs
against different storage substrates; see
`vignettes/articles/roadmap.Rmd` for the public direction.

This vignette is about `gDirSource` specifically. It walks through:

- the **project directory layout** a `gDirSource` manages
- the **`source*` verbs** for writing, locating, adopting, and pruning
  artifacts (`sourceWrite`, `sourceContains`, `sourceAdopt`,
  `sourcePrune`)
- the **`snapshot*` verbs** for object-level lifecycle (`snapshotSave`,
  `snapshotLoad`, `snapshotDelete`)
- common **deployment patterns** (local interactive, NFS/HPC, shared
  analyses)

For the *why* — manifest design, dump lifecycle, multi-analysis
guarantees, deployment trade-offs — see `vignettes/articles/design.Rmd`,
“Project Management (gsource)”.

## Quick start

``` r

library(GiottoDisk)

# Create / attach to a project directory
src <- gDirSource("/path/to/my_project")

# Write a store; the source allocates a uid and records the artifact
pts <- terra::vect(data.frame(x = 1:5, y = 1:5), geom = c("x", "y"))
pgs <- parquetGeomStore() |> storeWrite(pts)
pgs <- sourceWrite(src, pgs)   # now lives in the vault

# Snapshot a giotto object — adopts any external stores along the way
gobject@source <- src
snapshotSave(src, gobject, name = "qc_pass")

# Reload in a later session
src2 <- gDirSource("/path/to/my_project")
g <- snapshotLoad(src2, "qc_pass")
```

## Project directory layout

`gDirSource(path)` manages a directory with this structure:

    <project_dir>/
      giottodir.json       # manifest of all tracked artifacts
      artifacts/<uid>/...  # vault — each artifact gets its own uid-named subdir
      _pending/            # write-ahead-log staging area for manifest edits
      giottosave/          # `.rds` / `.qs` snapshots of giotto objects
      dump/                # (optional) persistent staging area

The directory is created on first use. Concurrent writers are safe —
each `sourceWrite` queues an atomic entry in `_pending/` rather than
racing on `giottodir.json` directly; the next read consolidates pending
entries into the manifest.

You can inspect the manifest at any time:

``` r

src["abc123"]              # full record for one uid
src["abc123", "store_type"]  # one field
src[, "store_type"]        # one column across all uids
as.data.frame(src)         # full manifest as a data.frame
```

## Source verbs

### `sourceWrite(src, store)` — register a new artifact

Allocates a fresh uid, writes the store under `artifacts/<uid>/`, hashes
the contents, and queues a pending manifest entry. Subsequent reads
through the source see the artifact.

``` r

pgs <- parquetGeomStore() |> storeWrite(my_points)
pgs <- sourceWrite(src, pgs)   # uid assigned, files moved into vault
```

After `sourceWrite`, `pgs@uid` matches a manifest entry and `pgs@path`
points into the vault.

### `sourceContains(src, store)` — is this artifact managed?

Returns `TRUE` if the store is already vault-resident and registered.
Dispatch varies by store kind:

- `fileStore` (parquetStore, parquetExprStore, etc.): uid lookup in the
  manifest.
- `SpatRaster`: checks whether `terra::sources(store)` paths sit under
  the vault prefix. In-memory rasters always return `FALSE`.
- Union stores: `TRUE` only if **all** substores are contained.
- `IterableMatrix` (BPCells, via the `ANY` fallback): checks all leaf
  `@dir` paths against the vault prefix.

### `sourceAdopt(src, store)` — move/register an unmanaged store

Brings a store under source management. The exact behaviour depends on
where the store currently lives and what type it is.

**`fileStore` (parquet stores)**: preserves the store’s existing `@uid`,
moves files into `artifacts/<uid>/`, writes the manifest entry.

**`SpatRaster`**: allocates a fresh uid. On-disk rasters are
`file.copy`’d into the vault
([`terra::rast()`](https://rspatial.github.io/terra/reference/rast.html)
re-opens at the new path); in-memory rasters are written as COG.

**`IterableMatrix`**: the leaf `@dir` paths are walked via
`.im_map_leaves` and each leaf is moved into the vault. Compound
structures (`@matrix_list` from cbind/rbind, lazy transforms wrapping a
leaf) are rebuilt with the new leaf paths.

**Union stores**: delegate to per-substore adoption.

``` r

# An overlay raster you computed outside the project
r <- terra::rast("/scratch/overlay.tif")
r <- sourceAdopt(src, r)   # now in the vault, registered, safe across sessions
```

#### Vault-resident but unregistered

If a path already sits under the vault (common when using
`setArtifactDumpDir(src)` for HPC, see below), `sourceAdopt` extracts
the uid from the path, checks `src[uid]`, and writes the manifest entry
if missing — no file movement, just registration.

#### Shared leaves (raw / normalized matrices)

`raw` and `normalized` expression matrices often share the same
underlying BPCells leaf directory (`@dir`). Adopting `raw` would move
that directory; adopting `normalized` would then fail with ENOENT trying
to read the old path.

`snapshotSave` mitigates this with a session-scoped path map
(`.adopt_session_map`): the first adoption of a leaf records old → new,
subsequent adoptions of the same leaf redirect to the cached vault path
without attempting a second move. The map resets at the start of each
`snapshotSave`.

#### External-path guard

By default, `sourceAdopt` only moves data from the **dump directory**.
Paths outside the dump are skipped with a warning to avoid silently
moving data the user owns. Override with:

``` r

options(giottodisk.adopt_external = TRUE)
```

### `sourcePrune(src)` — clean up unreferenced artifacts

Removes artifacts not tagged by any `giottosave` snapshot. Protection is
transitive: if snapshot A references store X, and X depends on Y, both X
and Y are kept (BFS over the manifest’s `depends` graph).

``` r

sourcePrune(src)  # safe — only deletes things no snapshot needs
```

## Snapshot verbs

A snapshot is a serialized giotto object (`.rds` or `.qs`) under
`giottosave/<name>.<ext>` plus the manifest tags that pin its artifacts.
Snapshots are how a giotto object survives across sessions without
breaking when the artifacts it references get relocated or pruned.

### `snapshotSave(src, gobject, name)`

Five-step lifecycle, in order:

1.  **Session reset** — `.adopt_session_reset()` clears the per-leaf
    path cache.
2.  **Adopt external images** — any `SpatRaster` slot pointing outside
    the vault is `sourceAdopt`’d.
3.  **Adopt external expression matrices** — `IterableMatrix` /
    `HDF5Array` not in the vault are adopted (sharing the session cache
    so shared leaves move once).
4.  **Write snapshot atomically** — temp file → `file.rename` to
    `giottosave/<name>.rds` (or `.qs`).
5.  **Tag artifacts** — every store referenced by the gobject (or
    reachable via `depends`) gets `name` appended to its `giottosave`
    manifest field. Tagging is hash-based via `.ss_hash_expr_base`,
    which strips lazy `@ops` so a filtered store and its unfiltered
    parent map to the same artifact.

After save, `sourcePrune` is safe — every artifact reachable from any
snapshot stays.

### `snapshotLoad(src, name)`

Returns the giotto object snapshot. Stores inside the object retain
their source backing — subsequent saves and prunes work the same way.

``` r

g <- snapshotLoad(src, "qc_pass")
spatPlot(g)
```

### `snapshotDelete(src, name)`

Removes the `giottosave/<name>` file **and** untags `name` from all
artifacts in the manifest. After delete, artifacts that were tagged
*only* by this snapshot become prune-eligible.

``` r

snapshotDelete(src, "old_attempt")
sourcePrune(src)
```

## Common patterns

### Shared `gDirSource` across analyses

Multiple gobjects can share one `gDirSource`. This is the recommended
pattern when analyses share underlying data (same experiment, different
processing pipelines):

``` r

src <- gDirSource("/projects/visium_exp1")
g1@source <- src   # initial qc pipeline
g2@source <- src   # re-analysis with different params

# Both snapshots can reference the same `raw` matrix in the vault.
snapshotSave(src, g1, name = "qc_v1")
snapshotSave(src, g2, name = "qc_v2")
```

The vault holds the shared artifacts once; `sourcePrune` keeps any
artifact referenced by *at least one* snapshot across all gobjects in
the project.

**Granularity:** per-experiment or per-sample-cohort. Large enough to
get deduplication benefits across related analyses, small enough that
`sourcePrune` keeps the manifest tidy. A single `gDirSource` shared
across an entire lab would accumulate artifacts unboundedly.

### HPC / NFS: persistent dump

On HPC systems, the default
[`tempdir()`](https://rdrr.io/r/base/tempfile.html)-based dump dies with
the session, which is wrong for long-running jobs. Point the dump at the
`gDirSource`’s own vault:

``` r

src <- gDirSource("/projects/visium_exp1")
setArtifactDumpDir(src)   # dump = artifacts/

# Artifacts now write directly to their final vault location.
# sourceAdopt becomes a registration step — no file movement at all.
```

| Setting | Recommended dump |
|----|----|
| Local / workstation | [`tempdir()`](https://rdrr.io/r/base/tempfile.html) (default) |
| NFS / HPC (SLURM) | `setArtifactDumpDir(src)` |
| SLURM array jobs | One `gDirSource` per sample (no shared-manifest contention) |
| Multi-analysis | Single shared `gDirSource` |

### What survives a crash?

- `_pending/<uid>.json` writes are atomic (write to temp → rename). A
  crashed writer leaves nothing behind.
- `giottodir.json` consolidation is atomic (write to temp → rename over
  the old file). Either the old manifest or the new manifest is always
  intact, never a partial write.
- Snapshot writes are atomic (`giottosave/<name>.tmp` → rename to
  `<name>.rds`).
- The vault is append-only via `sourceWrite` / `sourceAdopt`. Nothing
  inside `artifacts/` ever gets modified in place.

If a process crashes mid-`snapshotSave` *after* artifact adoption but
*before* the snapshot file lands, the adopted artifacts are “abandoned”
— present in the vault and the manifest but not yet tagged by any save.
The next `sourcePrune` will reclaim them.

## Further reading

- `vignettes/articles/design.Rmd` — Project Management (gsource):
  manifest design rationale, dump lifecycle, IterableMatrix adoption
  algorithm, multi-analysis architecture, deployment trade-offs.
- [`?gDirSource`](https://giotto-suite.github.io/GiottoDisk/reference/gDirSource.md),
  [`?sourceWrite`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite.md),
  [`?snapshotSave`](https://giotto-suite.github.io/GiottoDisk/reference/snapshotSave.md)
  etc. — generic-level reference docs.
