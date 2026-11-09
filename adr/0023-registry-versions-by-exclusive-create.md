# 0023. The registry is committed as numbered versions by exclusive create; edits only create records or add tags

- **Status:** Proposed
- **Date:** 2026-11-09
- **Supersedes:** —
- **Superseded by:** —

## Context

`giottodir.json` was one file, rewritten whole at each consolidation under an
advisory `filelock`. Pending edits were `_pending/<uid>.json`, merged into a
record field by field with `modifyList()`. Snapshot tagging read an
artifact's tag list, appended its tag and wrote the whole list back.

Two writers saving into one project lost tags: each wrote a list built from
its own, possibly stale, read, and `modifyList()` replaces a vector instead of
merging it. With `filelock` installed, two processes adding 20 tags each to
one artifact kept 37, 40 and 22 of 40 in three runs. A tag lost from an
artifact that only that snapshot protects lets `sourcePrune()` delete data the
snapshot still references. The lock could not help: it serialized
consolidations, not the stale reads. Besides, `sourcePrune()` rewrote the
registry without taking the lock, a lock timeout carried on unlocked, and
`filelock` was only suggested. The viewer's working states, jobs and second
sessions will save into shared projects far more often than before.

## Decision

Edits carry one of two operations: `create` (a record's fields, written once;
on an existing record only missing fields are filled) and `add` (values
joined into the set fields `giottosave` and `depends`). Each edit is its own
uniquely named file. Edits then give the same result in any order and when
applied twice. No field is changed after creation; `src[uid, field] <-`
errors for anything but `giottosave`, which adds.

The registry is a series of whole files, `giottodir.<n>.json`. A commit
hard-links a complete temp file to the next number, which fails if that name
exists, so of writers that read the same version one commits and the others
re-read and retry. Readers take the highest number whose header (`seq`)
agrees with its name, falling back to the one below. Each version records the
edit files it applied; those are deleted one version later, and two versions
are kept. `giottodir.json` becomes a stub stamped with schema v3, so an older
GiottoDisk refuses the project. `filelock` is gone. `snapshotSave()` queues its
tags before renaming its folder in, and `sourcePrune()` commits the same way,
deletes vault folders only after its commit, and drops tags whose folder is
gone.

## Consequences

- No lock: correctness rests on create-if-absent, which hard links give on
  local filesystems and NFS. Without hard links (FAT, some SMB shares, synced
  folders) commits fail with an explanation;
  `options(giottodisk.commit_by_rename = TRUE)` commits by rename, safe only
  for a single writer.
- A remote backend needs the same two primitives, create-if-absent (S3
  `If-None-Match: *`, GCS generation 0) and a listing, in place of
  `.gdsrc_publish()` and `.gdsrc_registry_versions()`; the format is unchanged.
- `sourcePrune()` is still not safe against a save running alongside it: a
  folder renamed in after pruning listed the folders can lose artifacts only
  it protects, and dropping dead tags makes that worse. This is left to
  discipline (do not prune while saving), not closed with grace periods.
- A field cannot change after creation. A feature that needs to update one
  must add last-writer-wins with a stored time per field, and bump the schema.
- Old-format pending edits are refused after the upgrade; `sourceUpdate()`
  folds them in first.

## Alternatives considered

- **Keep the lock, make tags additive:** fixes the tag race but not two
  consolidations dropping each other's edits, nor unlocked pruning, and locks
  are no-ops on some network filesystems.
- **Fixed file name with a check before rename:** the check and the rename are
  separate steps, so a writer can still overwrite another's commit. Safe only
  with a conditional write, which local filesystems do not offer.
- **Delta-style log of changes with checkpoints:** more files, replay and
  checkpoint policy for no gain while full versions stay small.
- **SQLite, DuckDB, DuckLake, Iceberg/Delta, a Rust `object_store` commit
  layer:** each gave up something major (network filesystems, a single
  writer, tables only, or Rust on the core path).
- **Per-save protection entries in the registry:** write cost was not the
  problem; folders stay the record of which saves exist.

## References

`R/class-gsource.R` (registry versions, edits, pruning),
`R/methods-sourceSchema.R` (`.gdsrc_upgrade_2_to_3`),
`R/methods-snapshotSave.R`; `tests/testthat/test-registry-commits.R`;
design.Rmd "Manifest design".
