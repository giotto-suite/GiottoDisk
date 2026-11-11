# 0025. A prune keeps what recently deleted snapshots reached

- **Status:** Proposed
- **Date:** 2026-11-11
- **Supersedes:** —
- **Superseded by:** —

## Context

Several processes writing to one project is intended: a console session, a
viewer, workers running jobs. `sourcePrune()` deletes every artifact no live
snapshot reaches, so a prune in one process can remove data another is
using:

- a save that replaces a snapshot of the same name deletes the older one,
  which another session may have loaded;
- a session's cleanup deletes a snapshot another session has just loaded.

In both cases the data belonged to a snapshot that existed moments before.
Artifacts that no snapshot ever reached are different: without a snapshot
there is no object that refers to them, so they have little use, and work
worth keeping is saved as a snapshot.

## Decision

A deleted snapshot records when it was deleted (`deleted_at` in
`snapshot.json`), written before its object file is removed. `sourcePrune()`
treats a snapshot deleted within a grace period as live: what it reaches is
kept, and its tags stay on the records, so the next prune keeps it too. The
period is on by default, 24 hours, set by `sourcePrune(grace =)` or the
`giottodisk.grace` option; `grace = FALSE` turns it off. A snapshot deleted
with no recorded time, by a GiottoDisk from before this, gets no grace
period. Artifacts no snapshot reached get none either.

`sourcePrune(dry_run = TRUE)` reports what a prune would remove and what the
grace period keeps, by which snapshot, and changes nothing. A real prune
returns the same table.

## Consequences

- Space from a deleted snapshot is freed by the first prune more than a day
  after the delete, or at once with `grace = FALSE`.
- Deleting by name stays immediate: the snapshot is gone from listings and
  loads at once; only its data lingers.
- A process that keeps work only as unsaved artifacts can still lose it to
  another process's prune. Processes that run jobs save the results as
  snapshots when they finish.
- Bulk cleanup of snapshots (a session's unkept saves) is to skip snapshots
  created within the same period, using `created`, when it is added.

## Alternatives considered

- **A grace period on artifacts, by creation time:** protects untagged work,
  but needs a time and session on every artifact record and keeps data that
  nothing refers to, while leaving a replaced snapshot's old data
  unprotected.
- **Leases or markers for data in use:** every process would have to renew
  them, and a crashed one would hold data until they expire.
- **No grace period, prune only by discipline:** the current rule, which
  several writers cannot follow.

## References

`R/class-gsource.R` (`.gdsrc_artifact_prune`, `.gdsrc_prune_plan`,
`.gdsrc_graced_tags`), `R/methods-snapshotSidecars.R`
(`.gdsrc_info_mark_deleted`), `R/methods-sourcePrune.R`;
`tests/testthat/test-prune-grace.R`; GiottoDisk#159.
