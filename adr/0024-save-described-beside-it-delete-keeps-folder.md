# 0024. A save is described in a file beside it; deleting one keeps its folder

- **Status:** Proposed
- **Date:** 2026-11-09
- **Supersedes:** —
- **Superseded by:** —

## Context

A snapshot folder held the object and two files derived from it,
`manifest.json` (what it is) and `history.ndjson` (why it is). Nothing recorded
facts about the *save*: who made it, from which snapshot, whether it should
last. GiottoLens is to save working states it may clean up at the end of a
session, highlighted apart from the user's own saves, with the user deciding
what to keep; it also wants a tree of states. Deleting a snapshot removed its
folder, so a child's link to it would dangle unless every child were rewritten.

## Decision

Each snapshot folder holds `snapshot.json`, written in the staging folder with
the object: fixed at save (`schema`, `name`, `tag`, `created`, `creator`,
`session`, `object`, `parent`) and changeable with `snapshotMetadata()` (`keep`,
`description`, `labels`, `extra`). `snapshotInfo()` lists them without loading
anything. Deleting a snapshot, or dropping an overwritten one, removes only its
object file; the folder and its description files stay. A folder without an
object file was already not a snapshot, so protection is released exactly as
before.

## Consequences

- Children always name an existing parent; nothing is repointed on delete, and
  a deleted step can still be shown, or diffed through its `manifest.json`.
- Deleted folders accumulate, a few kilobytes each. Removing those no live
  snapshot descends from can come later if it matters.
- `creator` and `keep` are separate: who made a save never changes, whether it
  lasts can. `session` is a plain id, never a viewer's access token.
- Facts about saves stay out of the registry (adr/0023): the folder is the one
  record of a save, and editing a description never contends with artifact
  writes. Concurrent edits of one save's `snapshot.json` resolve as last write
  wins.
- `parent` is passed by the caller. Recording it automatically on
  `snapshotLoad()` needs a field on the object (GiottoClass) and is left for
  then.

## Alternatives considered

- **Keep save facts in the object** (`@versions`, the manifest or history):
  labelling would change the object and its fingerprints, and a working state
  loaded and saved again as kept would still say "working".
- **A section of the registry per save:** a second record of which saves
  exist, able to disagree with the folders, and more writes to the shared file.
- **Delete the folder and repoint children:** every delete rewrites other
  saves' files.
- **Store the whole ancestor chain:** no repointing either, but a deleted step
  becomes a gap instead of a node that can still be shown.

## References

`R/methods-snapshotSidecars.R` (`snapshotInfo`, `snapshotMetadata`,
`.ss_gdsrc_info`), `R/methods-snapshotSave.R`
(`.gdsrc_snapshot_remove_object`), `R/methods-snapshotDelete.R`;
`tests/testthat/test-snapshot-info.R`; design.Rmd "Snapshots and tags".
