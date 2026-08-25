# Giotto Directory Source

Extends
[gsource](https://giotto-suite.github.io/GiottoDisk/reference/gsource.md).
Used to designate and manage a directory as a giotto project directory
containing on-disk data artifacts. Contained artifacts are tracked via a
manifest with editable metadata.

Edits to the `gDirSource` manifest are atomic. Additionally, a WAL
pattern is used to ensure concurrency safety. Edits are only made as
loose files specific to single artifacts (and by extension usually
single-process). Loose edits are later consolidated into the main
manifest on-read (preferably only from the main process). filelock
integration is available for extra concurrency safety to avoid edge
cases where two processes attempt consolidation at the same time but see
different `_pending` states.

## Slots

- `path`:

  character. Directory to use

## Artifact tracking

Data can be written to the backing directory using
[`sourceWrite()`](https://giotto-suite.github.io/GiottoDisk/reference/sourceWrite.md).
This serializes the data to disk using the preferred store types for the
`gsource` class. The artifact is also assigned an
[artifact_uid](https://giotto-suite.github.io/GiottoDisk/reference/artifact_uid.md)
at this point.

## Directory Structure

- `giottodir.json` - json file manifest of the contained artifacts.

- `giottodir.json.lock` - (optional) lockfile for filelock to control
  writes to `giottodir.json` for extra concurrency safety

- `artifacts` - vault directory containing the actual data artifacts
  which are within subdirectories named by their `uid`.

- `_pending` - directory containing pending edits to the manifest.

- `giottosave` - directory containing .RDS of giotto projects that
  reference assets controlled by the `gDirSource`.
