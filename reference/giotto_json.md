# Giotto Directory JSON

Internal utilities for editing the
[gDirSource](https://giotto-suite.github.io/GiottoDisk/reference/gDirSource.md)
managed `giottodir.json`.

## Usage

``` r
.json_atomic_write(
  x,
  file,
  temp_file = tempfile(tmpdir = dirname(file), fileext = ".tmp"),
  cleanup = TRUE
)

.json_read(file)

.gdsrc_json_header()

.gdsrc_json_write(p, x)

.gdsrc_json_read(p, consolidate = FALSE)

.gdsrc_json_edit(p, uid, x)

.gdsrc_json_add_artifact(
  p,
  store_type,
  uid,
  hash,
  meta = NULL,
  giottosave = NA_character_,
  depends = NULL
)

.gdsrc_json_consolidate(p)
```

## Arguments

- uid:

  character. Unique ID for artifact tracking

- store_type:

  character. Type of file format storage (e.g. \\h5\\)

- hash:

  character. Hash of the in-memory store object to track changes

- meta:

  list (optional). Additional list of atomic object(s) that can

- giottosave:

  `character` tagged giotto save(s) if any. be attached as further
  metadata to the particular uid

## Functions

- `.json_atomic_write()`: Safe atomic writes of information to json
  files so that there is no instant at which the data does not exist on
  disk.

- `.json_read()`: json reading with specific params to work with
  GiottoDisk's usecase where values are generally atomic vectors

- `.gdsrc_json_header()`: general header to apply at top of json files
  from GiottoDisk. Stamps with information for:

  - GiottoDisk version

  - manifest version

- `.gdsrc_json_write()`: Write function for the main `giottodir.json`.

- `.gdsrc_json_read()`: Read function for the main `giottodir.json`.
  This function does not check or consolidate pending edits by default.

- `.gdsrc_json_edit()`: Queue an edit to the manifest by generating a
  pending edit that is written to disk as a loose json. This can later
  be consolidated into the main `giottodir.json`.

- `.gdsrc_json_add_artifact()`: Add a tracked artifact to the directory.
  This adds the tracking metadata for the object and places it in a
  location specified by `.gdsrc_json_pending_dir()` to await
  consolidation into the main `giottodir.json`.

- `.gdsrc_json_consolidate()`: Consolidate pending edits into central
  `giottodir.json` manifest. Scans for pending edits then applies the
  changes in a for loop on the manifest content before writing back out.
