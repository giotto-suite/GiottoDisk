# Snapshots live in giottosave/<name>.<id>/ and are tagged `<name>.<id>`;
# tests address them by name.

.snap_gsdir <- function(proj) file.path(proj, "giottosave")

# tag (folder name) of the current snapshot called `name`
.snap_tag <- function(proj, name) {
    .gdsrc_find_snapshot(.snap_gsdir(proj), name = name)$tag
}

# object file of the current snapshot called `name`, NULL if there is none
.snap_file <- function(proj, name) {
    .gdsrc_find_snapshot(.snap_gsdir(proj), name = name)$file
}

# write `obj` as a snapshot called `name` without going through snapshotSave
.write_snapshot_fixture <- function(proj, name, obj) {
    d <- file.path(.snap_gsdir(proj), .gdsrc_snapshot_tag(name))
    dir.create(d, recursive = TRUE)
    f <- file.path(d, "gobject.RDS")
    saveRDS(obj, f)
    invisible(f)
}
