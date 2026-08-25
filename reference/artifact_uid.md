# Artifact Unique Identifier

Artifacts are tracked in
[gDirSource](https://giotto-suite.github.io/GiottoDisk/reference/gDirSource-class.md)
using an unique identifier. These are randomized alphanumeric strings of
length 8 that can additionally be tagged with the process ID (pid) and
node ID that generated it. Only the pid is tagged by default, but the
node can also be tagged to help prevent collisions on multi- node setups
by setting the appropriate option below:

- `giottodisk.uid_include_pid` (default = TRUE)

- `giottodisk.uid_include_node` (default = FALSE)

These uids are not affected by user seed setting and instead use a
temporary random seed based on the time and a process-specific counter.
