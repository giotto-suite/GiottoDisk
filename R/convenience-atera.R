# Disk-backed Atera reader.
#
# Atera output currently has the same layout as Xenium output, so
# `AteraDiskReader` subclasses `XeniumDiskReader` and overrides nothing. It
# exists so Atera has its own entry point and its own place to diverge later.
#
# Because it subclasses rather than copies, it inherits the `create_gobject`
# orchestration that `XeniumDiskReader` already defines -- so this adds no
# further copy of that duplicated block (see convenience-xenium.R).
#
# `@include` is required: Collate is alphabetical, so without it this file is
# sourced before `convenience-xenium.R` and `contains =` fails at build time.

#' @include convenience-xenium.R
NULL

#' @title Atera disk reader
#' @name AteraDiskReader-class
#' @description
#' Disk-backed reader for Atera output. Inherits `XeniumDiskReader`, since the
#' two output layouts are presently identical.
#' @keywords internal
setClass("AteraDiskReader",
    contains = "XeniumDiskReader",
    # only override: path-detection messages should name Atera
    prototype = list(platform = "Atera")
)

#' @title Import an Atera assay (disk-backed)
#' @name importAteraDisk
#' @description
#' Disk-backed counterpart to [Giotto::importAtera()]. Produces an
#' `AteraDiskReader` which writes transcripts and boundaries to a
#' `gDirSource`-managed project vault, exactly as the Xenium disk reader
#' does -- including zarr-only output directories (the only format Atera
#' will ship); see [importXeniumDisk()] for how zarr input is handled.
#' @param atera_dir Atera output directory
#' @param backend a `gsource` (typically `gDirSource`) project backend.
#' @param qv_threshold minimum Phred-scaled quality score retained. Only
#'   applies when transcript-level data is present.
#' @returns `AteraDiskReader` object
#' @seealso [Giotto::importAtera()] for the in-memory variant
#' @export
importAteraDisk <- function(atera_dir = NULL, backend, qv_threshold = 20) {
    if (missing(backend)) {
        stop("[importAteraDisk] `backend` is required", call. = FALSE)
    }
    a <- list(Class = "AteraDiskReader", backend = backend, qv = qv_threshold)
    # the inherited slot is named `xenium_dir`; only the public argument differs
    if (!is.null(atera_dir)) a$xenium_dir <- atera_dir
    do.call(new, args = a)
}
