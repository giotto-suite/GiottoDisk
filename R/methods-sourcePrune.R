#' @name sourcePrune
#' @title Prune a project source
#' @description
#' On-disk processing can generate many intermediate
#' artifacts. Manually call this function in order to
#' drop artifacts that are not tagged with a
#' giottosave that depends on them (the giottosave .rds
#' file must also exist).
#' @param src `gsource` object. Used to provide default save locations
#' and save formats for a managed Giotto backend.
#' @export
setMethod("sourcePrune", signature("gDirSource"),
    function(src, ...) {
        .gdsrc_artifact_prune(src@path)
    })