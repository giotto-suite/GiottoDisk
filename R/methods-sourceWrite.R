# docs ####

#' @name sourceWrite
#' @title Write Data to a Source
#' @description
#' Write data to a [gsource] that manages data formats and file
#' organization. Docs for specific `gsource` implementations have further
#' information.
#'
#' * [sourceWrite-gDirSource]
#'
#' @param src `gsource` object. Used to provide default save locations
#' and save formats for a managed Giotto backend.
#' @param data data to write
#' @param meta `list`. Additional metadata to attach to this Giotto backend
#' @param store_type `character`. Store type to write as
#' @param ... additional params passed to `storeWrite` method
#' managed artifact.
#' @returns written `store` object
#' @family sourceWrite methods
NULL

#' @name sourceWrite-gDirSource
#' @title Write Data to a `gDirSource`
#' @inheritParams sourceWrite
#' @description
#' Write data to a [gDirSource] managed project directory. Defaults for
#' different data types are settable via global options:
#' 
#' * *sparse matrices:* `giotto.gdsrc_sparsematrix_format` (default = "bpcells")
#' * *dense matrices:* `"giotto.gdsrc_densematrix_format"` (default = "h5")
#' * *dataframes:* `"giotto.gdsrc_dataframe_format"` (default = "parquet")
#' * *spatvector:* `"giotto.gdsrc_spatvector_format"` (default = "parquetGeom")
#'
#' @examples
#' testdir <- file.path(tempdir(), "testdirsource")
#' gsrc <- sourceCreate(path = testdir)
#' @returns written `store` object
#' @family sourceWrite methods
NULL

# definitions ####

#' @rdname sourceWrite-gDirSource
#' @examples
#' # arrays ----------------------------------------- #
#' m <- matrix(1:9, nrow = 3)
#' m_written <- sourceWrite(gsrc, m)
#' storeRead(m_written)
#' @export
setMethod("sourceWrite", signature("gDirSource", "memoryMatrix"),
    function(src, data, meta = NULL, store_type = NULL, ...) {
        if (is.null(store_type)) {
            store_type <- if (.is_sparse_matrix(data)) {
                getOption("giotto.gdsrc_sparsematrix_format", "bpcells")
            } else {
                getOption("giotto.gdsrc_densematrix_format", "h5")
            }
        }

        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })

.is_sparse_matrix <- function(mat) {
    if (inherits(mat, "sparseMatrix")) return(TRUE)
    if (inherits(mat, "matrix")) return(sum(mat != 0) / length(mat) < 0.1)
    FALSE
}

#' @rdname sourceWrite-gDirSource
#' @examples
#' # spatvector ----------------------------------------- #
#' sv <- terra::vect(system.file("ex/lux.shp", package="terra"))
#' sv_written <- sourceWrite(gsrc, sv)
#' storeRead(sv_written) # default is output = "query"
#' storeRead(sv_written, output = "sf")
#' storeRead(sv_written, output = "terra")
#' storeRead(sv_written, output = "tibble")
#' @export
setMethod("sourceWrite", signature("gDirSource", "SpatVector"),
    function(src, data, meta = NULL, store_type = getOption("giotto.gdsrc_spatvector_format", "parquetGeom"), ...) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })

#' @rdname sourceWrite-gDirSource
#' @examples
#' # tables ----------------------------------------- #
#' store <- sourceWrite(gsrc, iris)
#' force(store)
#' storeRead(store) # default output is an arrow query
#' storeRead(store, output = "tibble") # pull into memory as tibble
#' @export
setMethod("sourceWrite", signature("gDirSource", "data.frame"),
    function(src, data, meta = NULL, store_type = getOption("giotto.gdsrc_dataframe_format", "parquet"), ...) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })

setMethod("sourceWrite", signature("gDirSource", "fileStore"),
    function(src, data, meta = NULL, store_type, ...) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })


# internals ####

# return allocated artifact save path.
# value is named with the uid
.gdsrc_allocate_artifact_dir <- function(p, uid = NULL, basename = "data", create = TRUE) {
    checkmate::assert_character(p)
    checkmate::assert_character(uid, null.ok = TRUE)
    # setup save info
    if (is.null(uid)) {
        uid <- .make_uid()
    }
    savedir <- .gdsrc_artifact_dir(p, uid = uid)
    if (create && !dir.exists(savedir)) dir.create(savedir, recursive = TRUE)
    savepath <- normalizePath(file.path(savedir, basename), mustWork = FALSE)
    savepath
}

.gdsrc_write_artifact <- function(p, data, store_type, meta = NULL, ...) {
    dstore <- storeCreate(type = store_type)
    uid <- dstore@uid
    savepath <- .gdsrc_allocate_artifact_dir(p, uid = uid, create = TRUE)
    dstore@path <- savepath # update save path
    dstore <- storeWrite(dstore, data, ...)
    # record hash of delayed representation
    hash <- .hash(storeRead(dstore))
    # record artifact entry
    .gdsrc_json_add_artifact(p,
        store_type = store_type,
        uid = uid,
        hash = hash,
        meta = meta
    )
    dstore
}
