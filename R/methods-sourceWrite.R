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
#' * *matrices:* `"giotto.gdsrc_matrix_format"` (default = "h5")
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
    function(src, data, meta = NULL, ...) {
        store_type <- getOption("giotto.gdsrc_matrix_format", "h5")
        .gdsrc_write_artifact(src,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })

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
    function(src, data, meta = NULL, ...) {
        store_type <- getOption("giotto.gdsrc_spatvector_format", "parquetGeom")
        .gdsrc_write_artifact(src,
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
    function(src, data, meta = NULL, ...) {
        store_type <- getOption("giotto.gdsrc_dataframe_format", "parquet")
        .gdsrc_write_artifact(src,
            data = data,
            store_type = store_type,
            meta = meta,
            ...
        )
    })


# internals ####

.gdsrc_write_artifact <- function(gsrc, data, store_type, meta = NULL, ...) {
    # ensure store_type exists in catalog
    gsrc <- .gdsrc_json_add_store(
        gsrc = gsrc,
        store_type = store_type,
        path = store_type
    )

    # setup save info
    uid <- .make_uid()
    dirpath <- gsrc[store_type]
    savepath <- file.path(dirpath, uid)
    if (!dir.exists(dirpath)) dir.create(dirpath, recursive = TRUE)
    # write to disk
    dstore <- storeCreate(type = store_type, path = savepath)
    dstore <- storeWrite(dstore, data, ...)
    # record hash of delayed representation
    hash <- .hash(storeRead(dstore))
    # record artifact entry
    .gdsrc_json_add_artifact(
        gsrc = gsrc,
        store_type = store_type,
        uid = uid,
        hash = hash,
        meta = meta
    )
    dstore
}
