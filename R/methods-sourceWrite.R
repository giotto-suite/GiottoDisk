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
#' * *sparse matrices:* `giotto.gdsrc_sparsematrix_format` (default = "parquetExpr")
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
                getOption("giotto.gdsrc_sparsematrix_format", "parquetExpr")
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

#' @rdname sourceWrite-gDirSource
#' @export
setMethod("sourceWrite", signature("gDirSource", "ANY"),
    function(src, data, meta = NULL, ...) {
    if (inherits(data, "IterableMatrix")) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = "bpcells",
            meta = meta,
            ...
        )
    } else {
        stop("[sourceWrite] no method for class: ", class(data), call. = FALSE)
    }
})

#' @rdname sourceWrite-gDirSource
#' @export
# In-memory network → parquetEdgeStore. Used by GiottoClass's network
# setters when a giotto object has a gsource backend attached.
setMethod("sourceWrite", signature("gDirSource", "igraph"),
    function(src, data, meta = NULL, ...) {
        .gdsrc_write_artifact(src@path,
            data = data,
            store_type = "parquetEdgeStore",
            meta = meta,
            ...
        )
    })

#' @rdname sourceWrite-gDirSource
#' @export
#' @description
#' giotto dispatch: promotes an in-memory `giotto` object to a
#' disk-backed one. Builds a fresh `giotto` with `src` attached as the
#' backend, then re-attaches every data subobject via [setGiotto()].
#' Each backend-aware setter (setExpression / setPolygonInfo /
#' setFeatureInfo / setNearestNetwork / setSpatialNetwork) routes its
#' in-memory subobject to the corresponding disk-backed store
#' (parquetExprStore / parquetGeomStore / parquetEdgeStore).
#'
#' Subobjects whose setters don't (yet) have backend-aware write
#' logic stay in-memory but are still attached to the new gobject.
#' @return `giotto` with `@source = src` and backable subobjects
#'   disk-backed.
setMethod("sourceWrite", signature("gDirSource", "giotto"),
    function(src, data, ...) {
        # Cross-source warning: if `data` is already backed by a
        # different gsource, the rebuild flow will re-attach existing
        # disk-backed subobjects to `src` (because each setter's
        # `inherits(dataStore)` guard skips re-writing them), but those
        # references still point at the old vault.
        if (!is.null(data@source)) {
            same_src <- tryCatch(
                identical(
                    normalizePath(data@source@path, mustWork = FALSE),
                    normalizePath(src@path, mustWork = FALSE)
                ),
                error = function(e) FALSE
            )
            if (!isTRUE(same_src)) {
                warning(
                    "[sourceWrite(gDirSource, giotto)] data is already ",
                    "backed by a different source. Disk-backed ",
                    "subobjects will be re-attached without moving ",
                    "artifacts — references may be stale.\n",
                    "  data@source: ", data@source@path, "\n",
                    "  src:          ", src@path,
                    call. = FALSE
                )
            }
        }

        # Idempotency is handled at the subobject level: every
        # backend-aware setter checks `inherits(x@<slot>, "dataStore")`
        # and skips re-writing if already disk-backed. So a fully-backed
        # input gobject round-trips unchanged. If a subobject snuck
        # through in-memory (e.g. because some setter didn't yet have
        # auto-write at the time), it gets disk-backed on this pass.

        # Flat list of every data subobject. Mirrors sliceGiotto's
        # pattern — `[[slot_names]]` returns a list of subobjects from
        # the requested giotto slots.
        dataslots <- c(
            "spatial_info", "spatial_locs", "spatial_network",
            "feat_info", "expression", "cell_metadata",
            "feat_metadata", "spatial_enrichment", "nn_network",
            "dimension_reduction", "multiomics"
        )
        datalist <- data[[dataslots]]

        # Fresh gobject — copy static slots over, attach src as backend.
        # Avoids @h5_file (deprecated) and @misc (transient).
        g_new <- GiottoClass::giotto(
            expression_feat = data@expression_feat,
            images          = data@images,
            parameters      = data@parameters,
            instructions    = data@instructions,
            offset_file     = data@offset_file,
            versions        = data@versions,
            join_info       = data@join_info,
            source          = src
        )

        # Re-attach all subobjects. Backend-aware setters (expression,
        # polygon, points, NN net, spatial net) auto-write through to
        # the gsource'd vault on the way in.
        g_new <- GiottoClass::setGiotto(g_new, datalist,
            initialize = FALSE, verbose = FALSE)

        methods::initialize(g_new)
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
    stats::setNames(savepath, uid)
}

.gdsrc_write_artifact <- function(p, data, store_type, meta = NULL, ...) {
    dstore <- storeCreate(type = store_type)
    uid <- dstore@uid
    savepath <- .gdsrc_allocate_artifact_dir(p, uid = uid, create = TRUE)
    dstore@path <- savepath # update save path
    dstore <- storeWrite(dstore, data, ...)
    # record hash of delayed representation
    hash <- .hash(storeRead(.store_nostate(dstore)))
    # record artifact entry
    .gdsrc_json_add_artifact(p,
        store_type = store_type,
        uid = uid,
        hash = hash,
        meta = meta
    )
    dstore
}
