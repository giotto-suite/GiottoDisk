#' @include h5ad.R reader-shared.R
NULL

# AnnData .h5ad ingest for `gDirSource`-managed projects.
#
# A reader over the h5ad building blocks in h5ad.R: each `load_*` call is one
# of them with the reader's file and backend filled in, and create_gobject()
# attaches them in the order the other disk readers use (expression, spatial
# locations, metadata, then dimension reductions and networks). No AnnData
# reader exists in memory, so this one has no in-memory parent class.



# CLASS ####

setClass("H5adDiskReader",
    slots = list(
        h5ad_path = "character",
        layer = "ANY",
        calls = "list",
        backend = "ANY"
    ),
    prototype = list(
        layer = NULL,
        backend = NULL
    )
)

# * init ####
setMethod(
    "initialize", signature("H5adDiskReader"),
    function(.Object, ..., h5ad_path, layer, backend) {
        obj <- callNextMethod(.Object, ...)

        if (!missing(h5ad_path)) {
            checkmate::assert_file_exists(h5ad_path)
            obj@h5ad_path <- normalizePath(path.expand(h5ad_path))
        }
        if (!missing(layer)) obj@layer <- layer
        if (!missing(backend)) {
            if (is.character(backend)) backend <- gDirSource(path = backend)
            checkmate::assert_class(backend, "gsource")
            obj@backend <- backend
        }
        if (is.null(obj@backend)) {
            stop("[H5adDiskReader] `backend` is required", call. = FALSE)
        }
        if (!length(obj@h5ad_path)) return(obj)

        h5ad <- obj@h5ad_path
        gsrc <- obj@backend
        lyr <- obj@layer

        ## expression (disk) ####
        obj@calls$load_expression <- function(
            path = h5ad,
            layer = lyr,
            name = "raw",
            spat_unit = "cell",
            feat_type = "rna",
            batch_records = NULL,
            output = c("exprObj", "store"),
            verbose = NULL
        ) {
            output <- match.arg(output)
            GiottoUtils::vmsg(.v = verbose, "[H5adDiskReader] expression <-",
                if (is.null(layer)) "X" else layer)
            pe <- sourceWrite(gsrc,
                h5adInput(path, layer = layer, batch_records = batch_records),
                store_type = "parquetExpr")
            if (output == "store") return(pe)
            GiottoClass::createExprObj(pe, name = name, spat_unit = spat_unit,
                feat_type = feat_type, provenance = spat_unit)
        }

        ## other slots (in memory, or edges streamed to the vault) ####
        obj@calls$load_spatlocs <- function(path = h5ad, key = "spatial",
                                            name = "raw", flip_y = FALSE,
                                            spat_unit = "cell") {
            h5adSpatLocs(path, key = key, name = name, flip_y = flip_y,
                         spat_unit = spat_unit)
        }
        obj@calls$load_cellmeta <- function(path = h5ad, columns = NULL,
                                            spat_unit = "cell",
                                            feat_type = "rna") {
            h5adCellMeta(path, columns = columns, spat_unit = spat_unit,
                         feat_type = feat_type)
        }
        obj@calls$load_featmeta <- function(path = h5ad, layer = lyr,
                                            columns = NULL,
                                            spat_unit = "cell",
                                            feat_type = "rna") {
            h5adFeatMeta(path, layer = layer, columns = columns,
                         spat_unit = spat_unit, feat_type = feat_type)
        }
        obj@calls$load_dimreds <- function(path = h5ad, keys = NULL,
                                           spat_unit = "cell",
                                           feat_type = "rna",
                                           verbose = NULL) {
            h5adDimReds(path, keys = keys, spat_unit = spat_unit,
                        feat_type = feat_type, verbose = isTRUE(verbose))
        }
        obj@calls$load_networks <- function(path = h5ad, keys = NULL,
                                            spat_unit = "cell",
                                            feat_type = "rna",
                                            batch_records = NULL,
                                            verbose = NULL) {
            h5adNetworks(path, gsrc, keys = keys, spat_unit = spat_unit,
                         feat_type = feat_type, batch_records = batch_records,
                         verbose = isTRUE(verbose))
        }

        ## giotto object ####
        obj@calls$create_gobject <- function(
            load_expression = TRUE,
            load_spatlocs = TRUE,
            load_cellmeta = TRUE,
            load_featmeta = TRUE,
            load_dimreds = TRUE,
            load_networks = TRUE,
            layer = lyr,
            spatial_key = "spatial",
            flip_y = FALSE,
            spat_unit = "cell",
            feat_type = "rna",
            gobject = NULL,
            instructions = NULL,
            verbose = NULL
        ) {
            funs <- obj@calls
            g <- .backed_gobject_init(gobject, gsrc, instructions)
            set <- function(g, x) GiottoClass::setGiotto(g, x, verbose = FALSE)

            if (load_expression) {
                g <- set(g, funs$load_expression(layer = layer,
                    spat_unit = spat_unit, feat_type = feat_type,
                    verbose = verbose))
            }
            if (load_spatlocs) {
                g <- set(g, funs$load_spatlocs(key = spatial_key,
                    flip_y = flip_y, spat_unit = spat_unit))
            }
            if (load_cellmeta) {
                g <- set(g, funs$load_cellmeta(spat_unit = spat_unit,
                    feat_type = feat_type))
            }
            if (load_featmeta) {
                g <- set(g, funs$load_featmeta(layer = layer,
                    spat_unit = spat_unit, feat_type = feat_type))
            }
            if (load_dimreds) {
                for (x in funs$load_dimreds(spat_unit = spat_unit,
                    feat_type = feat_type, verbose = verbose)) g <- set(g, x)
            }
            if (load_networks) {
                for (x in funs$load_networks(spat_unit = spat_unit,
                    feat_type = feat_type, verbose = verbose)) g <- set(g, x)
            }
            g
        }

        obj
    }
)



# access ####

#' @export
setMethod("$", signature("H5adDiskReader"), function(x, name) {
    if (name %in% c("h5ad_path", "layer")) {
        return(methods::slot(x, name))
    }
    x@calls[[name]]
})

#' @export
setMethod("$<-", signature("H5adDiskReader"), function(x, name, value) {
    if (name %in% c("h5ad_path", "layer")) {
        methods::slot(x, name) <- value
        return(initialize(x))
    }
    stop("Only 'h5ad_path' and 'layer' can be set", call. = FALSE)
})

#' @export
`.DollarNames.H5adDiskReader` <- function(x, pattern) {
    dn <- c("h5ad_path", "layer")
    if (length(x@calls)) dn <- c(dn, paste0(names(x@calls), "()"))
    utils::.DollarNames(dn, pattern)
}

setMethod("show", signature("H5adDiskReader"), function(object) {
    cat(sprintf("GiottoDisk <%s>\n", class(object)[1L]))
    f <- object@h5ad_path
    cat("file  :", if (length(f)) GiottoUtils::str_abbreviate(f) else "", "\n")
    cat("layer :", if (is.null(object@layer)) "X" else object@layer, "\n")
    if (length(object@calls)) {
        cat("funs  :", paste0(names(object@calls), "()", collapse = ", "),
            "\n")
    }
})



# CREATE READER ####

#' @title Import an AnnData .h5ad file (disk-backed)
#' @name importH5adDisk
#' @description
#' Disk-backed reader for AnnData `.h5ad` files, without python. The
#' expression matrix (`X`, `raw/X` or a layer, sparse CSR or CSC) streams into
#' a `parquetExprStore` in the `gDirSource`-managed project vault, and the
#' `obsp` graphs named in `uns` stream into `parquetEdgeStore`s. Cell and
#' feature metadata, spatial locations and dimension reductions are read into
#' memory. Each `load_*()` function is one of the h5ad building blocks
#' ([h5adInput()], [h5adSpatLocs()], [h5adCellMeta()], [h5adFeatMeta()],
#' [h5adDimReds()], [h5adNetworks()]); `create_gobject()` attaches them all,
#' optionally into an existing backed `gobject`.
#' @param h5ad_path path to the `.h5ad` file
#' @param backend a `gsource` (typically `gDirSource`) project backend, or a
#'   path to one
#' @param layer character or NULL. Expression matrix to read: `NULL` or `"X"`
#'   for `X`, `"raw"` for `raw/X`, otherwise a name under `layers/`. Dense
#'   matrices are refused; many files keep their counts as a sparse layer.
#' @returns `H5adDiskReader` object
#' @family h5ad readers
#' @examples
#' \dontrun{
#' rdr <- importH5adDisk("data.h5ad", backend = "my_project", layer = "counts")
#' g <- rdr$create_gobject()
#' }
#' @export
importH5adDisk <- function(h5ad_path = NULL, backend, layer = NULL) {
    if (missing(backend)) {
        stop("[importH5adDisk] `backend` is required", call. = FALSE)
    }
    a <- list(Class = "H5adDiskReader", backend = backend, layer = layer)
    if (!is.null(h5ad_path)) a$h5ad_path <- h5ad_path
    do.call(methods::new, args = a)
}
