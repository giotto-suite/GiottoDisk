# tilework integration to use fileStore classes
#' @include utils-spatial.R
NULL

# methods ####

#' @name getBoundedData
#' @title Get Bounded Data
#' @description
#' Get the region of the data defined by a bounding `SpatExtent`.
#' `crop()` operations on `fileStore` objects return `dataStore`,
#' but this function directly produces a usable interface.
#' @param x `queryableStore` inheriting object
#' @param bound `SpatExtent`
#' @param sdimx `character`. Which column contains spatial x values
#' to filter on
#' @param sdimy `character`. Which column contains spatial y values
#' to filter on.
#' @param group_col `character`. Only used when `envelope = TRUE`.
#' Which column contains grouping IDs that are used to calculate the
#' spatial envelope.
#' @param envelope `logical` (default = FALSE). Whether to perform 
#' bound filtering based on the envelope centroid instead of 
#' individual x and y values so that grouped rows are not
#' separated by the filtering.
#' 
#' This is generally not needed for `parquetGeomStore` inheriting
#' objects since each row is its own geometry (point or polygon, etc)
#' and xy filtering is performed on the point or the centroid.
#' @param inclusive `logical` (length 1 or 4. Default = TRUE).
#' Whether bounds should be inclusive (inclusive is >=/<= vs >/<). 
#' Ordering of inputs is bottom, left, top, right. If a single logical
#' is provided, it will be replicated and used for all 4 bounds.
#' Default is fully inclusive.
#' @param output `character`. Output format passed to [storeRead()].
#' @param ... additional params to pass to [storeRead()]
#' @family tile methods
#' @returns output as determined by `output` param -- see [storeRead()]
#' @export
setMethod("getBoundedData", signature("queryableStore", "SpatExtent"),
    function(x, bound, 
        sdimx = "x",
        sdimy = "y",
        inclusive = TRUE,
        group_col = "id",
        envelope = FALSE,
        output = "query",
        ...) {
    crop_callback <- function(atab) {
        if (!all(c(sdimx, sdimy) %in% names(atab))) {
            stop("[GiottoDisk][getBoundedData] ",
                "sdimx and sdimy cols not found in data\n", 
                call. = FALSE)
        }
      
        if (isTRUE(envelope)) {
            checkmate::assert_character(group_col)
            use_data <- .dplyr_xy_envelopes(atab,
                sdimx = sdimx,
                sdimy = sdimy,
                group_col = group_col
            )
            res <- .dplyr_crop(use_data,
                sdimx = "ecentroid_x",
                sdimy = "ecentroid_y",
                extent = bound,
                inclusive = inclusive
            )
            dplyr::semi_join(atab, res,
                by = stats::setNames("id", group_col)
            )
        } else {
            .dplyr_crop(atab,
                sdimx = sdimx,
                sdimy = sdimy,
                extent = bound,
                inclusive = inclusive
            )
        }
    }  
      
    storeRead(x, 
        output = output, 
        callback = crop_callback, 
        ...
    )
})

#' @rdname getBoundedData
#' @export
setMethod("getBoundedData", signature("parquetGeomStore", "SpatExtent"),
    function(x, bound, sdimx = "x_index", sdimy = "y_index", inclusive = TRUE, ...) {
    callNextMethod(x, bound, # to queryableStore method
        sdimx = sdimx,
        sdimy = sdimy,
        inclusive = inclusive,
        ...
    )
})

setMethod("getBoundedData", signature("fileStore", "ANY"), function(x, bound) {
    stop("[getBoundedData] fileStore input:",
    "\n * Consider first coercing to `queryableStore` if `storeRead()` outputs {dplyr} accessible interface.",
    "\n * Otherwise, first `storeWrite()` as `parquetStore`")
})

#' @name getTile
#' @title Get Tile Data
#' @description
#' Get the region of data defined by a {tilework} tile* object. This function
#' determines the tile(s) requested based on the tilePlan then calls
#' [getBoundedData()] for the actual retrieval step.
#' 
#' # Additional params:
#' 
#' * `advance` **`tileIterator` only** `logical` (default = TRUE). Whether to
#'   advance the iterator.
#' @param x `fileStore` object (currently `parquetStore` inheriting
#' only)
#' @param tiles {tilework} `tile*` object
#' @param i **ANY `tile*` except tileIterator** tile vector index or row index 
#' if `j` is also provided
#' @param j **ANY `tile*` except tileIterator** tile col index
#' @param contiguous `logical` (default = FALSE). Whether to retrieve tiles
#' using inclusivity rules (see [getBoundedData]) that prevent double counts
#' between neighboring tiles for values that fall on the boundary line.
#' 
#' When used, `get_params$inclusive` settings will be ignored.
#' 
#' This is mainly relevant only for padding = 0 situations with `tilePlan`
#' inheriting `tile` inputs and operations that will process all tiles in the
#' `tilePlan`. Inclusivity rules are calculated based on the entire set of tiles
#' defined in the `tilePlan`.
#' @param pad `numeric` (optional) additional padding to apply before tile retrieval.
#'   Useful for temporarily increasing padding without affecting `tile*` object.
#' @param ... additional named params passed through to [getBoundedData()]
#'   (e.g. `sdimx`, `sdimy`, `inclusive`, `envelope`, `output`).
#' @seealso [tilework::tilePlan]
#' @family tile methods
#' @returns a lazy arrow/dplyr query
NULL

#' @rdname getTile
#' @export
setMethod("getTile", signature("queryableStore", "spatialTilePlan"),
    function(x, tiles, i = NULL, j, contiguous = FALSE, pad = NULL, ...) {
    vmsg(.is_debug = TRUE, "[getTile] contiguous: ", contiguous)
    if (!is.null(pad)) {
        checkmate::assert_numeric(pad)
        tiles <- tiles + pad
    }

    if (!isTRUE(contiguous)) {
        if (missing(j)) {
            return(callNextMethod(x, tiles, i = i, ...))
        } else {
            return(callNextMethod(x, tiles, i = i, j = j, ...))
        }
    }

    # contiguous workflow
    if (missing(j)) { # ensure ij indexing
        ij <- .tile_idx_to_ij(tiles, i)
        i <- ij[[1L]]
        j <- ij[[2L]]
    } else {
        ij_tab <- expand.grid(j, i)
        i <- ij_tab[[2L]]
        j <- ij_tab[[1L]]
    }

    bounds_list <- tiles[i, j, expand_grid = FALSE]
    extra <- list(...)

    mapply(function(b, ri, ci) {
        inclusivity <- c(
            ri == 1, # bottom
            TRUE, # left
            TRUE, # top
            ci == ncol(tiles) # right
        )
        extra$inclusive <- inclusivity
        do.call(getBoundedData, c(list(x, b), extra))
    }, bounds_list, i, j, SIMPLIFY = FALSE)
})

#' @rdname getTile
#' @export
setMethod("getTile", signature("queryableStore", "freeTilePlan"), 
    function(x, tiles, i = NULL, contiguous = FALSE, pad = NULL, ...) {
    vmsg(.is_debug = TRUE, "[getTile] contiguous: ", contiguous)
    if (!is.null(pad)) {
        checkmate::assert_numeric(pad)
        tiles <- tiles + pad
    }

    if (!isTRUE(contiguous)) {
        return(callNextMethod(x, tiles, i = i, ...))
    }
      
    extra <- list(...)
      
    blims <- c(
        max(tiles$bounds[, 2L]), # xmax (right)
        min(tiles$bounds[, 3L]) # ymin (bottom)
    )
      
    lapply(i, function(ii) {
        inclusivity <- c(
            tiles$bounds[ii, 3L] == blims[2L], # bottom: tile ymin == global ymin
            TRUE, # left
            TRUE, # top
            tiles$bounds[ii, 2L] == blims[1L] # right: tile xmax == global xmax
        )
        extra$inclusive <- inclusivity
        b <- tiles[ii][[1L]]
        do.call(getBoundedData, c(list(x, b), extra))
    })
})

# internals ####

# Convert a vector i index to ij matrix indexing. Note that when using
# with {tilework} `[` subsetting, `expand_grid = FALSE` should be set.
# `x` accepts a `tilePlan`
# `i` is the index of the tile 
.tile_idx_to_ij <- function (x, i) {
    i_idx <- floor(i/ncol(x)) + 1L
    no_resid <- i%%ncol(x) == 0L
    i_idx[no_resid] <- i_idx[no_resid] - 1L
    j_idx <- i%%ncol(x)
    j_idx[j_idx == 0L] <- ncol(x)
    list(i_idx, j_idx)
}

# convert ij matrix indexing to vector i index.
# `i,j` i and (maybe) j indices to use
# `dims`` c(nrows, ncols) of tilePlan
# `expand_grid` - whether to use expand.grid()
.tile_ij_to_idx <- function(i = NULL, j = NULL, dims, expand_grid = TRUE) {
    if (is.null(i)) {
        if (is.null(j)) {
            stop("no i or j indices given", call. = FALSE)
        }
        if (isFALSE(expand_grid)) {
            stop("expand_grid may not be false when i is missing", call. = FALSE)
        }
        i <- seq_len(dims[[1L]])
    }
    if (is.null(j)) return(i) # already vector index, directly use
    if (expand_grid) {
        ij_tab <- expand.grid(j, i)
        i <- ij_tab[[2L]]
        j <- ij_tab[[1L]]
    }
    ((i - 1L) * dims[[2L]]) + j
}

# portions of .tile_crop and .tile_annotate are parallel implementations
# of `parquetStore` `getTile()` and `getBoundedData()`. This is because
# they are used in `parquetGeomTileStore` internals that may not
# always start with `parquetStore` source data (which currently requires
# a row index)

# envelope = TRUE means cropping by envelope centroids
.tile_crop <- function(tiles, data, i,
        sdimx = "x",
        sdimy = "y",
        group_col = "id",
        envelope = FALSE) {
    if (isTRUE(envelope)) {
        use_data <- .dplyr_xy_envelopes(data,
            sdimx = sdimx,
            sdimy = sdimy,
            group_col = group_col
        )
        sdimx = "ecentroid_x"
        sdimy = "ecentroid_y"
    } else {
        use_data <- data
    }

    ij <- .tile_idx_to_ij(tiles, i)
    right_edge <- ij[[2]] == ncol(tiles)
    bottom_edge <- ij[[1]] == nrow(tiles)
    crop_inclusion <- c(bottom_edge, TRUE, TRUE, right_edge)
    res <- use_data %>%
        .dplyr_crop(
            extent = tiles[i][[1]],
            sdimx = sdimx,
            sdimy = sdimy,
            inclusive = crop_inclusion
        )
    if (!isTRUE(envelope)) return(res)

    # select actual rows to pull based on envelopes
    dplyr::semi_join(data, res, by = stats::setNames("id", group_col))
}

# setup a tilePlan and annotate with tile metadata
# `$n_records` - count of items per tile
.tile_annotate <- function(tiles, data,
        n_tiles = 100,
        sdimx,
        sdimy,
        group_col = NULL,
        envelope = TRUE) {
    checkmate::assert_class(tiles, "tilePlan")
    checkmate::assert_character(sdimx, len = 1L)
    checkmate::assert_character(sdimy, len = 1L)

    ext(tiles) <- .dplyr_ext(data, sdimx = sdimx, sdimy = sdimy)
    length(tiles) <- n_tiles
    n_tiles <- length(tiles) # actual length may be different
    tiles$n_records <- NA_integer_

    # get spatial envelope centroids
    if (isTRUE(envelope)) {
        checkmate::assert_character(group_col, len = 1L)
        data <- .dplyr_xy_envelopes(data,
            sdimx = sdimx,
            sdimy = sdimy,
            group_col = group_col
        )
        sdimx = "ecentroid_x"
        sdimy = "ecentroid_y"
    }
  
    count <- vapply(seq_len(n_tiles), function(i) {
        tile_data <- .tile_crop(tiles,
            data = data,
            i = i,
            sdimx = sdimx,
            sdimy = sdimy,
            envelope = FALSE
        )
        .dplyr_nrow(tile_data)
    }, FUN.VALUE = double(1L))
    tiles$n_records <- count
    tiles
}
