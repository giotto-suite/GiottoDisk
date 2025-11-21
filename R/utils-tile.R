#' @include utils-spatial.R
NULL


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

.tile_idx_to_ij <- function (x, i) {
    i_idx <- floor(i/ncol(x)) + 1L
    no_resid <- i%%ncol(x) == 0L
    i_idx[no_resid] <- i_idx[no_resid] - 1L
    j_idx <- i%%ncol(x)
    j_idx[j_idx == 0L] <- ncol(x)
    list(i_idx, j_idx)
}

.annotate_tileiterator <- function(tiles, data,
        n_tiles = 100,
        sdimx = "x",
        sdimy = "y",
        poly_id = "id",
        envelope = TRUE) {
    checkmate::assert_class(tiles, "tileIterator")
    ext(tiles) <- .dplyr_ext(data, sdimx = sdimx, sdimy = sdimy)
    length(tiles) <- n_tiles
    n_tiles <- length(tiles) # actual length may be different
    tiles$n_records <- NA_integer_
    tiles$row_offset <- NA_integer_
    tiles$row_offset[1] <- 0L

    # get spatial envelope centroids
    if (isTRUE(envelope)) {
        data <- .dplyr_xy_envelopes(data,
            sdimx = sdimx,
            sdimy = sdimy,
            group_col = poly_id
        )
        sdimx = "ecentroid_x"
        sdimy = "ecentroid_y"
    }
    for (i in seq_len(n_tiles)) {
        # determine offsets
        tile_data <- .tile_crop(tiles,
            data = data,
            i = i,
            sdimx = sdimx,
            sdimy = sdimy,
            envelope = FALSE
        )
        tile_rows <- .dplyr_nrow(tile_data)
        tiles$n_records[i] <- tile_rows
        if (i < n_tiles) {
            tiles$row_offset[i + 1L] <- tile_rows + tiles$row_offset[i]
        }
    }
    tiles
}

