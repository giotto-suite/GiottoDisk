# Helper function to convert a SpatExtent object to a simple numeric vector
# and strip the names
.ext_to_num_vec <- function(x) {
    out <- x[]
    names(out) <- NULL
    out
}

.dplyr_ext <- function(data, sdimx = "x_index", sdimy = "y_index") {
    ranges <- data %>%
        dplyr::select(dplyr::all_of(c(sdimx, sdimy))) %>%
        dplyr::summarize(
            x_min = min(!!as.name(sdimx), na.rm = TRUE),
            x_max = max(!!as.name(sdimx), na.rm = TRUE),
            y_min = min(!!as.name(sdimy), na.rm = TRUE),
            y_max = max(!!as.name(sdimy), na.rm = TRUE)
        ) %>%
        dplyr::collect()

    # Create the extent object
    ext(c(ranges$x_min, ranges$x_max), c(ranges$y_min, ranges$y_max))
}

# find xy bounds + centroids of bounds
.dplyr_xy_envelopes <- function(data,
        sdimx = "x_index",
        sdimy = "y_index",
        group_col = "poly_ID"
    ) {
    centroids <- data %>%
        dplyr::group_by(!!as.name(group_col)) %>%
        dplyr::summarize(
            xmin = min(!!as.name(sdimx), na.rm = TRUE),
            xmax = max(!!as.name(sdimx), na.rm = TRUE),
            ymin = min(!!as.name(sdimy), na.rm = TRUE),
            ymax = max(!!as.name(sdimy), na.rm = TRUE),
            .groups = "drop"
        ) %>%
        # Calculate envelope centroids
        dplyr::mutate(
            ecentroid_x = (xmin + xmax) / 2,
            ecentroid_y = (ymin + ymax) / 2
        ) %>%
        dplyr::rename(id = !!as.name(group_col))

    return(centroids)
}

# bottom, left, top, right
.dplyr_crop <- function(data, extent,
        sdimx = "x_index",
        sdimy = "y_index",
        inclusive = c(FALSE, TRUE, TRUE, FALSE)
    ) {
    if (length(inclusive) == 1L) inclusive <- rep(inclusive, 4)
    e <- .ext_to_num_vec(ext(extent))
    if (inclusive[[1]]) {
        data <- dplyr::filter(data, !!as.name(sdimy) >= e[3])
    } else {
        data <- dplyr::filter(data, !!as.name(sdimy) > e[3])
    }

    if (inclusive[[2]]) {
        data <- dplyr::filter(data, !!as.name(sdimx) >= e[1])
    } else {
        data <- dplyr::filter(data, !!as.name(sdimx) > e[1])
    }

    if (inclusive[[3]]) {
        data <- dplyr::filter(data, !!as.name(sdimy) <= e[4])
    } else {
        data <- dplyr::filter(data, !!as.name(sdimy) < e[4])
    }

    if (inclusive[[4]]) {
        data <- dplyr::filter(data, !!as.name(sdimx) <= e[2])
    } else {
        data <- dplyr::filter(data, !!as.name(sdimx) < e[2])
    }
    data
}
