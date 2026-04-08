# Helper function to convert a SpatExtent object to a simple numeric vector
# and strip the names
.ext_to_num_vec <- function(x) {
    out <- x[]
    names(out) <- NULL
    out
}

.dplyr_ext <- function(data, sdimx = "x_index", sdimy = "y_index") {
    ranges <- data |>
        dplyr::summarize(
            x_min = min(!!as.name(sdimx), na.rm = TRUE),
            x_max = max(!!as.name(sdimx), na.rm = TRUE),
            y_min = min(!!as.name(sdimy), na.rm = TRUE),
            y_max = max(!!as.name(sdimy), na.rm = TRUE)
        ) |>
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
    centroids <- data |>
        dplyr::group_by(!!as.name(group_col)) |>
        dplyr::summarize(
            xmin = min(!!as.name(sdimx), na.rm = TRUE),
            xmax = max(!!as.name(sdimx), na.rm = TRUE),
            ymin = min(!!as.name(sdimy), na.rm = TRUE),
            ymax = max(!!as.name(sdimy), na.rm = TRUE),
            .groups = "drop"
        ) |>
        # Calculate envelope centroids
        dplyr::mutate(
            ecentroid_x = (xmin + xmax) / 2,
            ecentroid_y = (ymin + ymax) / 2
        ) |>
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

# create the geoparquet json metadata
.geoparquet_metadata <- function(geom_col = "geom",                  
    geomtype = NULL,                                                 
    crs = NULL,                                                      
    extent = NULL) {
                                                                      
    # terra geomtype -> GeoParquet OGC geometry type string
    .geomtype_map <- c(                                              
        "points"       = "Point",
        "lines"        = "LineString",
        "polygons"     = "Polygon",
        "multipoints"  = "MultiPoint",
        "multilines"   = "MultiLineString",
        "multipolygons" = "MultiPolygon"
    )

    col_meta <- list(encoding = "WKB")

    # empty list [] = unknown/mixed, which is valid per spec
    if (!is.null(geomtype) && nzchar(geomtype)) {
        mapped <- .geomtype_map[tolower(geomtype)]
        col_meta$geometry_types <- if (!is.na(mapped)) {
            list(unname(mapped)) 
        } else {
            list()
        }
        
    } else {
        col_meta$geometry_types <- list()
    }

    # omit crs field entirely if not set (valid per spec)
    if (!is.null(crs) && nzchar(crs)) {
        col_meta$crs <- crs
    }

    # terra/store@extent order: xmin, xmax, ymin, ymax
    # GeoParquet bbox order: [xmin, ymin, xmax, ymax]
    if (!is.null(extent) && length(extent) == 4L) {
        bbox <- c(extent[[1L]], extent[[3L]], extent[[2L]], extent[[4L]])
        col_meta$bbox <- bbox
    }

    geo_meta <- list(
        version = "1.0.0",
        primary_column = geom_col,
        columns = stats::setNames(list(col_meta), geom_col)
    )

    as.character(jsonlite::toJSON(geo_meta, auto_unbox = TRUE))
}
