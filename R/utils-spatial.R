# Morton (Z-order) encode 2D coordinates to a sortable integer key.
# ext: numeric(4) c(xmin, xmax, ymin, ymax) — global bounds for normalisation.
.morton_encode <- function(x, y, ext) {
    xi <- as.integer((x - ext[[1L]]) / (ext[[2L]] - ext[[1L]]) * (2^16 - 1))
    yi <- as.integer((y - ext[[3L]]) / (ext[[4L]] - ext[[3L]]) * (2^16 - 1))
    spread <- function(v) {
        v <- bitwAnd(v, 0x0000FFFF)
        v <- bitwAnd(bitwOr(v, bitwShiftL(v, 8L)), 0x00FF00FF)
        v <- bitwAnd(bitwOr(v, bitwShiftL(v, 4L)), 0x0F0F0F0F)
        v <- bitwAnd(bitwOr(v, bitwShiftL(v, 2L)), 0x33333333)
        v <- bitwAnd(bitwOr(v, bitwShiftL(v, 1L)), 0x55555555)
        v
    }
    bitwOr(spread(xi), bitwShiftL(spread(yi), 1L))
}

# Strip names from a SpatExtent, returning a plain numeric vector
.ext_to_num_vec <- function(x) {
    out <- x[]
    names(out) <- NULL
    out
}

# Replace +/-Inf extent bounds with a large finite value so terra can build polygons.
.clamp_ext_infinite <- function(e, big = 1e15) {
    v <- .ext_to_num_vec(e)
    v[1L] <- max(v[1L], -big)  # xmin
    v[2L] <- min(v[2L],  big)  # xmax
    v[3L] <- max(v[3L], -big)  # ymin
    v[4L] <- min(v[4L],  big)  # ymax
    ext(v[seq(2L)], v[seq(3L, 4L)])
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

    ext(c(ranges$x_min, ranges$x_max), c(ranges$y_min, ranges$y_max))
}

# Extent of affine-transformed coordinates computed directly in Arrow (post-multiply).
# Must always query -- disk_extent is in intrinsic space.
.dplyr_ext_affine <- function(data, aff,
        sdimx = "x_index",
        sdimy = "y_index") {
    m <- aff@affine  # 3x3 homogeneous
    a11 <- m[1L, 1L]; a21 <- m[2L, 1L]; tx <- m[1L, 3L]
    a12 <- m[1L, 2L]; a22 <- m[2L, 2L]; ty <- m[2L, 3L]

    x_sym <- as.name(sdimx)
    y_sym <- as.name(sdimy)

    # inline coefficients so Arrow sees numeric literals, not R env vars
    xout <- rlang::expr(!!x_sym * !!a11 + !!y_sym * !!a21 + !!tx)
    yout <- rlang::expr(!!x_sym * !!a12 + !!y_sym * !!a22 + !!ty)

    ranges <- data |>
        dplyr::summarize(
            x_min = min(!!xout, na.rm = TRUE),
            x_max = max(!!xout, na.rm = TRUE),
            y_min = min(!!yout, na.rm = TRUE),
            y_max = max(!!yout, na.rm = TRUE)
        ) |>
        dplyr::collect()

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

# Convert a data.frame of points to list(geom = SpatVector, meta = data.table).
# Metadata cols are separated so terra never ingests them.
.df_to_terra_pts <- function(x,
    id_col = NULL,
    sdimx = NULL,
    sdimy = NULL,
    feat_ID_colname = NULL, # sink: prevent duplicate arg if passed via ... 
    x_colname = NULL, # sink: prevent duplicate arg if passed via ... 
    y_colname = NULL, # sink: prevent duplicate arg if passed via ... 
    ...) {
    checkmate::assert_data_frame(x)
    data.table::setDT(x)
    
    geom_cols <- c(id_col, sdimx, sdimy)
    meta_cols <- c(id_col, setdiff(colnames(x), geom_cols))
    x_geom <- x[, geom_cols, with = FALSE]
    x_meta <- x[, meta_cols, with = FALSE]
  
    x_sv <- GiottoClass::createGiottoPoints(x_geom,
        feat_ID_colname = id_col, # not actually geom id
        x_colname = sdimx,
        y_colname = sdimy,
        split_keyword = NULL, # not done at this level
        ...
    )[] # drop to spatvector
    list(
        geom = x_sv,
        meta = x_meta
    )
}

.df_to_terra_poly <- function(x,
    id_col = NULL,
    part_col = NULL,
    sdimx = NULL,
    sdimy = NULL,
    hole_col = NULL,
    ...) {
    checkmate::assert_data_frame(x)
    data.table::setDT(x)

    geom_cols <- c(id_col, sdimx, sdimy, part_col, hole_col)
    meta_cols <- c(id_col, setdiff(colnames(x), geom_cols))
    x_geom <- x[, geom_cols, with = FALSE]
    x_meta <- x[, meta_cols, with = FALSE]

    x_sv <- GiottoClass::createGiottoPolygon(x_geom,
        part_col = part_col,
        calc_centroids = FALSE,
        make_valid = TRUE,
        ...
    )[] # drop to spatvector

    # createGiottoPolygon may drop polygons (invalid geom under
    # make_valid = TRUE, degenerate parts, etc.) — align meta to the
    # surviving poly_IDs so `cbind(data, meta)` in
    # `.terra_to_parquet_format` doesn't mismatch row counts. SV's id
    # column is always "poly_ID" after createGiottoPolygon; x_meta still
    # uses the user-provided `id_col` name.
    sv_ids <- terra::values(x_sv)[["poly_ID"]]
    meta_dedup <- unique(x_meta)
    meta_aligned <- meta_dedup[
        match(sv_ids, meta_dedup[[id_col]]), ,
        drop = FALSE
    ]

    list(
        geom = x_sv,
        meta = meta_aligned
    )
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

# affine transform helpers ####

# TRUE when the 2x2 linear submatrix has non-zero off-diagonal elements,
# meaning the transform includes rotation or shear.
.affine_has_rotation <- function(aff) {
    m <- .aff_linear_2d(aff)
    !isTRUE(all.equal(m[1L, 2L], 0, tolerance = .Machine$double.eps^0.5)) ||
        !isTRUE(all.equal(m[2L, 1L], 0, tolerance = .Machine$double.eps^0.5))
}

# Axis-aligned bounding box of an nx2 corner matrix.
# Returns numeric(4): c(xmin, xmax, ymin, ymax).
.affine_aabb <- function(corners) {
    c(min(corners[, 1L]), max(corners[, 1L]),
      min(corners[, 2L]), max(corners[, 2L]))
}

# Build a compound R call expressing the interior of a convex parallelogram
# as 4 half-plane inequalities on x_index and y_index.
# corners: 4x2 numeric matrix of (x, y) polygon corners in order.
# Returns: a1*x_index + b1*y_index >= c1 & ... (4 conditions)
.affine_halfplane_expr <- function(corners) {
    # Ensure CCW winding via signed area (shoelace formula)
    x <- corners[, 1L]
    y <- corners[, 2L]
    n <- nrow(corners)
    i_next <- c(seq(2L, n), 1L)
    signed_area <- sum(x * y[i_next] - x[i_next] * y) / 2
    if (signed_area < 0) corners <- corners[rev(seq_len(n)), ]

    conditions <- vector("list", n)
    for (i in seq_len(n)) {
        p1 <- corners[i, ]
        p2 <- corners[(i %% n) + 1L, ]
        d <- p2 - p1
        # Inward normal for CCW polygon: rotate edge 90deg CCW = (-dy, dx)
        nx <- -d[2L]
        ny <-  d[1L]
        c_val <- nx * p1[1L] + ny * p1[2L]
        conditions[[i]] <- call(">=",
            call("+",
                call("*", nx, as.name("x_index")),
                call("*", ny, as.name("y_index"))
            ),
            c_val
        )
    }
    Reduce(function(a, b) call("&", a, b), conditions)
}

# Internal accessor for 2x2 linear submatrix -- avoids importing unexported
# GiottoClass internals.
.aff_linear_2d <- function(x) {
    if (inherits(x, "affine2d")) x <- x@affine
    x[seq(2L), seq(2L)]
}

# Fast metadata-based extent estimate for parquetGeomBase -- no row scan.
# Precedence: @crop > disk_extent > live scan; @window intersected on top.
# Pending affine projects AABB corners: exact for axis-aligned, conservative
# overestimate for rotation/shear. Row-level ops are not reflected.
.pgeom_ext_estimate <- function(x, aff = .pgeom_pending_transform(x)) {
    # Tightest available intrinsic upper bound
    if (length(x@crop) > 0L) {
        # @crop already intersects disk_extent at record time -- always tighter
        e <- ext(x@crop)
    } else {
        de <- .pstore_disk_extent(x)
        if (!is.null(de)) {
            e <- ext(de)
        } else {
            # No metadata to work from -- fall back to live scan
            return(.pgeom_ext_intrinsic(x))
        }
    }

    if (length(x@window) > 0L) {
        e <- terra::intersect(e, ext(x@window))
        if (is.null(e)) return(NULL)
    }

    if (is.null(aff)) return(e)

    # Project corners through affine and return the AABB.
    ext(affine(terra::as.polygons(e), aff))
}

# Apply @post_ops to a materialized result in order.
# Supported types: "transform" (affine2d), "geom_filter" (future exact polygon mask).
# Unknown op types are silently skipped.
.apply_post_ops <- function(result, post_ops, output) {
    for (op in post_ops) {
        if (op$type == "transform") {
            aff <- op$affine2d
            result <- switch(output,
                "terra" = affine(result, aff),
                "sf" = {
                    sv <- terra::vect(result)
                    sv <- affine(sv, aff)
                    sf::st_as_sf(sv)
                },
                "tibble" = {
                    if (all(c("x_index", "y_index") %in% names(result))) {
                        # post-multiply convention (pre_multiply = FALSE),
                        # matching GiottoClass .affine_matrix() default
                        m <- aff@affine          # 3x3 homogeneous
                        A <- m[seq(2L), seq(2L)] # 2x2 linear
                        t_vec <- m[seq(2L), 3L]  # translation
                        xy <- as.matrix(result[, c("x_index", "y_index")])
                        xy_new <- xy %*% A
                        xy_new[, 1L] <- xy_new[, 1L] + t_vec[1L]
                        xy_new[, 2L] <- xy_new[, 2L] + t_vec[2L]
                        result$x_index <- xy_new[, 1L]
                        result$y_index <- xy_new[, 2L]
                    }
                    result
                }
            )
        }
    }
    result
}
