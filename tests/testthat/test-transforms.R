make_grid_pts <- function(path = tempfile()) {
    df <- expand.grid(x = 1:3, y = 1:3)
    df$cell_ID <- paste0("c", seq_len(nrow(df)))
    sv <- terra::vect(df, geom = c("x", "y"), crs = "")
    s <- parquetGeomStore(path)
    storeWrite(s, sv)
}

test_that("affine: terra output has transformed coordinates", {
    s <- make_grid_pts()
    aff <- GiottoClass::affine(diag(c(2, 2)))
    s2 <- affine(s, aff)
    sv <- storeRead(s2, output = "terra")
    e <- terra::ext(sv)
    expect_equal(as.numeric(e$xmax), 6, tolerance = 1e-6)
    expect_equal(as.numeric(e$ymax), 6, tolerance = 1e-6)
})

test_that("affine -> crop (axis-aligned): output within crop bounds", {
    s <- make_grid_pts()
    aff <- GiottoClass::affine(diag(c(2, 2)))
    s2 <- affine(s, aff)
    # grid after scale-2 spans x/y in [2, 6]; crop to [0, 5] excludes x=3 points
    s3 <- crop(s2, terra::ext(0, 5, 0, 5))
    sv <- storeRead(s3, output = "terra")
    expect_gt(nrow(sv), 0)
    expect_lt(nrow(sv), 9)
    coords <- terra::crds(sv)
    expect_true(all(coords[, 1] <= 5 + 1e-6))
    expect_true(all(coords[, 2] <= 5 + 1e-6))
})

test_that("affine -> crop -> affine: coordinates correct after chain", {
    s <- make_grid_pts()  # 3x3 grid, x/y in {1,2,3}
    aff1 <- GiottoClass::affine(diag(c(1, 1)))  # identity
    s2 <- affine(s, aff1)
    # crop to x <= 2.5 in output space (= input space for identity) -- keeps x in {1,2}
    s3 <- crop(s2, terra::ext(0, 2.5, 0, 4))
    aff2 <- GiottoClass::affine(diag(c(2, 2)))  # scale-2
    s4 <- affine(s3, aff2)
    sv <- storeRead(s4, output = "terra")
    # identity then scale-2 = scale-2; x in {1,2} -> {2,4}, y in {1,2,3} -> {2,4,6}
    expect_equal(nrow(sv), 6L)
    coords <- terra::crds(sv)
    expect_true(all(coords[, 1] %in% c(2, 4)))
    expect_true(all(coords[, 2] %in% c(2, 4, 6)))
})

test_that("affine -> crop -> affine: composes to single transform op", {
    s <- make_grid_pts()
    aff2 <- GiottoClass::affine(diag(c(2, 2)))
    s2 <- affine(s, aff2)
    # crop covers entire transformed extent -- no rows excluded
    s3 <- crop(s2, terra::ext(0, 7, 0, 7))
    aff_half <- GiottoClass::affine(diag(c(0.5, 0.5)))
    s4 <- affine(s3, aff_half)
    post_types <- vapply(s4@post_ops, function(op) op$type, character(1L))
    op_types <- vapply(s4@ops, function(op) op$type, character(1L))
    # composition must collapse to exactly one "transform" entry in @post_ops
    expect_equal(sum(post_types == "transform"), 1L)
    expect_false("filter" %in% op_types)
    sv <- storeRead(s4, output = "terra")
    expect_gt(nrow(sv), 0)
})

test_that("ext(): exact=FALSE is AABB containing exact=TRUE; half-plane filter reflected by exact", {
    s <- make_grid_pts()
    theta <- pi / 4  # 45 degrees
    rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
    aff <- GiottoClass::affine(rot)
    s2 <- affine(s, aff)
    e_exact <- ext(s2, exact = TRUE)
    e_fast  <- ext(s2, exact = FALSE)
    expect_s4_class(e_exact, "SpatExtent")
    expect_s4_class(e_fast,  "SpatExtent")
    # fast AABB must contain (or equal) the exact scanned extent
    expect_lte(e_fast$xmin, e_exact$xmin + 1e-6)
    expect_gte(e_fast$xmax, e_exact$xmax - 1e-6)
    expect_lte(e_fast$ymin, e_exact$ymin + 1e-6)
    expect_gte(e_fast$ymax, e_exact$ymax - 1e-6)

    # Without affine: crop to x <= 2.5 keeps x in {1,2} (actual xmax=2).
    # exact=TRUE scans actual data; exact=FALSE uses @crop AABB (xmax=2.5).
    s3 <- crop(s, terra::ext(0, 2.5, 0, 4))
    e_exact_filtered <- ext(s3, exact = TRUE)
    e_fast_filtered  <- ext(s3, exact = FALSE)
    expect_equal(as.numeric(e_exact_filtered$xmax), 2, tolerance = 1e-6)
    expect_lte(e_exact_filtered$xmax, e_fast_filtered$xmax + 1e-6)
})

test_that("rotation -> crop: half-plane filter injected into @ops", {
    s <- make_grid_pts()
    theta <- pi / 4
    rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
    aff <- GiottoClass::affine(rot)
    s2 <- affine(s, aff)
    # crop in rotated space; grid after 45-deg rotation spans x ~[1.4, 4.2], y ~[-1.4, 1.4]
    s3 <- crop(s2, terra::ext(1, 4, -1, 1))
    op_types <- vapply(s3@ops, function(op) op$type, character(1L))
    post_types <- vapply(s3@post_ops, function(op) op$type, character(1L))
    expect_true("filter" %in% op_types)
    expect_true("transform" %in% post_types)
    tbl <- storeRead(s3, output = "tibble")
    expect_gt(nrow(tbl), 0)
})

test_that("affine: tibble output has updated x_index/y_index", {
    s <- make_grid_pts()
    aff <- GiottoClass::affine(diag(c(2, 2)))
    s2 <- affine(s, aff)
    tbl <- storeRead(s2, output = "tibble")
    expect_true("x_index" %in% names(tbl))
    expect_true("y_index" %in% names(tbl))
    # original x_index/y_index values are 1, 2, 3; scale-2 doubles them
    expect_true(all(tbl$x_index %in% c(2, 4, 6)))
    expect_true(all(tbl$y_index %in% c(2, 4, 6)))
})

test_that("spatShift: translates coordinates", {
    s <- make_grid_pts()
    s2 <- spatShift(s, dx = 10, dy = 20)
    sv <- storeRead(s2, output = "terra")
    coords <- terra::crds(sv)
    # original x/y in {1,2,3}; shifted by +10, +20
    expect_true(all(coords[, 1] %in% c(11, 12, 13)))
    expect_true(all(coords[, 2] %in% c(21, 22, 23)))
})

test_that("rescale: scales coordinates around centroid", {
    s <- make_grid_pts()
    # scale-2 with default centroid (centroid of {1,2,3} x {1,2,3} = (2,2))
    # x' = 2 + 2*(x-2); y' = 2 + 2*(y-2)
    s2 <- rescale(s, fx = 2)
    sv <- storeRead(s2, output = "terra")
    coords <- terra::crds(sv)
    expect_equal(as.numeric(terra::ext(sv)$xmin), 0, tolerance = 1e-6)
    expect_equal(as.numeric(terra::ext(sv)$xmax), 4, tolerance = 1e-6)
    expect_equal(nrow(sv), 9L)
})

test_that("spin: rotates coordinates around centroid", {
    s <- make_grid_pts()
    # 180-degree rotation around centroid (2,2): output extent same as input
    s2 <- spin(s, angle = 180)
    sv <- storeRead(s2, output = "terra")
    e <- terra::ext(sv)
    expect_equal(nrow(sv), 9L)
    expect_equal(as.numeric(e$xmin), 1, tolerance = 1e-6)
    expect_equal(as.numeric(e$xmax), 3, tolerance = 1e-6)
    expect_equal(as.numeric(e$ymin), 1, tolerance = 1e-6)
    expect_equal(as.numeric(e$ymax), 3, tolerance = 1e-6)
})

test_that("flip: reflects coordinates vertically", {
    s <- make_grid_pts()
    # vertical flip around y=0: y' = -y; x unchanged
    s2 <- flip(s, direction = "vertical", y0 = 0)
    sv <- storeRead(s2, output = "terra")
    coords <- terra::crds(sv)
    expect_true(all(coords[, 1] %in% c(1, 2, 3)))
    expect_true(all(coords[, 2] %in% c(-1, -2, -3)))
})

test_that("t: swaps x and y coordinates", {
    # non-square grid makes the swap observable: x in [1,2], y in [1,4]
    sv <- terra::vect(expand.grid(x = 1:2, y = 1:4), geom = c("x", "y"), crs = "")
    sv$cell_ID <- paste0("c", seq_len(nrow(sv)))
    s <- parquetGeomStore(tempfile())
    storeWrite(s, sv)
    s2 <- t(s)
    sv2 <- storeRead(s2, output = "terra")
    e2 <- terra::ext(sv2)
    expect_equal(nrow(sv2), 8L)
    # after transpose: x in [1,4], y in [1,2]
    expect_equal(as.numeric(e2$xmax), 4, tolerance = 1e-6)
    expect_equal(as.numeric(e2$ymax), 2, tolerance = 1e-6)
})

test_that("shear: shears coordinates", {
    s <- make_grid_pts()
    # x-shear fx=1 with default centre (2,2): x' = x + (y-2)*1
    # x=1,y=1 -> x'=0; x=1,y=2 -> x'=1; x=1,y=3 -> x'=2
    s2 <- shear(s, fx = 1)
    sv <- storeRead(s2, output = "terra")
    expect_equal(nrow(sv), 9L)
    e <- terra::ext(sv)
    # x range must be wider than [1,3] after shear
    expect_lt(as.numeric(e$xmin), 1 + 1e-6)
    expect_gt(as.numeric(e$xmax), 3 - 1e-6)
})

test_that("crop -> spin: centroid is live-scanned from cropped extent", {
    s <- make_grid_pts()  # 3x3 grid, centroid (2,2)
    # crop to x <= 2.5 keeps x in {1,2}; cropped centroid is (1.5,2) not (2,2)
    s2 <- crop(s, terra::ext(0, 2.5, 0, 4))
    s3 <- spin(s2, angle = 180)
    sv <- storeRead(s3, output = "terra")
    # 180-deg spin around (1.5,2): x' = 3-x, y' = 4-y
    # x in {1,2} -> {2,1}; result xmax=2, not 3 (which would indicate wrong centroid)
    expect_equal(nrow(sv), 6L)
    expect_equal(as.numeric(terra::ext(sv)$xmax), 2, tolerance = 1e-6)
})

test_that("crop -> rescale: centroid is live-scanned from cropped extent", {
    s <- make_grid_pts()  # 3x3 grid, centroid (2,2)
    # crop to x <= 2.5 keeps x in {1,2}; cropped centroid is (1.5,2)
    s2 <- crop(s, terra::ext(0, 2.5, 0, 4))
    s3 <- rescale(s2, fx = 2)
    sv <- storeRead(s3, output = "terra")
    # rescale-2 around (1.5,2): x' = 1.5 + 2*(x-1.5) = 2x-1.5
    # x=1 -> 0.5, x=2 -> 2.5; xmin=0.5, xmax=2.5
    expect_equal(nrow(sv), 6L)
    expect_equal(as.numeric(terra::ext(sv)$xmin), 0.5, tolerance = 1e-6)
    expect_equal(as.numeric(terra::ext(sv)$xmax), 2.5, tolerance = 1e-6)
})

test_that("window -> spin: centroid is live-scanned from windowed extent", {
    s <- make_grid_pts()  # 3x3 grid, centroid (2,2)
    # window to x <= 2.5 keeps x in {1,2}; windowed centroid is (1.5,2)
    window(s) <- terra::ext(0, 2.5, 0, 4)
    s2 <- spin(s, angle = 180)
    sv <- storeRead(s2, output = "terra")
    # same geometry as crop -> spin: xmax=2
    expect_equal(nrow(sv), 6L)
    expect_equal(as.numeric(terra::ext(sv)$xmax), 2, tolerance = 1e-6)
})

test_that("verb chain: spatShift -> spin -> rescale composes to single transform op", {
    s <- make_grid_pts()
    s2 <- spatShift(s, dx = 1, dy = 1)
    s3 <- spin(s2, angle = 90)
    s4 <- rescale(s3, fx = 0.5)
    post_types <- vapply(s4@post_ops, function(op) op$type, character(1L))
    expect_equal(sum(post_types == "transform"), 1L)
    sv <- storeRead(s4, output = "terra")
    expect_equal(nrow(sv), 9L)
})
