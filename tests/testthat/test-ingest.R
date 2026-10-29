# The optional Rust ingest backend (GiottoIngest): every tiled transcript
# write it takes over must give the R writer's store (plan, tile dirs, rows
# in every column but the not-reproducible row_index, geoparquet metadata,
# store slots), and the R path must run whenever the gate says no.

options("tilework.warn_sequential" = FALSE)

.skip_no_ingest <- function() {
    skip_if_not_installed("GiottoIngest")
    skip_if_not_installed("duckdb")
    skip_if_not(.gd_ingest("points_parquet"))
}

# Run `write(gsource)` with the Rust path off, then on; check that the second
# run did go through GiottoIngest's writers and the first did not.
.both_paths <- function(write, rust = TRUE) {
    calls <- 0L
    count <- function(f) {
        force(f)  # the original, before it is mocked
        function(...) { calls <<- calls + 1L; f(...) }
    }
    testthat::local_mocked_bindings(
        gi_write_point_tiles = count(GiottoIngest::gi_write_point_tiles),
        gi_atera_tx_cells = count(GiottoIngest::gi_atera_tx_cells),
        .package = "GiottoIngest")
    r <- withr::with_options(list(giottodisk.use_ingest = FALSE),
        write(gDirSource(path = tempfile("r_"))))
    expect_identical(calls, 0L)
    rs <- withr::with_options(list(giottodisk.use_ingest = TRUE),
        write(gDirSource(path = tempfile("rs_"))))
    if (rust) expect_gt(calls, 0L) else expect_identical(calls, 0L)
    list(r = r, rs = rs)
}

.tile_dirs <- function(st) {
    d <- list.dirs(.idpath(st@path, st@uid), recursive = FALSE)
    sort(d[grepl("^tile_index=", basename(d))])
}

# rows of a tile, ordered by every column (no unique key in general)
.tile_rows <- function(d) {
    x <- as.data.frame(arrow::open_dataset(d))
    x <- x[, setdiff(names(x), c("row_index", "geom")), drop = FALSE]
    x <- x[do.call(order, unname(as.list(x))), , drop = FALSE]
    rownames(x) <- NULL
    x
}

# the file's own key-value metadata, where GeoParquet readers (DuckDB spatial,
# SedonaDB) look; arrow would also find a key kept only in ARROW:schema
.geo <- function(d) {
    f <- list.files(d, pattern = "parquet$", full.names = TRUE)[1L]
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    kv <- DBI::dbGetQuery(con, sprintf(
        "SELECT decode(key) AS key, decode(value) AS value FROM parquet_kv_metadata(%s)",
        DBI::dbQuoteString(con, f)))
    kv$value[kv$key == "geo"]
}

.expect_same_store <- function(b) {
    r <- b$r; rs <- b$rs
    expect_s4_class(rs, "parquetGeomTileStore")
    expect_identical(rs@tiles@bounds, r@tiles@bounds)
    expect_identical(rs@tiles@metadata, r@tiles@metadata)
    expect_identical(rs@geomtype, r@geomtype)
    expect_identical(rs@params, r@params)
    rd <- .tile_dirs(r); sd <- .tile_dirs(rs)
    expect_identical(basename(sd), basename(rd))
    for (k in seq_along(rd)) {
        expect_identical(.tile_rows(sd[k]), .tile_rows(rd[k]), label = basename(rd[k]))
        expect_identical(.geo(sd[k]), .geo(rd[k]))
    }
    # the Rust store reads back like the R one
    expect_identical(nrow(rs), nrow(r))
    expect_identical(as.vector(terra::ext(rs)), as.vector(terra::ext(r)))
}

.tx_parquet <- function(n = 3000L) {
    set.seed(11)
    df <- data.frame(
        transcript_id = seq_len(n) * 7,
        cell_id = sprintf("c%03d", seq_len(n) %% 97L),
        feature_name = sample(c("A", "B", "NegControlProbe_1"), n, TRUE),
        x_location = round(stats::runif(n, 0, 100), 3),
        y_location = round(stats::runif(n, 0, 50), 3),
        qv = round(stats::runif(n, 0, 40), 1),
        fov_name = "A1"
    )
    tbl <- arrow::arrow_table(df)
    tbl$transcript_id <- tbl$transcript_id$cast(arrow::uint64())
    tbl$x_location <- tbl$x_location$cast(arrow::float32())
    tbl$y_location <- tbl$y_location$cast(arrow::float32())
    f <- tempfile(fileext = ".parquet")
    arrow::write_parquet(tbl, f, chunk_size = 400L)
    f
}

test_that("Xenium transcripts: the Rust store is the R store", {
    .skip_no_ingest()
    f <- .tx_parquet()
    b <- .both_paths(function(src) .xenium_transcript_disk(f, src,
        qv_threshold = 20, dropcols = "fov_name", output = "store",
        threshold = 300))
    expect_gt(nrow(b$r@tiles@bounds), 3L)
    .expect_same_store(b)
})

test_that("Xenium zarr transcripts (converted to parquet shards) match", {
    .skip_no_ingest()
    skip_if_no_zarr_deps()
    fx <- make_zarr_fixture()
    withr::local_options(giottodisk.artifact_dump = withr::local_tempdir())
    b <- .both_paths(function(src) .xenium_transcript_disk(fx$paths$transcripts,
        src, qv_threshold = 20, output = "store", threshold = 5))
    .expect_same_store(b)
})

test_that("MERSCOPE transcripts (csv and parquet, FOV filter) match", {
    .skip_no_ingest()
    for (fmt in c("csv", "parquet")) {
        d <- .mk_merscope_disk_dir(tx_format = fmt)
        tx <- file.path(d, paste0("detected_transcripts.", fmt))
        withr::local_options(giottodisk.artifact_dump = withr::local_tempdir())
        for (fovs in list(NULL, 1L)) {
            b <- .both_paths(function(src) .merscope_transcript_disk(tx, src,
                FOVs = fovs, dropcols = c("x", "y", "barcode_id"),
                output = "store", threshold = 4))
            .expect_same_store(b)
        }
    }
})

test_that("bin points (a plain parquet source) match", {
    .skip_no_ingest()
    set.seed(3)
    f <- tempfile(fileext = ".parquet")
    arrow::write_parquet(arrow::arrow_table(feat_ID = 1:2000,
        x = stats::runif(2000, 0, 80), y = -stats::runif(2000, 0, 40)), f)
    write <- function(src) {
        qs <- as(fileStore(path = f, read_fun = .overlap_read_fun), "queryableStore")
        qs@params$points_spec <- .ingest_points_spec(f)
        sourceWrite(src, qs, store_type = "parquetGeomTile", type = "points",
            id_col = "feat_ID", sdimx = "x", sdimy = "y", threshold = 250)
    }
    .expect_same_store(.both_paths(write))
})

test_that("a given plan with an empty tile, and an `i` subset, match", {
    .skip_no_ingest()
    f <- .tx_parquet()
    qs_for <- function() {
        qs <- as(fileStore(path = f, read_fun = .overlap_read_fun), "queryableStore")
        qs@params$points_spec <- .ingest_points_spec(f)
        qs
    }
    fp <- tilework::freeTilePlan()
    fp$bounds <- rbind(c(0, 40, 0, 50), c(40, 41, -10, -5), c(41, 100, 0, 50))
    fp$n_records <- c(1, 0, 1)
    for (i in list(NULL, c(1L, 3L), 3L)) {
        b <- .both_paths(function(src) sourceWrite(src, qs_for(),
            store_type = "parquetGeomTile", type = "points", tiles = fp, i = i,
            id_col = "feature_name", sdimx = "x_location", sdimy = "y_location"))
        .expect_same_store(b)
    }
})

test_that("points on one vertical line (zero-width extent) match", {
    .skip_no_ingest()
    f <- tempfile(fileext = ".parquet")
    arrow::write_parquet(data.frame(feat_ID = 1:50, x = 5, y = seq(0, 49) + 0.5), f)
    write <- function(src) {
        qs <- as(fileStore(path = f, read_fun = .overlap_read_fun), "queryableStore")
        qs@params$points_spec <- .ingest_points_spec(f)
        sourceWrite(src, qs, store_type = "parquetGeomTile", type = "points",
            id_col = "feat_ID", sdimx = "x", sdimy = "y")
    }
    b <- .both_paths(write)
    expect_identical(nrow(b$rs@tiles@bounds), 1L)
    .expect_same_store(b)
})

test_that("Atera v1 zarr transcripts: the Rust decoder writes the R store", {
    .skip_no_ingest()
    skip_if_not(.gd_ingest("points_atera_zarr"))
    av6 <- testthat::test_path("fixtures", "atera_v6", "transcripts.zarr.zip")
    skip_if_not(file.exists(av6))
    write <- function(src) sourceWrite(src, GiottoDisk:::tenxZarrTxInput(av6),
        store_type = "parquetGeomTile", threshold = 100, qv_threshold = 30,
        workers = 1L, window_rows = 300L, buffer_rows = 500L)
    b <- .both_paths(write)
    expect_identical(length(b$r@tiles@metadata$n_records), 16L)
    .expect_same_store(b)
    # part files: the same per tile (several flushes per tile here)
    for (d in basename(.tile_dirs(b$r))) {
        expect_identical(
            sort(list.files(file.path(.idpath(b$rs@path, b$rs@uid), d))),
            sort(list.files(file.path(.idpath(b$r@path, b$r@uid), d))))
    }
})

# gate ####

test_that("a missing GiottoIngest is not installed, silently (private repo)", {
    fake <- paste0("NoIngest", sample.int(1e6, 1L))
    expect_no_message(r <- .gd_ingest_installed(fake))
    expect_false(r)
})

test_that("giottodisk.use_ingest = FALSE skips the install check", {
    local_mocked_bindings(.gd_ingest_installed = function(...) stop("checked"))
    withr::local_options(giottodisk.use_ingest = FALSE)
    expect_false(.gd_ingest("points_parquet"))
})

test_that("without the capability, or for writes Rust does not reproduce, R writes", {
    f <- .tx_parquet(200L)
    qs <- as(fileStore(path = f, read_fun = .overlap_read_fun), "queryableStore")
    qs@params$points_spec <- .ingest_points_spec(f)
    st <- storeCreate(path = tempfile(), type = "parquetGeomTile")
    called <- FALSE
    local_mocked_bindings(.gd_ingest = function(cap) { called <<- TRUE; FALSE })
    expect_null(.store_write_points_ingest(st, qs, NULL, NULL, NULL, "feature_name",
        "x_location", "y_location", TRUE, FALSE, list(), NULL))
    expect_true(called)
    # write_param, dry_run, non-contiguous, another codec: no gate check at all
    called <- FALSE
    expect_null(.store_write_points_ingest(st, qs, NULL, NULL, NULL, "feature_name",
        "x_location", "y_location", TRUE, FALSE, list(chunk_size = 10), NULL))
    expect_null(.store_write_points_ingest(st, qs, NULL, NULL, NULL, "feature_name",
        "x_location", "y_location", TRUE, TRUE, list(), NULL))
    expect_null(.store_write_points_ingest(st, qs, NULL, NULL, NULL, "feature_name",
        "x_location", "y_location", FALSE, FALSE, list(), NULL))
    withr::with_options(list(giottodisk.parquet_compression = "zstd"),
        expect_null(.store_write_points_ingest(st, qs, NULL, NULL, NULL,
            "feature_name", "x_location", "y_location", TRUE, FALSE, list(), NULL)))
    expect_false(called)
})

test_that("input GiottoIngest reports as unsupported is written in R", {
    .skip_no_ingest()
    df <- data.frame(x = stats::runif(300, 0, 10), y = stats::runif(300, 0, 10),
        gene = "g", when = as.POSIXct("2020-01-01", tz = "UTC") + 1:300)
    f <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, f)
    write <- function(src) {
        qs <- as(fileStore(path = f, read_fun = .overlap_read_fun), "queryableStore")
        qs@params$points_spec <- .ingest_points_spec(f)
        sourceWrite(src, qs, store_type = "parquetGeomTile", type = "points",
            id_col = "gene", sdimx = "x", sdimy = "y")
    }
    st <- withr::with_options(list(giottodisk.use_ingest = TRUE),
        write(gDirSource(path = withr::local_tempdir())))
    expect_identical(nrow(st), 300)
    expect_true("when" %in% names(storeRead(st)))
})

test_that("a source with every row filtered out writes no store, as in R", {
    .skip_no_ingest()
    f <- .tx_parquet(300L)
    b <- .both_paths(function(src) {
        tryCatch(.xenium_transcript_disk(f, src, qv_threshold = 1e6, output = "store"),
            error = function(e) conditionMessage(e))
    }, rust = FALSE)
    expect_identical(b$rs, b$r)
})
