# Optional Rust ingest backend (GiottoIngest). Like the kernels: a caller asks
# .gd_ingest(cap) and runs its own R path when the answer is FALSE, so the
# package writes the same stores whether or not GiottoIngest is installed.
# `giottodisk.use_ingest = FALSE` forces the R paths.
.gd_ingest <- function(cap) {
    isTRUE(getOption("giottodisk.use_ingest", TRUE)) &&
        .gd_ingest_installed() &&
        GiottoIngest::has_capability(cap)
}

# A missing or outdated GiottoIngest is mentioned once per session, since the
# R path is only slower. `pkg` is a seam for the tests.
.gd_ingest_installed <- function(pkg = "GiottoIngest") {
    GiottoUtils::package_check(pkg,
        repository = sprintf("github:giotto-suite/%s>=0.1.0", pkg),
        on_missing = "message",
        once = TRUE,
        custom_msg = sprintf(paste0(
            "%s is not installed, so tiled transcript writes run in R: same ",
            "stores, slower and with more memory. It is installed from ",
            "GitHub (giotto-suite/%s) and needs a Rust toolchain."
        ), pkg, pkg)
    )
}

# Rust threads are light (one process, shared memory), so more of them than
# R workers: the session's workers, or up to 8 physical cores.
.gd_ingest_threads <- function() {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores)) cores <- 1L
    max(.par_workers(), min(8L, as.integer(cores)))
}

# The Rust writers use snappy, the default, without a level.
.gd_ingest_codec_ok <- function(compression = .parquet_compression()) {
    identical(compression, "snappy") && is.null(.parquet_compression_level())
}

# A point source described for GiottoIngest: the files and query a reader's
# read_fun applies (filters, y flip, string casts, dropped columns), set on
# `qs@params$points_spec` by the readers that have one.
.ingest_points_spec <- function(path, flip_y = FALSE, filters = list(),
    to_string = character(), drop = character()) {
    files <- if (dir.exists(path)) {
        sort(list.files(path, pattern = "\\.parquet$", recursive = TRUE,
            full.names = TRUE), method = "radix")
    } else {
        path
    }
    list(files = files, flip_y = isTRUE(flip_y), filters = filters,
        to_string = as.character(to_string), drop = as.character(drop))
}

# storeWrite(parquetGeomTileStore, queryableStore) in Rust: the same plan
# (tilework's quadtreePlan, counted by GiottoIngest), tiles and metadata.
# NULL when the R path should run instead: no spec, gate off, a request the
# Rust writer does not reproduce, or input it reports as unsupported.
.store_write_points_ingest <- function(store, data, threshold, tiles, i,
    id_col, sdimx, sdimy, contiguous, dry_run, write_param, verbose) {
    spec <- data@params$points_spec
    if (is.null(spec) || !isTRUE(contiguous) || isTRUE(dry_run) ||
        length(write_param) > 0L || !.gd_ingest_codec_ok() ||
        !.gd_ingest("points_parquet")) {
        return(NULL)
    }
    threads <- .gd_ingest_threads()
    tryCatch({
        src <- GiottoIngest::gi_source(spec$files, sdimx, sdimy,
            flip_y = spec$flip_y, filters = spec$filters)
        ex <- GiottoIngest::gi_extent(src, threads)
        if (ex$n == 0) return(list(empty = TRUE))
        threshold <- threshold %||% .auto_threshold(ex$n, type = "points")
        if (is.null(tiles)) {
            tiles <- .points_seed_tiles(terra::ext(ex$ext))
        }
        plans <- .tile_plan_quadtree(store, src,
            threshold = threshold, tiles = tiles, i = i, dry_run = FALSE,
            verbose = verbose, qtree_args = list(threads = threads))
        store <- plans$store
        store@geomtype <- "points"
        root <- .idpath(store@path, store@uid)
        GiottoIngest::gi_write_point_tiles(src,
            store@tiles@bounds[plans$tile_indices, , drop = FALSE],
            tile_dirs = basename(.tilepath(root, plans$tile_indices)),
            id_col = id_col, to_string = spec$to_string, drop = spec$drop,
            root = root, geo_json = .arrow_meta_add_geoparquet(store)$geo,
            threads = threads)
        store@params$crs <- ""
        vmsg(.v = verbose, sprintf(
            "[storeWrite] %.0f points into %d tiles (GiottoIngest)",
            ex$n, length(plans$tile_indices)))
        initialize(store)
    }, giottoingest_unsupported = function(e) {
        vmsg(.v = verbose, "[storeWrite]", conditionMessage(e),
            "; writing in R")
        NULL
    })
}

# The Atera grid cells in Rust (GiottoIngest::gi_atera_tx_cells), on the
# same plan, with the same per-cell results as .atera_tx_cell_worker().
# NULL when the R workers should run instead.
.atera_tx_cells_ingest <- function(zarr_path, cells, plan, root,
    qv_threshold, flip, window_rows, buffer_rows, compression, dropcols,
    verbose = NULL) {
    if (length(dropcols) > 0L || !.gd_ingest_codec_ok(compression) ||
        !.gd_ingest("points_atera_zarr")) {
        return(NULL)
    }
    src <- .zarr_open(zarr_path)
    lk <- .atera_tx_lookups(src)
    .zarr_close(src)
    sz <- .atera_tx_sizes(.gd_ingest_threads(), length(cells), window_rows,
        buffer_rows)
    tryCatch(
        GiottoIngest::gi_atera_tx_cells(zarr_path, cells,
            ords = match(cells, unique(plan$tile)) - 1L,
            s = vapply(cells, function(t) plan$s[plan$tile == t][1L],
                numeric(1L)),
            tile_index = lapply(cells, function(t)
                plan$tile_index[plan$tile == t]),
            gene_names = lk$gene_names, cw_gene = lk$cw_gene, root = root,
            qv_threshold = qv_threshold, flip = flip,
            window_rows = sz$window_rows, buffer_rows = sz$buffer_rows,
            threads = sz$workers),
        giottoingest_unsupported = function(e) {
            vmsg(.v = verbose, "[zarr]", conditionMessage(e),
                "; decoding in R")
            unlink(root, recursive = TRUE)
            NULL
        })
}
