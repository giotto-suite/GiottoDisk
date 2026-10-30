# atera-zarr-v6 transcripts -> parquetGeomTileStore, without an
# intermediate table.
#
# The archive stores transcripts in a regular spatial grid
# (`grid/<i,j>`, `grid_size` microns per cell), gene-sorted within each
# cell. Each grid cell is split into s x s store tiles, s sized from the
# cell's row count so a tile holds about `threshold` points. Tissue rarely
# fills a cell, so the split is data-driven: x quantile bands, then y
# quantile cells within each band, cut from a strided sample of the cell's
# locations. One worker owns one grid cell: it samples, streams the cell
# in row windows, routes rows to tile buffers and appends them to
# `tile_index=<k>/` part files, keeping `row_index` contiguous per tile.
# Every row is read and written once.

#' @title 10x Atera zarr transcripts input
#' @name tenxZarrTxInput-class
#' @description
#' Marks an Atera v1 `transcripts.zarr.zip` archive as input to
#' [sourceWrite()] / [storeWrite()] for a `parquetGeomTileStore` of
#' points. Writing streams the archive straight into tile partitions.
#' @keywords internal
setClass("tenxZarrTxInput", contains = "fileStore")

#' @rdname tenxZarrTxInput-class
#' @param zarr_path path to `transcripts.zarr.zip` (or unzipped `.zarr`)
#' @keywords internal
tenxZarrTxInput <- function(zarr_path) {
    layout <- .zarr_detect_product(zarr_path, "transcripts")
    if (!identical(layout, "atera-zarr-v6")) {
        stop("[tenxZarrTxInput] ", zarr_path, " is a ", layout,
            " archive; only atera-zarr-v6 is streamed to a tile store",
            call. = FALSE)
    }
    new("tenxZarrTxInput", path = normalizePath(zarr_path))
}

# Tile slots for every grid cell: s x s per cell, s from the cell's row
# count. `tile_index` is the row number; bounds are filled in by the
# workers, which cut them from the data.
.atera_tx_tile_plan <- function(src, threshold) {
    tiles <- .atera_tx_tiles(src)
    n <- vapply(tiles, function(t) {
        unlist(.zarr_meta(src, paste0("grid/", t, "/location"))$shape)[1L]
    }, numeric(1L))
    s <- if (is.finite(threshold)) pmax(1, ceiling(sqrt(n / threshold)))
        else rep(1, length(n))
    plan <- data.frame(
        tile = rep(tiles, s^2),
        s = rep(s, s^2),
        slot = unlist(lapply(s, function(k) seq_len(k^2)))
    )
    plan$tile_index <- seq_len(nrow(plan))
    attr(plan, "rows") <- stats::setNames(n, tiles)
    plan
}

# Quantile cuts for an s x s split of points (x, y): `xc` holds the s - 1
# interior x cuts, `yc[[b]]` the s - 1 interior y cuts of x band b.
.quantile_split <- function(x, y, s) {
    probs <- seq_len(s - 1) / s
    xc <- unique(stats::quantile(x, probs, names = FALSE, type = 1))
    band <- findInterval(x, xc) + 1L
    yc <- lapply(seq_len(s), function(b) {
        yb <- y[band == b]
        if (length(yb)) unique(stats::quantile(yb, probs, names = FALSE,
            type = 1)) else numeric(0)
    })
    list(xc = xc, yc = yc, s = s)
}

# Slot (1-based, band-major) of each point under a split.
.quantile_slot <- function(sp, x, y) {
    band <- findInterval(x, sp$xc) + 1L
    cell <- integer(length(x))
    for (b in unique(band)) {
        hit <- band == b
        cell[hit] <- findInterval(y[hit], sp$yc[[b]]) + 1L
    }
    (band - 1L) * sp$s + cell
}

# Bounds of each slot of a split over the extent `e` (xmin, xmax, ymin,
# ymax), as a slots x 4 matrix.
.quantile_bounds <- function(sp, e) {
    s <- sp$s
    xb <- c(e[1L], sp$xc, e[2L])
    out <- matrix(NA_real_, s * s, 4L)
    for (b in seq_len(s)) {
        if (b > length(xb) - 1L) next
        yb <- c(e[3L], sp$yc[[b]], e[4L])
        for (k in seq_len(length(yb) - 1L)) {
            out[(b - 1L) * s + k, ] <- c(xb[b], xb[b + 1L], yb[k], yb[k + 1L])
        }
    }
    out
}

# Strided sample of a grid cell's locations: the whole cell when small,
# else `n_ranges` evenly spaced runs. Rows are gene-sorted, so spreading
# the runs samples many genes.
.atera_tx_sample_xy <- function(src, tile, n, sample_rows = 1048576,
    n_ranges = 16L) {
    pref <- paste0("grid/", tile, "/location")
    if (n <= sample_rows) {
        loc <- .zarr_columns(src, pref, range = c(1, n))
        return(list(x = loc[[1L]], y = loc[[2L]]))
    }
    len <- floor(sample_rows / n_ranges)
    starts <- floor(seq(1, n - len + 1, length.out = n_ranges))
    parts <- lapply(starts, function(r0) {
        .zarr_columns(src, pref, range = c(r0, r0 + len - 1))[1:2]
    })
    list(x = unlist(lapply(parts, `[[`, 1L)),
        y = unlist(lapply(parts, `[[`, 2L)))
}

# One worker: stream one grid cell into its tile partitions. Returns the
# tile indices, row counts and bounds (y flipped when `flip`).
# Mirrored in GiottoIngest (Rust): a change here needs its parity tests.
#' @keywords internal
#' @noRd
.atera_tx_cell_worker <- function(zarr_path, tile, plan, root, uid,
    qv_threshold, flip, window_rows, buffer_rows, compression,
    dropcols = character(0L)) {
    src <- .zarr_open(zarr_path)
    on.exit(.zarr_close(src), add = TRUE)
    lk <- .atera_tx_lookups(src)
    sub <- plan[plan$tile == tile, , drop = FALSE]
    n <- attr(plan, "rows")[[tile]]
    ord <- match(tile, unique(plan$tile)) - 1L

    smp <- .atera_tx_sample_xy(src, tile, n)
    sp <- .quantile_split(smp$x, smp$y, sub$s[1L])
    # raw-space bounds from the sample extent; widened below by the data
    bounds <- .quantile_bounds(sp,
        c(min(smp$x), max(smp$x), min(smp$y), max(smp$y)))

    counts <- numeric(nrow(sub))
    ext <- matrix(c(Inf, -Inf, Inf, -Inf), nrow(sub), 4L, byrow = TRUE)
    # per tile: list of (arrow slice, x, y) pieces awaiting a write
    buf <- vector("list", nrow(sub))
    buffered <- 0
    n_part <- 0L

    flush <- function() {
        for (i in which(lengths(buf) > 0L)) {
            pieces <- buf[[i]]
            buf[[i]] <<- list()
            n_part <<- n_part + 1L
            tb <- do.call(arrow::concat_tables, lapply(pieces, `[[`, "tbl"))
            x <- unlist(lapply(pieces, `[[`, "x"))
            y <- unlist(lapply(pieces, `[[`, "y"))
            ext[i, ] <<- c(min(ext[i, 1L], x), max(ext[i, 2L], x),
                min(ext[i, 3L], y), max(ext[i, 4L], y))
            if (flip) y <- -y
            keep <- setdiff(names(tb), c("x_location", "y_location",
                setdiff(dropcols, "feature_name")))
            meta <- tb[, c("feature_name", setdiff(keep, "feature_name"))]
            names(meta)[1L] <- "feat_ID"
            tile_store <- storeCreate(
                path = .tilepath(root, sub$tile_index[i]),
                type = "parquetGeom", uid = uid
            )
            .store_write_points_arrow(tile_store, x = x, y = y, meta = meta,
                row_offset = counts[i], uid_partition = FALSE,
                compression = compression,
                basename_template = sprintf("part-%05d-{i}.parquet", n_part)
            )
            counts[i] <<- counts[i] + length(x)
        }
        buffered <<- 0
    }

    for (r0 in seq.int(1, n, by = window_rows)) {
        u <- list(tile = tile, ord = ord, r0 = r0,
            r1 = min(r0 + window_rows - 1, n))
        tbl <- .atera_tx_unit_table(src, u, lk)
        tbl$transcript_id <- tbl$transcript_id$cast(arrow::utf8())
        x <- as.vector(tbl$x_location)
        y <- as.vector(tbl$y_location)
        key <- .quantile_slot(sp, x, y)
        if (!is.null(qv_threshold)) {
            key[as.vector(tbl$qv) < qv_threshold] <- NA_integer_
        }
        # group rows by tile: one reorder, then zero-copy slices
        o <- order(key, na.last = NA)
        if (!length(o)) next
        tbl <- tbl$Take(o - 1L)
        key <- key[o]
        x <- x[o]
        y <- y[o]
        runs <- rle(key)
        ends <- cumsum(runs$lengths)
        starts <- ends - runs$lengths + 1L
        for (r in seq_along(starts)) {
            idx <- starts[r]:ends[r]
            i <- runs$values[r]
            buf[[i]][[length(buf[[i]]) + 1L]] <- list(
                tbl = tbl$Slice(starts[r] - 1L, runs$lengths[r]),
                x = x[idx], y = y[idx]
            )
        }
        buffered <- buffered + length(o)
        if (buffered >= buffer_rows) flush()
    }
    flush()

    hit <- counts > 0
    bounds[hit, 1L] <- pmin(bounds[hit, 1L], ext[hit, 1L])
    bounds[hit, 2L] <- pmax(bounds[hit, 2L], ext[hit, 2L])
    bounds[hit, 3L] <- pmin(bounds[hit, 3L], ext[hit, 3L])
    bounds[hit, 4L] <- pmax(bounds[hit, 4L], ext[hit, 4L])
    if (flip) bounds[, 3:4] <- -bounds[, 4:3]
    list(tile_index = sub$tile_index, counts = counts, bounds = bounds)
}

# Per-worker memory of the transcripts writer (RSS above the session),
# measured on the largest Atera kidney grid cell: a fixed part (arrow pool,
# decode buffers, allocator retention) plus bytes per decoded window row
# and per buffered row.
.atera_tx_mem <- list(fixed = 2.1 * 2^30, window = 355, buffer = 180)

# Workers, window and buffer within the ingest budget: as many workers as
# fit at the smallest window (1M rows) and buffer (2M rows), then the
# rest of each worker's share split evenly between window and buffer.
# Explicit `window_rows` / `buffer_rows` are kept.
.atera_tx_sizes <- function(workers, n_cells, window_rows = NULL,
    buffer_rows = NULL) {
    m <- .atera_tx_mem
    workers <- max(1L, min(as.integer(workers), as.integer(n_cells)))
    budget <- .ingest_budget()
    if (is.na(budget)) budget <- workers * (m$fixed + 2^31)
    least <- m$fixed + 1e6 * m$window + 2e6 * m$buffer
    workers <- max(1L, min(workers, floor(budget / least)))
    spare <- max(0, budget / workers - m$fixed)
    list(
        workers = workers,
        window_rows = window_rows %||%
            max(1e6, floor(spare / 2 / m$window)),
        buffer_rows = buffer_rows %||%
            max(2e6, floor(spare / 2 / m$buffer))
    )
}

# parallel::mclapply, forking a fresh child per element even with one
# worker (mclapply runs in the session itself when mc.cores = 1).
.mclapply_fresh <- function(X, FUN, workers) {
    if (workers > 1L) {
        return(parallel::mclapply(X, FUN, mc.cores = workers,
            mc.preschedule = FALSE))
    }
    lapply(X, function(x) {
        job <- parallel::mcparallel(FUN(x))
        parallel::mccollect(job)[[1L]]
    })
}

#' @rdname storeWrite-parquetGeomTileStore
#' @param qv_threshold minimum quality score kept (`NULL` keeps all)
#' @param flip_vertical negate y, as the Xenium/Atera readers do
#' @param workers parallel workers, one grid cell each
#' @param window_rows transcripts decoded per read. Default: derived
#'   from the memory budget (`giottodisk.chunk_ram_frac`) and `workers`
#' @param buffer_rows transcripts buffered per worker before writing; a
#'   grid cell is written as about (cell rows / `buffer_rows`) files per
#'   tile. Default: derived like `window_rows`
#' @param dropcols attribute columns not to write
#' @export
setMethod(
    "storeWrite",
    signature("parquetGeomTileStore", "tenxZarrTxInput"),
    function(store, data, threshold = NULL, qv_threshold = 20,
        flip_vertical = TRUE, workers = .par_workers(),
        window_rows = NULL, buffer_rows = NULL,
        compression = .parquet_compression(), dropcols = character(0L),
        verbose = NULL, ...) {
        vmsg(.is_debug = TRUE,
            "[storeWrite] parquetGeomTileStore,tenxZarrTxInput")
        src <- .zarr_open(data@path)
        n_total <- as.numeric(.zarr_attrs(src, "")$number_rnas)
        threshold <- threshold %||% .auto_threshold(n_total, type = "points")
        plan <- .atera_tx_tile_plan(src, threshold)
        .zarr_close(src)

        root <- .idpath(store@path, store@uid)
        # largest grid cells first, so they do not finish last
        cells <- names(sort(attr(plan, "rows"), decreasing = TRUE))
        # sizes as asked, for the Rust path's own sizing below
        window_arg <- window_rows
        buffer_arg <- buffer_rows
        sz <- .atera_tx_sizes(workers, length(cells), window_rows,
            buffer_rows)
        workers <- sz$workers
        window_rows <- sz$window_rows
        buffer_rows <- sz$buffer_rows
        run_one <- function(cell) {
            .atera_tx_cell_worker(data@path, cell, plan, root, store@uid,
                qv_threshold = qv_threshold, flip = flip_vertical,
                window_rows = window_rows, buffer_rows = buffer_rows,
                compression = compression, dropcols = as.character(dropcols))
        }
        # Rust path (GiottoIngest): every grid cell on threads in this
        # process; otherwise, on unix every grid cell runs in a forked child,
        # also with one worker, so each cell's memory is returned when it
        # finishes
        res <- .atera_tx_cells_ingest(data@path, cells, plan, root,
            qv_threshold = qv_threshold, flip = flip_vertical,
            window_rows = window_arg, buffer_rows = buffer_arg,
            compression = compression, dropcols = as.character(dropcols),
            verbose = verbose)
        if (is.null(res)) res <- if (.Platform$OS.type == "unix") {
            .mclapply_fresh(cells, run_one, workers)
        } else if (workers > 1L) {
            GiottoUtils::lapply_flex(cells, run_one, cores = workers,
                future.seed = NULL)
        } else {
            lapply(cells, run_one)
        }
        errs <- vapply(res, inherits, logical(1L), "try-error")
        if (any(errs)) {
            stop("[zarr] transcript worker failed: ",
                attr(res[[which(errs)[1L]]], "condition")$message,
                call. = FALSE)
        }

        counts <- numeric(nrow(plan))
        bounds <- matrix(NA_real_, nrow(plan), 4L,
            dimnames = list(NULL, c("xmin", "xmax", "ymin", "ymax")))
        for (r in res) {
            counts[r$tile_index] <- r$counts
            bounds[r$tile_index, ] <- r$bounds
        }
        # an empty slot (fewer distinct cuts than planned) takes a
        # degenerate box at the store origin; it holds no rows
        bounds[is.na(bounds)] <- 0
        if (!any(counts > 0)) return(NULL)
        fp <- methods::new("freeTilePlan", bounds = bounds)
        fp@metadata <- data.frame(tile = seq_along(counts), n_records = counts)
        store@tiles <- fp
        nz <- bounds[counts > 0, , drop = FALSE]
        .pstore_disk_extent(store) <- c(min(nz[, 1L]), max(nz[, 2L]),
            min(nz[, 3L]), max(nz[, 4L]))
        store@geomtype <- "points"
        store@params$crs <- ""
        vmsg(.v = verbose, sprintf(
            "[storeWrite] %.0f transcripts into %d tiles (%d grid cells)",
            sum(counts), sum(counts > 0), length(cells)))
        initialize(store)
    }
)
