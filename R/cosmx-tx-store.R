#' @include methods-storeWrite.R
NULL

# CosMx transcripts -> parquetGeomTileStore, one FOV per worker; each FOV is
# split into s x s tiles as the Atera grid cells are (`.slot_tile_plan()`).

#' @title CosMx transcripts input
#' @name cosmxTxInput-class
#' @description
#' Marks a FOV-partitioned parquet copy of a CosMx tx_file as input to
#' [sourceWrite()] / [storeWrite()] for a `parquetGeomTileStore` of points.
#' Writing streams each FOV straight into its tile partitions.
#' @keywords internal
setClass("cosmxTxInput", contains = "fileStore")

#' @rdname cosmxTxInput-class
#' @param path directory of `fov=<k>/` parquet partitions
#' @keywords internal
cosmxTxInput <- function(path) {
    if (!.cosmx_fov_partitioned(path)) {
        stop("[cosmxTxInput] ", path, " has no fov=<k> partitions",
            call. = FALSE)
    }
    new("cosmxTxInput", path = normalizePath(path))
}

.cosmx_fov_partitioned <- function(path) {
    dir.exists(path) && any(startsWith(
        basename(list.dirs(path, recursive = FALSE)), "fov="))
}

.cosmx_tx_fov_files <- function(path, fovs = NULL) {
    d <- list.dirs(path, recursive = FALSE)
    d <- d[startsWith(basename(d), "fov=")]
    if (!is.null(fovs)) {
        fov <- suppressWarnings(as.integer(sub("^fov=", "", basename(d))))
        d <- d[fov %in% fovs]
    }
    stats::setNames(lapply(d, list.files, pattern = "[.]parquet$",
        full.names = TRUE), basename(d))
}

#' @keywords internal
#' @noRd
.cosmx_tx_fov_worker <- function(files, group, plan, root, uid,
    window_rows, buffer_rows, compression, dropcols = character(0L),
    sdimx = "x_global_px", sdimy = "y_global_px", id_col = "target") {
    sub <- plan[plan$tile == group, , drop = FALSE]
    n <- attr(plan, "rows")[[group]]
    fov <- suppressWarnings(as.integer(sub("^fov=", "", group)))
    readers <- lapply(files, arrow::ParquetFileReader$create)
    groups <- do.call(rbind, lapply(seq_along(readers), function(f) {
        g <- readers[[f]]$num_row_groups
        if (g > 0) cbind(f, seq_len(g) - 1L)
    }))
    nm <- names(readers[[1L]]$GetSchema())
    keep <- setdiff(nm, c(sdimx, sdimy, id_col, dropcols))
    with_fov <- !"fov" %in% c(nm, dropcols)

    # split from the FOV's locations: all of them when they fit a window,
    # else every k-th row group
    step <- max(1L, ceiling(n / window_rows))
    ixy <- match(c(sdimx, sdimy), nm) - 1L
    smp <- lapply(seq(1L, nrow(groups), by = step), function(r) {
        readers[[groups[r, 1L]]]$ReadRowGroup(groups[r, 2L], ixy)
    })
    sx <- unlist(lapply(smp, function(t) as.vector(t[[sdimx]])))
    sy <- unlist(lapply(smp, function(t) as.vector(t[[sdimy]])))
    sp <- .quantile_split(sx, sy, sub$s[1L])
    bounds <- .quantile_bounds(sp, c(min(sx), max(sx), min(sy), max(sy)))
    rm(smp, sx, sy)

    counts <- numeric(nrow(sub))
    ext <- matrix(c(Inf, -Inf, Inf, -Inf), nrow(sub), 4L, byrow = TRUE)
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
            meta <- tb[, c(id_col, keep)]
            names(meta)[1L] <- "feat_ID"
            if (with_fov) {
                meta$fov <- arrow::Array$create(rep(fov, length(x)),
                    type = arrow::int32())
                meta <- meta[, c("feat_ID", "fov", keep)]
            }
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

    route <- function(tbl) {
        x <- as.vector(tbl[[sdimx]])
        y <- as.vector(tbl[[sdimy]])
        key <- .quantile_slot(sp, x, y)
        # group rows by tile: one reorder, then zero-copy slices
        o <- order(key, na.last = NA)
        if (!length(o)) return(invisible())
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
            buf[[i]][[length(buf[[i]]) + 1L]] <<- list(
                tbl = tbl$Slice(starts[r] - 1L, runs$lengths[r]),
                x = x[idx], y = y[idx]
            )
        }
        buffered <<- buffered + length(o)
        if (buffered >= buffer_rows) flush()
    }

    window <- list()
    window_n <- 0
    for (r in seq_len(nrow(groups))) {
        t <- readers[[groups[r, 1L]]]$ReadRowGroup(groups[r, 2L])
        window[[length(window) + 1L]] <- t
        window_n <- window_n + t$num_rows
        if (window_n >= window_rows || r == nrow(groups)) {
            route(do.call(arrow::concat_tables, window))
            window <- list()
            window_n <- 0
        }
    }
    flush()

    hit <- counts > 0
    bounds[hit, 1L] <- pmin(bounds[hit, 1L], ext[hit, 1L])
    bounds[hit, 2L] <- pmax(bounds[hit, 2L], ext[hit, 2L])
    bounds[hit, 3L] <- pmin(bounds[hit, 3L], ext[hit, 3L])
    bounds[hit, 4L] <- pmax(bounds[hit, 4L], ext[hit, 4L])
    list(tile_index = sub$tile_index, counts = counts, bounds = bounds)
}

#' @rdname storeWrite-parquetGeomTileStore
#' @param fovs CosMx FOVs to write (`NULL` writes all)
#' @export
setMethod(
    "storeWrite",
    signature("parquetGeomTileStore", "cosmxTxInput"),
    function(store, data, threshold = NULL, fovs = NULL,
        workers = .par_workers(), window_rows = NULL, buffer_rows = NULL,
        compression = .parquet_compression(), dropcols = character(0L),
        verbose = NULL, ...) {
        vmsg(.is_debug = TRUE,
            "[storeWrite] parquetGeomTileStore,cosmxTxInput")
        # parquet reads hang with fewer than 2 arrow IO threads
        io <- arrow::io_thread_count()
        if (io < 2L) {
            arrow::set_io_thread_count(2L)
            on.exit(suppressWarnings(arrow::set_io_thread_count(io)),
                add = TRUE)
        }
        files <- .cosmx_tx_fov_files(data@path, fovs)
        rows <- vapply(files, function(f) {
            sum(vapply(f, function(p) {
                as.numeric(arrow::ParquetFileReader$create(p)$num_rows)
            }, numeric(1L)))
        }, numeric(1L))
        rows <- rows[rows > 0]
        if (!length(rows)) return(NULL)
        threshold <- threshold %||% .auto_threshold(sum(rows), type = "points")
        plan <- .slot_tile_plan(rows, threshold)

        root <- .idpath(store@path, store@uid)
        # largest FOVs first, so they do not finish last
        groups <- names(sort(rows, decreasing = TRUE))
        sz <- .atera_tx_sizes(workers, length(groups), window_rows,
            buffer_rows)
        run_one <- function(group) {
            .cosmx_tx_fov_worker(files[[group]], group, plan, root,
                store@uid, window_rows = sz$window_rows,
                buffer_rows = sz$buffer_rows, compression = compression,
                dropcols = as.character(dropcols))
        }
        res <- if (.Platform$OS.type == "unix") {
            .mclapply_fresh(groups, run_one, sz$workers)
        } else if (sz$workers > 1L) {
            GiottoUtils::lapply_flex(groups, run_one, cores = sz$workers,
                future.seed = NULL)
        } else {
            lapply(groups, run_one)
        }
        errs <- vapply(res, inherits, logical(1L), "try-error")
        if (any(errs)) {
            stop("[storeWrite] CosMx transcript worker failed: ",
                attr(res[[which(errs)[1L]]], "condition")$message,
                call. = FALSE)
        }

        store <- .slot_tiles_store(store, res, nrow(plan))
        if (is.null(store)) return(NULL)
        counts <- store@tiles@metadata$n_records
        vmsg(.v = verbose, sprintf(
            "[storeWrite] %.0f transcripts into %d tiles (%d FOVs)",
            sum(counts), sum(counts > 0), length(groups)))
        initialize(store)
    }
)
