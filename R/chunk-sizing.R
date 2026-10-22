#' @include class-parquetExprStore.R
NULL

# chunk-sizing ####
# How a streaming pass decides how much to read at once.
#
# A window is derived per read rather than stored on the object: it depends on
# free RAM, which is a property of the machine doing the reading, so baking one
# into a store at write time made it stale the moment the store moved. What IS
# machine-independent is the payload size, and that is what `@stats` caches.
#
# The model is rows x payload-per-row against a fraction of free RAM. Nothing
# in it is matrix-specific -- a tabular store can use it by supplying its own
# shape and `bytes_per_nz`.
#
# Contents:
#   .recommend_chunk_size   the sizing computation
#   .chunk_bytes_per_row    payload per row, shared with the report
#   .chunk_budget           RAM budget less the persistent PCA sketch
#   .sc_free_ram            per-OS free-RAM probe
#   storeChunkInfo          exported report over both read shapes
#
# Options that steer it:
#   giottodisk.chunk_ram_frac   fraction of free RAM to budget (default 0.25);
#                               also sizes the DuckDB memory limit
#   giottodisk.chunk_size       pin an absolute window; also the fallback
#                               when the derivation cannot run
#   giottodisk.chunk_min_rows   smallest derived window, in cells (10000)
#   giottodisk.chunk_min_frac   smallest budget, as a fraction of free RAM
#                               (0.05)
#   giottodisk.gef_batch_rows   pin the records per Stereo-seq GEF ingest
#                               batch (`.gef_batch_rows()`)
#   giottodisk.h5_batch_rows    pin the records per 10x .h5 ingest batch
#                               (`.tenxh5_ranges()`)

# Recommend a streaming chunk size from free RAM.
#
# Internal: the value has no user-facing destination -- there is no
# `chunk_size` slot or constructor argument to pass it to. Users steer the
# window through `giottodisk.chunk_size` (pin an absolute value) and
# `giottodisk.chunk_ram_frac` (scale the budget); this is what those options
# steer. `.pe_window_cells()` is the only caller.
#
# Estimates the largest window, in rows, whose materialized payload stays
# inside a target fraction of currently free RAM.
#
#   n_cells       rows in the current view (cells, for an expression store)
#   n_genes       columns in the current view
#   density       fill fraction of the view
#   k             PCA dimensions; reserves the persistent sketch from the
#                 budget. Pass 0 for a pass that builds no sketch.
#   bytes_per_nz  in-memory bytes per stored value; 12 for a chunk landing in
#                 a sparse matrix, ~48 for a collected triplet frame.
#   ram_frac      fraction of free RAM to budget per chunk. Defaults to
#                 `.ram_frac()`.
#
# The model is rows x payload-per-row, so nothing here is matrix-specific --
# a tabular store can use it by supplying its own shape and bytes_per_nz.
#
# Returns the recommended window (integer), invisibly.
.recommend_chunk_size <- function(
    n_cells,
    n_genes,
    density,
    k        = 50L,
    bytes_per_nz = 12,
    ram_frac = .ram_frac()
) {
    n_cells <- as.numeric(n_cells)
    n_genes <- as.numeric(n_genes)
    k       <- as.numeric(k)
    p       <- 10                               # fixed oversampling

    # Fallback when the machine's free RAM cannot be read: return the
    # configured default rather than a guessed budget, so a failed probe is
    # visibly a fallback instead of a silent 4 GB assumption.
    free_bytes <- tryCatch(.sc_free_ram(), error = function(e) NA_real_)
    if (!isTRUE(is.finite(free_bytes)) || free_bytes <= 0) {
        fallback <- as.integer(getOption("giottodisk.chunk_size", 250000L))
        return(as.integer(min(n_cells, fallback)))
    }

    as.integer(min(n_cells, max(.chunk_min_rows(),
        floor(.chunk_budget(free_bytes, ram_frac, n_cells, k, p) /
              .chunk_bytes_per_row(n_genes, density, bytes_per_nz)))))
}


# The two halves of the model, split out so the reporting verb can show the
# same arithmetic it is describing rather than reimplementing it.
#
# Per-row payload is rows x columns x fill x bytes-per-stored-value. Nothing
# matrix-specific about it -- `n_genes` is just the view's column count.
.chunk_bytes_per_row <- function(n_genes, density, bytes_per_nz) {
    as.numeric(n_genes) * as.numeric(density) * as.numeric(bytes_per_nz)
}

# Budget is a fraction of free RAM, less the persistent PCA sketch (which
# lives for the whole pass, not per chunk), floored at `.chunk_min_frac()` of
# free RAM (default 5%) so a large sketch cannot drive the budget to zero.
.chunk_budget <- function(free_bytes, ram_frac, n_cells, k = 50, p = 10) {
    sketch <- as.numeric(n_cells) * (as.numeric(k) + p) * 8
    max(free_bytes * ram_frac - sketch, free_bytes * .chunk_min_frac())
}


# Fraction of free RAM to budget
.ram_frac <- function() getOption("giottodisk.chunk_ram_frac", 0.25)

# Floors of the derivation. The window never drops below `chunk_min_rows`
# cells (per-window overhead dominates below it) and the budget never below
# `chunk_min_frac` of free RAM. Lower them on a machine where even the floor
# does not fit.
.chunk_min_rows <- function() {
    as.numeric(getOption("giottodisk.chunk_min_rows", 10000))
}
.chunk_min_frac <- function() {
    as.numeric(getOption("giottodisk.chunk_min_frac", 0.05))
}

# detect free RAM — macOS, Linux, Windows
.sc_free_ram <- function() {
    sysname <- Sys.info()[["sysname"]]
    if (sysname == "Darwin") {
        pg <- as.numeric(system("sysctl -n hw.pagesize", intern = TRUE))
        vm <- system("vm_stat", intern = TRUE)
        .val <- function(pat) {
            ln <- grep(pat, vm, value = TRUE, fixed = TRUE)[1L]
            if (is.na(ln)) return(0)
            as.numeric(regmatches(ln, regexpr("[0-9]+", ln))) * pg
        }
        .val("Pages free:") + .val("Pages inactive:") + .val("Pages speculative:")
    } else if (sysname == "Linux") {
        mem <- readLines("/proc/meminfo")
        ln  <- grep("^MemAvailable", mem, value = TRUE)[1L]
        as.numeric(regmatches(ln, regexpr("[0-9]+", ln))) * 1024
    } else if (.Platform$OS.type == "windows") {
        out <- system("wmic OS get FreePhysicalMemory /value",
                      intern = TRUE, ignore.stderr = TRUE)
        ln  <- grep("FreePhysicalMemory", out, value = TRUE)[1L]
        as.numeric(regmatches(ln, regexpr("[0-9]+", ln))) * 1024
    } else {
        stop("unsupported OS")
    }
}




# storeChunkInfo ####

#' @name storeChunkInfo
#' @title Report how a store's streaming windows are chosen
#' @description
#' Streaming passes read a store in windows sized to stay within a fraction of
#' free RAM. The window is **derived per read** rather than stored, from the
#' store's shape and its cached `@stats` marginals, so it adapts to the machine
#' doing the reading. This reports what those windows come out to and why.
#'
#' Two read shapes are shown because the package has two, and they differ by
#' roughly 4x in memory per stored value:
#'
#' \describe{
#'   \item{`sparse matrix`}{chunks land in a `dgCMatrix` (4-byte index +
#'     8-byte value). Used by the PCA band loops and the `storeWrite()` bake.}
#'   \item{`triplet frame`}{chunks are collected as `row_id` / `col_id` /
#'     `value` / `source_id` plus arrow and data.table overhead, measured at
#'     47-54 bytes per stored value. Used by R-side statistic accumulation
#'     when the op chain cannot be lowered to Acero.}
#' }
#'
#' @section Steering the window:
#' \describe{
#'   \item{`giottodisk.chunk_ram_frac`}{fraction of free RAM to budget per
#'     window (default 0.25), which also sizes the DuckDB memory limit. Lower
#'     it on a busy machine; raise it on a dedicated one. Scales every pass without knowing any store's shape. This
#'     is the one to reach for first.}
#'   \item{`giottodisk.chunk_size`}{pins an absolute window in cells, overriding
#'     the derivation entirely. An escape hatch — for a machine whose free
#'     memory cannot be read, for a shape the model sizes badly, and for tests
#'     that need to force several windows. Also the fallback value when the
#'     derivation cannot run.}
#'   \item{`giottodisk.chunk_min_rows`}{the smallest window the derivation
#'     returns, in cells (default 10,000). On a machine where a window of this
#'     size does not fit, lower it together with `giottodisk.chunk_ram_frac`.}
#'   \item{`giottodisk.chunk_min_frac`}{the smallest budget, as a fraction of
#'     free RAM (default 0.05), which applies when a large PCA sketch would
#'     otherwise leave no room for a window.}
#'   \item{`giottodisk.gef_batch_rows`}{pins the records per batch when a
#'     Stereo-seq `.gef` is written to a store. Derived from the same budget
#'     otherwise.}
#'   \item{`giottodisk.h5_batch_rows`}{pins the records per batch when a 10x
#'     `.h5` matrix is written to a store. Derived from the same budget
#'     otherwise.}
#' }
#'
#' Neither switches streaming on or off, and neither is a performance dial.
#' Windowing is not a mode: the loop always runs, and a budget covering the whole
#' view yields a single window. A smaller window lowers peak memory and adds
#' per-window overhead; a larger one does the reverse. The derived value is
#' already the largest window the budget allows, which is the fastest one that
#' stays bounded — so pinning `giottodisk.chunk_size` large does not make a pass
#' "unchunked", it removes the bound that keeps it from failing on a large store.
#'
#' @section Reproducibility:
#' The window count is derived from free RAM **at call time**, so it is not a
#' property of the store and can differ between runs on one machine.
#'
#' Combining window parts is exact for counts. For floating-point sums it
#' reassociates the addition, so a float statistic can differ by a few units in
#' the last place between runs with different window counts (measured 2-3 ULP).
#' Results are reproducible to tolerance, not bitwise: compare them with a
#' tolerance rather than hashing them or snapshotting exact digits.
#'
#' @section Measured instead of modelled:
#' The table is a model: it assumes a fixed number of bytes per stored value.
#' [storeChunkProbe()] reads one real window on this machine and reports the
#' time and memory it took, and the window size that fits a given budget.
#'
#' @section Further reading:
#' Which passes window and when, what forces a window, and why it has to be the
#' cell axis: `vignette("expression_windows", package = "GiottoDisk")`. The
#' decision, and the comparison against a spilling engine, are recorded in
#' `adr/0011-cell-windowing-over-spill.md`.
#'
#' @param x a `parquetExprStore` or `unionParquetExprStore`.
#' @param ram_frac numeric. Fractions to tabulate. Defaults to a spread around
#'   the configured value.
#' @param verbose logical. Print the report (default `TRUE`).
#' @param ... unused.
#' @returns A `data.table` of one row per (pass, ram_frac), invisibly.
#' @examples
#' \dontrun{
#' storeChunkInfo(pe)
#' options(giottodisk.chunk_ram_frac = 0.10)   # halve every streaming window
#' }
#' @export
setMethod("storeChunkInfo", signature("parquetExprBase"), function(
    x, ram_frac = c(0.10, 0.15, 0.20, 0.25, 0.35, 0.50),
    verbose = TRUE, ...
) {
    # NSE bindings
    pass <- chunk_rows <- n_chunks <- peak_mb <- current <- NULL

    n_cells <- as.numeric(x@n_cells)
    n_genes <- as.numeric(x@n_genes)
    if (n_cells <= 0 || n_genes <= 0) {
        stop("[storeChunkInfo] store has no cells or no features.",
             call. = FALSE)
    }
    dens <- .pe_view_density(x)
    if (!isTRUE(is.finite(dens)) || dens <= 0) {
        stop("[storeChunkInfo] could not determine density for this store.",
             call. = FALSE)
    }

    # Same two shapes the real consumers ask for; see .exprbase_chunk_size()
    # and .pe_accum_chunk_size().
    shapes <- list(
        list(pass = "sparse matrix", bytes_per_nz = 12, k = 50L),
        list(pass = "triplet frame", bytes_per_nz = 48, k = 0L)
    )
    cfg_frac <- .ram_frac()
    fracs    <- sort(unique(c(as.numeric(ram_frac), cfg_frac)))

    out <- data.table::rbindlist(lapply(shapes, function(sh) {
        bpr <- .chunk_bytes_per_row(n_genes, dens, sh$bytes_per_nz)
        rows <- vapply(fracs, function(f) {
            .recommend_chunk_size(n_cells, n_genes, dens, k = sh$k,
                bytes_per_nz = sh$bytes_per_nz, ram_frac = f)
        }, integer(1L))
        data.table::data.table(
            pass        = sh$pass,
            ram_frac    = fracs,
            chunk_rows  = rows,
            n_chunks    = as.integer(ceiling(n_cells / rows)),
            peak_mb     = as.integer(round(rows * bpr / 1e6)),
            current     = ifelse(fracs == cfg_frac, "<--", "")
        )
    }))

    if (isTRUE(verbose)) .print_chunk_info(x, out, dens, n_cells, n_genes)
    invisible(out)
})


.print_chunk_info <- function(x, tbl, dens, n_cells, n_genes) {
    fmt <- function(v) format(v, big.mark = ",", scientific = FALSE)
    free <- tryCatch(.sc_free_ram(), error = function(e) NA_real_)
    nnz  <- if (inherits(x, "unionParquetExprStore")) {
        sum(vapply(x@stores, .pestore_view_nnz, numeric(1L)))
    } else .pestore_view_nnz(x)
    pinned <- getOption("giottodisk.chunk_size")

    cat(sprintf("── storeChunkInfo %s\n", strrep("─", 51)))
    cat(sprintf("   Store    : %s\n", class(x)[1L]))
    cat(sprintf("   View     : %s features x %s cells (density %.3f)\n",
        fmt(n_genes), fmt(n_cells), dens))
    cat(sprintf("   Marginals: %s\n",
        if (isTRUE(is.finite(nnz)))
            sprintf("cached, %s stored values", fmt(round(nnz)))
        else "not cached — density was counted"))
    cat(sprintf("   Free RAM : %s\n",
        if (isTRUE(is.finite(free))) sprintf("%.1f GB", free / 1e9)
        else "undetectable — windows fall back to giottodisk.chunk_size"))
    if (!is.null(pinned)) {
        cat(sprintf("   PINNED   : giottodisk.chunk_size = %s — every window\n",
            fmt(as.integer(pinned))))
        cat("              below is overridden by this value.\n")
    }
    cat("\n")
    print(tbl, row.names = FALSE)
    cat(sprintf(
        "\n   ram_frac marked <-- is the configured value%s.\n\n",
        if (is.null(getOption("giottodisk.chunk_ram_frac"))) " (default)" else ""
    ))
    invisible(NULL)
}


# storeChunkProbe ####

#' @name storeChunkProbe
#' @title Measure a streaming window on this machine
#' @description
#' [storeChunkInfo()] models how large a streaming window is allowed to be.
#' `storeChunkProbe()` measures it instead: it reads one window of `cells`
#' cells into a sparse matrix, the read most streaming passes repeat window by
#' window (enrichment, deconvolution, the PCA band loops), and reports how long
#' it took and how much memory it needed. From the measured memory per stored
#' value it then works out the window that fits the budget, and how long a
#' full pass over the store takes.
#'
#' Run it once on a new machine, or when a streaming step needs more memory
#' than you can spare. Pin the suggested window with
#' `options(giottodisk.chunk_size = )`, or scale every window at once with
#' `giottodisk.chunk_ram_frac`.
#'
#' Memory is measured as the peak R heap above what was in use before the
#' read, plus the peak growth of arrow's memory pool. Memory the operating
#' system attributes to the process can be somewhat higher.
#'
#' @param x a `parquetExprStore` or `unionParquetExprStore`.
#' @param cells integer. Cells in the probed window (default 10,000, or all
#'   cells when there are fewer).
#' @param ram_frac numeric. Fraction of free RAM to budget per window. Defaults
#'   to the configured `giottodisk.chunk_ram_frac`.
#' @param verbose logical. Print the report (default `TRUE`).
#' @param ... unused.
#' @returns A one-row `data.table`, invisibly: the probed window (`cells`,
#'   `stored_values`, `seconds`, `peak_mb`, `bytes_per_value`) and the window
#'   that fits the budget (`free_mb`, `budget_mb`, `suggested_cells`,
#'   `n_windows`, `pass_seconds`).
#' @seealso [storeChunkInfo()] for the modelled windows and the options that
#'   steer them.
#' @examples
#' \dontrun{
#' pe <- getExpression(g, values = "normalized", output = "exprObj")[]
#' probe <- storeChunkProbe(pe)
#' options(giottodisk.chunk_size = probe$suggested_cells)
#' }
#' @export
setMethod("storeChunkProbe", signature("parquetExprBase"), function(
    x, cells = 10000L, ram_frac = .ram_frac(), verbose = TRUE, ...
) {
    n_cells <- as.numeric(x@n_cells)
    n_genes <- as.numeric(x@n_genes)
    if (n_cells <= 0 || n_genes <= 0) {
        stop("[storeChunkProbe] store has no cells or no features.",
             call. = FALSE)
    }
    checkmate::assert_number(cells, lower = 1)
    checkmate::assert_number(ram_frac, lower = 0, upper = 1)
    cells <- as.integer(min(cells, n_cells))
    window <- .pe_window_store(.pe_windows(x, cells)[[1L]])
    # warm-up: opening the dataset and the first scan are one-off costs that
    # every later window skips, so keep them out of the measured read
    invisible(storeRead(window[, seq_len(min(100L, cells))],
        output = "dgcmatrix", max_rows = Inf, max_cols = Inf))

    mb <- function(b) as.integer(round(b / 1e6))
    pool <- arrow::default_memory_pool()
    pool_before <- pool$max_memory
    heap_before <- .heap_bytes(gc(reset = TRUE), "used")
    t0 <- Sys.time()
    m <- storeRead(window, output = "dgcmatrix", max_rows = Inf,
        max_cols = Inf)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    heap_peak <- .heap_bytes(gc(), "max used") - heap_before
    pool_peak <- max(0, pool$max_memory - pool_before)
    nnz <- length(m@x)
    rm(m)

    peak <- heap_peak + pool_peak
    bpv <- if (nnz > 0) peak / nnz else NA_real_
    dens <- .pe_view_density(x)
    free <- tryCatch(.sc_free_ram(), error = function(e) NA_real_)
    budget <- if (isTRUE(is.finite(free))) {
        .chunk_budget(free, ram_frac, n_cells, k = 0L)
    } else NA_real_
    per_row <- .chunk_bytes_per_row(n_genes, dens, bpv)
    suggested <- if (isTRUE(is.finite(budget)) && isTRUE(per_row > 0)) {
        as.integer(min(n_cells, max(1, floor(budget / per_row))))
    } else NA_integer_
    n_win <- as.integer(ceiling(n_cells / suggested))

    out <- data.table::data.table(
        cells           = cells,
        stored_values   = nnz,
        seconds         = round(secs, 2),
        peak_mb         = mb(peak),
        bytes_per_value = round(bpv, 1),
        free_mb         = mb(free),
        budget_mb       = mb(budget),
        suggested_cells = suggested,
        n_windows       = n_win,
        pass_seconds    = round(secs / cells * n_cells, 1)
    )
    if (isTRUE(verbose)) .print_chunk_probe(x, out, ram_frac)
    invisible(out)
})

# R heap in bytes from a gc() table: Ncells are 56 bytes on 64-bit R, Vcells 8
.heap_bytes <- function(g, col) {
    sum(g[, col] * c(56, 8))
}

.print_chunk_probe <- function(x, out, ram_frac) {
    fmt <- function(v) format(v, big.mark = ",", scientific = FALSE)
    o <- as.list(out)
    cat(sprintf("── storeChunkProbe %s\n", strrep("─", 50)))
    cat(sprintf("   Store    : %s, %s features x %s cells\n", class(x)[1L],
        fmt(x@n_genes), fmt(x@n_cells)))
    cat(sprintf("   Probe    : %s cells, %s stored values\n",
        fmt(o$cells), fmt(o$stored_values)))
    cat(sprintf("   Measured : %.1f s, %s MB peak (%.0f bytes per value)\n",
        o$seconds, fmt(o$peak_mb), o$bytes_per_value))
    if (is.na(o$suggested_cells)) {
        cat("   Free RAM : undetectable; no window suggested\n\n")
        return(invisible(NULL))
    }
    cat(sprintf("   Budget   : %s MB (%s of %s MB free)\n",
        fmt(o$budget_mb), format(ram_frac), fmt(o$free_mb)))
    cat(sprintf("   Window   : %s cells fit the budget -> %s windows,\n",
        fmt(o$suggested_cells), fmt(o$n_windows)))
    cat(sprintf("              about %s s to read the whole store once\n",
        fmt(o$pass_seconds)))
    cat(sprintf("\n   options(giottodisk.chunk_size = %s) pins this window.\n\n",
        o$suggested_cells))
    invisible(NULL)
}
