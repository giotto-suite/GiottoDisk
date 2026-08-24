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
# Two options steer it:
#   giottodisk.chunk_ram_frac   fraction of free RAM to budget (default 0.25)
#   giottodisk.chunk_size       pin an absolute window; also the fallback
#                               when the derivation cannot run

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
#                 getOption("giottodisk.chunk_ram_frac", 0.25).
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
    ram_frac = getOption("giottodisk.chunk_ram_frac", 0.25)
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

    as.integer(min(n_cells, max(10000,
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
# lives for the whole pass, not per chunk), floored at 5% so a large sketch
# cannot drive the budget to zero.
.chunk_budget <- function(free_bytes, ram_frac, n_cells, k = 50, p = 10) {
    sketch <- as.numeric(n_cells) * (as.numeric(k) + p) * 8
    max(free_bytes * ram_frac - sketch, free_bytes * 0.05)
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


# pestore-chunking ####
#
# The one user-facing page on windowing. Everything that streams links here
# rather than restating a knob list that would drift. Prose lives in @sections
# so a topic that needs one inline can `@inheritSection` it instead of copying.

#' @name pestore-chunking
#' @aliases chunking parquetExprStore-chunking
#' @title How expression stores are streamed, and which knobs exist
#' @description
#' Passes over a `parquetExprStore` that cannot hold the whole store in memory
#' read it in **cell windows**. This page is the reference for how a window is
#' chosen, what the two options do, and what to expect of results. Individual
#' verbs link here instead of restating it.
#'
#' [storeChunkInfo()] reports the numbers for a given store on the current
#' machine.
#'
#' @section What a window is:
#' A window is a contiguous range of cells. A pass visits the windows in cell
#' order, does its work on each, and combines the parts; for a sum, a count or a
#' sum of squares that combination is addition, so the answer is the same one an
#' unwindowed pass would give.
#'
#' The axis is not a free choice. Stores are written cell-major, so a contiguous
#' cell range narrows to a range predicate that lets the parquet reader skip row
#' groups — each window reads only its own slice, and N windows still cost one
#' pass over the data. Narrowing by **feature** cannot do this: features are not
#' the sort key, so a feature batch reads the whole store and the cost grows with
#' the number of batches. If you are tempted to loop over gene blocks to stay
#' inside memory, ask for all the features at once instead and let the window do
#' the bounding.
#'
#' Statistics that need a global ordering along the feature axis — a rank, a
#' median, a Wilcoxon test — cannot be combined from windows at all, and have no
#' streaming path today.
#'
#' @section Not a mode:
#' There is no chunked mode and no unchunked mode, and neither option switches
#' one on. The window loop always runs; when the budget covers the whole view it
#' yields a single window, which is a single plan over the store. The options
#' change the *bound*, not the code path.
#'
#' In particular, pinning `giottodisk.chunk_size` to a large value does not make
#' a pass "unchunked" — it removes the memory bound that keeps the pass from
#' failing on a large store. Reach for it when you know the data fits.
#'
#' Which passes window at all is narrower than it looks. Grouped statistics
#' (marker detection, anything passing `groups`) window because the grouping
#' puts a join in front of an otherwise small aggregate. Ungrouped statistics on
#' a lowerable chain run as one streamed plan and do not window. PCA windows, but
#' accumulates in its own order.
#'
#' @section Steering the window:
#' \describe{
#'   \item{`giottodisk.chunk_ram_frac`}{fraction of free RAM to budget per
#'     window (default 0.25). Lower it on a busy machine; raise it on a
#'     dedicated one. Scales every pass without knowing any store's shape. This
#'     is the one to reach for first.}
#'   \item{`giottodisk.chunk_size`}{pins an absolute window in cells, overriding
#'     the derivation entirely. An escape hatch — for a machine whose free
#'     memory cannot be read, for a shape the model sizes badly, and for tests
#'     that need to force several windows. Also the fallback value when the
#'     derivation cannot run.}
#' }
#'
#' Neither is a performance dial. A smaller window lowers peak memory and adds
#' per-window overhead; a larger one does the reverse. The derived value is
#' already the largest window the budget allows, which is the fastest one that
#' stays bounded.
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
#' This is invisible to analysis — it is 15 significant figures in — with one
#' thing to keep in view: a **discrete** choice taken on a float statistic, such
#' as a top-N feature selection, can in principle land differently if two
#' candidates are separated by less than that. Ties that close are not
#' meaningful distinctions, but they can move a downstream result that is keyed
#' to the selection.
#'
#' @section Why the window is not stored:
#' It depends on free RAM, which is a property of the machine doing the reading,
#' not of the data. Baking it in at write time would make it stale as soon as the
#' store moved. What *is* machine-independent is the payload size per stored
#' value, and that is what the `@stats` marginals cache.
#'
#' A store with no cached marginals — one built by [parquetExprStore()] straight
#' from a path, where `scan_stats = FALSE` is the default — has to count its
#' nonzeros before it can size a window, which costs one extra pass. Stores from
#' `storeWrite()` and the convenience importers carry marginals already, and they
#' survive `snapshotSave()` / `snapshotLoad()`.
#'
#' @seealso [storeChunkInfo()] for the numbers on a given store.
NULL


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
#' @inheritSection pestore-chunking Steering the window
#' @inheritSection pestore-chunking Reproducibility
#'
#' @seealso [pestore-chunking] for what a window is, why the cell axis, and why
#'   windowing is not a mode you opt into.
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
    cfg_frac <- getOption("giottodisk.chunk_ram_frac", 0.25)
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
