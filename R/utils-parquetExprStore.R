#' @include class-parquetExprStore.R
NULL


# Disambiguate duplicate feature IDs the 10x way (suffix `--N`).
# Matches the convention used by Giotto::get10Xmatrix so the
# parquet path produces the same feat_ids as the in-memory dgCMatrix
# path. No-op when all names are already unique.
.disambiguate_feat_ids <- function(feat_ids) {
    feat_ids <- as.character(feat_ids)
    counts <- table(feat_ids)
    dups   <- names(counts)[counts > 1L]
    if (length(dups) == 0L) return(feat_ids)
    out <- feat_ids
    for (nm in dups) {
        idx <- which(feat_ids == nm)
        out[idx] <- paste0(nm, "--", seq_along(idx))
    }
    out
}


# Approximate density (nnz / (n_cells * n_genes)) of a parquetExprStore.
# Uses Arrow's count() — for an unfiltered store this resolves via
# parquet footer metadata (num_rows per row-group), no data scan.
# Subsetted stores pay a predicate-pushdown filter pass.
.pestore_density <- function(pe) {
    n_cells <- as.numeric(pe@n_cells)
    n_genes <- as.numeric(pe@n_genes)
    if (n_cells <= 0 || n_genes <= 0) return(0)
    # Cached marginals answer this without touching the data.
    cached <- .pestore_view_nnz(pe)
    if (isTRUE(is.finite(cached))) return(cached / (n_cells * n_genes))
    nnz <- storeRead(pe) |>
        dplyr::count() |>
        dplyr::collect() |>
        dplyr::pull(n)
    as.numeric(nnz) / (n_cells * n_genes)
}

# Streaming window, in cells, for a pass that materializes `bytes_per_nz`
# bytes per stored value.
#
# Derived per call rather than stored on the object. The window depends on
# free RAM, which is a property of the machine doing the reading -- baking it
# into the store at write time made it stale the moment the store moved. What
# IS machine-independent is the payload size, and that is what `@stats` caches.
#
# Two options steer this. `giottodisk.chunk_size` serves double duty: set
# explicitly it pins the window (an escape hatch for constrained environments
# and for tests), and unset its default is the value used when the derivation
# cannot run -- no cached marginals, or a machine whose free RAM cannot be
# read. `giottodisk.chunk_ram_frac` instead scales the budget the derivation
# works from, which is the knob to reach for when the shape is fine but the
# machine is busier than the default assumes.
#
# Union stores carry no marginals of their own; sum the substores'.
.pe_window_cells <- function(pe, bytes_per_nz = 12, k = 50L) {
    pinned <- getOption("giottodisk.chunk_size")
    n_cells <- as.numeric(pe@n_cells)
    n_genes <- as.numeric(pe@n_genes)
    if (!is.null(pinned)) {
        return(as.integer(max(1, min(n_cells, as.numeric(pinned)))))
    }
    # fallback value for unknown cases
    if (n_cells <= 0 || n_genes <= 0) return(250000L)

    dens <- tryCatch(.pe_view_density(pe), error = function(e) NA_real_)
    if (!isTRUE(is.finite(dens)) || dens <= 0) {
        return(as.integer(min(n_cells, 250000L)))
    }

    as.integer(.recommend_chunk_size(
        n_cells      = n_cells,
        n_genes      = n_genes,
        density      = dens,
        k            = k,
        bytes_per_nz = bytes_per_nz
    ))
}


# Fill fraction of the current view, for either store kind.
#
# `.pestore_density()` reads the cached marginals when they exist and counts
# when they do not -- a store handed a path directly (`storeCreate()`, or the
# constructor) never went through `storeWrite()` and has none, and one scan
# beats sizing a window blind.
.pe_view_density <- function(pe) {
    n_cells <- as.numeric(pe@n_cells)
    n_genes <- as.numeric(pe@n_genes)
    if (n_cells <= 0 || n_genes <= 0) return(0)
    if (inherits(pe, "unionParquetExprStore")) {
        return(sum(vapply(pe@stores, .pestore_view_nnz, numeric(1L))) /
               (n_cells * n_genes))
    }
    .pestore_density(pe)
}


# Finalizer: cache the Parquet's marginal nonzero counts on @stats.
#
# Two grouped counts over the freshly written file, one per axis. Replaces
# the old chunk_size finalizer, which spent a density pass at write time to
# bake a RAM-derived window into the store -- a value that goes stale the
# moment the store is opened on a different machine. Counting instead gives
# a machine-independent fact, and every later consumer derives its own window
# from it without touching the data again.
#
# Called on the store as written, so no subset or op chain is in play and the
# ids are the file's own. `[` never invalidates the result.
.pestore_finalize_stats <- function(pe, verbose = FALSE) {
    if (pe@n_cells <= 0 || pe@n_genes <= 0) return(pe)
    row_id <- col_id <- NULL   # NSE bindings

    .marginal <- function(key, n_out) {
        agg <- storeRead(pe, output = "query") |>
            dplyr::count(!!rlang::sym(key)) |>
            dplyr::collect()
        out <- integer(n_out)
        idx <- as.integer(agg[[key]])
        keep <- !is.na(idx) & idx >= 1L & idx <= n_out
        out[idx[keep]] <- as.integer(agg[["n"]][keep])
        out
    }

    pe@stats <- list(
        col_nnz = .marginal("col_id", as.integer(pe@n_genes)),
        row_nnz = .marginal("row_id", as.integer(pe@n_cells))
    )
    if (isTRUE(verbose)) {
        message("[storeWrite] cached marginals: ",
                format(sum(pe@stats$col_nnz), big.mark = ","), " nonzeros")
    }
    pe
}


# Total stored values in the CURRENT VIEW, from the cached marginals.
#
# Exact when at most one axis is narrowed. With both narrowed, the gene axis
# stays exact and the cell axis contributes its kept fraction -- nonzeros are
# far from uniform across features (HVG selection picks the densest rows on
# purpose) but reasonably uniform across cells, so the exact axis is the one
# that matters. Returns NA when no marginals are cached.
.pestore_view_nnz <- function(pe) {
    st <- pe@stats
    if (!length(st) || is.null(st$col_nnz) || is.null(st$row_nnz)) {
        return(NA_real_)
    }
    gi <- pe@gene_idx
    ci <- pe@cell_idx
    n_cells_file <- length(st$row_nnz)

    # Neither axis narrowed, or exactly one: read the answer off the marginal
    # for that axis. Only a two-axis subset needs the uniformity assumption,
    # and there the gene axis stays exact -- it is the one where nonzeros are
    # genuinely skewed, since feature selection targets the densest rows.
    if (length(gi) == 0L && length(ci) == 0L) {
        return(sum(as.numeric(st$col_nnz), na.rm = TRUE))
    }
    if (length(ci) == 0L) {
        return(sum(as.numeric(st$col_nnz[gi]), na.rm = TRUE))
    }
    if (length(gi) == 0L) {
        return(sum(as.numeric(st$row_nnz[ci]), na.rm = TRUE))
    }
    sum(as.numeric(st$col_nnz[gi]), na.rm = TRUE) *
        (length(ci) / max(n_cells_file, 1L))
}


# ---- op payloads -----------------------------------------------------------
#
# Same indexing convention as the `@stats` marginals above: a numeric vector
# whose POSITION is the on-disk row_id / col_id. Kept together because that
# convention is the thing to preserve -- it is what makes both invariant under
# `[`, and what lets the R-side executor index directly where arrow has to
# join.

# Resolve a payload to a full-length numeric vector for one substore.
.pe_axis_payload_vec <- function(payload, uid, n) {
    if (is.null(payload)) return(NULL)
    if (!is.list(payload)) return(rep_len(as.numeric(payload), n))
    v <- payload[[as.character(uid)]]
    if (is.null(v)) return(NULL)
    as.numeric(v)
}

# Build the joinable (source_id, <key>, w) table an Acero plan needs. Arrow
# cannot index an R vector from inside a query, so per-axis state has to arrive
# as a table -- rebuilt per call, which is what the previous executor did too.
.pe_axis_payload_table <- function(payload, key) {
    src <- names(payload)
    data.table::rbindlist(lapply(src, function(u) {
        v <- as.numeric(payload[[u]])
        data.table::data.table(
            source_id = rep_len(as.character(u), length(v)),
            key_id    = seq_along(v),
            w         = v
        )
    }))[!is.na(w)]
}


# ---- GEF ingest: cell-major batches ----------------------------------------
#
# Every expression store is written cell-major: each part file covers one
# ascending range of cells (AGENTS.md, "window the CELL axis"). A GEF stores
# its expression gene-major, so the readers either read a cell-major copy
# (cellbin `cellExp`) or reorder out of core through a spill (bin GEF, and
# cellbin files without `cellExp`). Duplicate gene names are summed by the
# final `(row_id, col_id)` aggregate: a cell never spans two batches, so no
# record for it can be left behind in another one.

# TRUE when the HDF5 file has an object at `name` (every level checked, since
# H5Lexists on a path through a missing group errors).
.h5_has <- function(path, name) {
    fid <- rhdf5::H5Fopen(path, flags = "H5F_ACC_RDONLY")
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    cur <- character(0L)
    for (p in strsplit(name, "/", fixed = TRUE)[[1L]]) {
        cur <- paste(c(cur, p), collapse = "/")
        if (!rhdf5::H5Lexists(fid, cur)) return(FALSE)
    }
    TRUE
}

# Records per ingest batch, from the same free-RAM model the streaming passes
# use. A GEF batch holds the compound read, the triplet frame and the
# aggregate's copy at once, about 96 B per record -- twice the collected
# triplet frame alone. `options(giottodisk.gef_batch_rows =)` pins it; when free
# RAM cannot be read, `.gef_batch_rows_fallback` records are used.
.gef_batch_rows_fallback <- 5e6

.gef_batch_rows <- function() {
    pinned <- getOption("giottodisk.gef_batch_rows")
    if (!is.null(pinned)) return(max(1, as.numeric(pinned)))
    free <- tryCatch(.sc_free_ram(), error = function(e) NA_real_)
    if (!isTRUE(is.finite(free)) || free <= 0) return(.gef_batch_rows_fallback)
    ram_frac <- getOption("giottodisk.chunk_ram_frac", 0.25)
    budget <- .chunk_budget(free, ram_frac, n_cells = 0, k = 0, p = 0)
    max(1e6, floor(budget / 96))
}

# Contiguous cell ranges of at most `max_rows` records each, cut greedily; a
# single cell with more records than that gets a range of its own. Returns a
# list of c(a, b).
.gef_cell_ranges <- function(n_per_cell, max_rows) {
    n <- length(n_per_cell)
    if (!n) return(list())
    cs <- c(0, cumsum(as.numeric(n_per_cell)))   # cs[i + 1] = records through cell i
    out <- list()
    a <- 1L
    while (a <= n) {
        b <- min(n, max(a, findInterval(cs[a] + max_rows, cs) - 1L))
        out[[length(out) + 1L]] <- c(a, b)
        a <- b + 1L
    }
    out
}

# Fine buckets for a spill of `n_records`: several per budget-sized batch, so
# that pass 2 can group consecutive buckets under the budget even when records
# are unevenly spread over the cell axis.
.gef_n_buckets <- function(n_records, max_rows, per_batch = 8L) {
    max(1L, as.integer(per_batch * ceiling(n_records / max_rows)))
}

# Two-pass iterator for gene-major records. Pass 1 (on the first
# `next_batch()`) calls `produce(i)` for i in seq_len(n_chunks); each returns a
# data.table with an integer bucket column `b` in 1..n_buckets and whatever
# `finish()` needs, and is spilled to one parquet file per (bucket, chunk)
# under the artifact dump directory. Pass 2 groups consecutive buckets while
# their records fit `max_rows` (a bucket larger than that goes alone) and
# returns `finish(dt)` for each group in ascending order: the batch's
# `(row_id, col_id, value)`, keyed. Buckets must partition the cell axis in
# ascending order for the output to be cell-major. The spill is removed on
# close.
.gef_spill_iterator <- function(produce, n_chunks, n_buckets, finish, max_rows) {
    b <- NULL  # data.table var
    spill <- .dump_tempfile()
    spilled <- FALSE
    groups <- list()
    g_i <- 0L
    closed <- FALSE

    run_pass1 <- function() {
        dir.create(spill, recursive = TRUE, showWarnings = FALSE)
        rows <- numeric(n_buckets)
        for (i in seq_len(n_chunks)) {
            dt <- produce(i)
            if (is.null(dt) || !nrow(dt)) next
            rows <- rows + tabulate(dt$b, nbins = n_buckets)
            parts <- split(dt, by = "b", keep.by = FALSE, sorted = TRUE)
            for (nm in names(parts)) {
                d <- file.path(spill, sprintf("b%06d", as.integer(nm)))
                dir.create(d, showWarnings = FALSE)
                .write_parquet_file(parts[[nm]], file.path(d, sprintf("c%06d.parquet", i)))
            }
        }
        # consecutive non-empty buckets, grouped under the budget
        nz <- which(rows > 0)
        cur <- integer(0L); acc <- 0
        for (k in nz) {
            if (length(cur) && acc + rows[k] > max_rows) {
                groups[[length(groups) + 1L]] <<- cur
                cur <- integer(0L); acc <- 0
            }
            cur <- c(cur, k); acc <- acc + rows[k]
        }
        if (length(cur)) groups[[length(groups) + 1L]] <<- cur
        spilled <<- TRUE
    }

    close_fn <- function() {
        closed <<- TRUE
        unlink(spill, recursive = TRUE)
        invisible(NULL)
    }

    next_batch <- function() {
        if (closed) return(NULL)
        if (!spilled) run_pass1()
        if (g_i >= length(groups)) return(NULL)
        g_i <<- g_i + 1L
        dirs <- file.path(spill, sprintf("b%06d", groups[[g_i]]))
        fs <- list.files(dirs, pattern = "[.]parquet$", full.names = TRUE)
        dt <- data.table::setDT(dplyr::collect(arrow::open_dataset(fs)))
        unlink(dirs, recursive = TRUE)
        finish(dt)
    }

    list(next_batch = next_batch, close = close_fn)
}


# Build chunk boundaries that respect duplicate-name groups: never split
# a run of raw geneDT rows that share the same name_to_row between two
# chunks. Non-adjacent duplicates are summed later, by the batch aggregate of
# `.gef_spill_iterator()`. Returns a list of c(g_lo, g_hi) integer pairs
# covering 1..n. `target_size` is the desired chunk size in raw geneDT rows.
.gef_safe_chunks <- function(name_to_row, target_size) {
    n <- length(name_to_row)
    if (n == 0L) return(list())
    target_size <- as.integer(target_size)

    # Boundary positions (1-based): position i is a "safe break" if
    # name_to_row[i] differs from name_to_row[i-1] (meaning the previous
    # group ends at i-1 and a new one starts at i). Position 1 is always
    # a starting point.
    nm <- name_to_row
    is_break <- c(TRUE, nm[-1] != nm[-n] |
                  (is.na(nm[-1]) != is.na(nm[-n])))
    is_break[is.na(is_break)] <- TRUE
    safe_starts <- which(is_break)

    # Walk safe_starts and group them into chunks of approximately
    # target_size raw genes each.
    out <- vector("list", length(safe_starts))
    k <- 0L
    last_start <- safe_starts[1L]
    for (s in safe_starts[-1L]) {
        if (s - last_start >= target_size) {
            k <- k + 1L
            out[[k]] <- c(last_start, s - 1L)
            last_start <- s
        }
    }
    k <- k + 1L
    out[[k]] <- c(last_start, n)
    out[seq_len(k)]
}
