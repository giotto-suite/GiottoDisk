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
    nnz <- storeRead(pe) |>
        dplyr::count() |>
        dplyr::collect() |>
        dplyr::pull(n)
    as.numeric(nnz) / (n_cells * n_genes)
}

# Finalizer: auto-tune store@chunk_size from dataset shape + density,
# using sc_recommend_chunk's defaults (n_hvg = min(2000L, n_genes),
# k = 50L, ram_frac = 0.25). Restores the auto-sizing legacy did via
# createGiottoFromParquet, but applied at storeWrite time so every
# Input pathway and the direct memoryMatrix path benefit.
.pestore_finalize_chunk_size <- function(pe, verbose = FALSE) {
    if (pe@n_cells <= 0 || pe@n_genes <= 0) return(pe)
    pe@chunk_size <- sc_recommend_chunk(
        n_cells = pe@n_cells,
        n_genes = pe@n_genes,
        density = .pestore_density(pe),
        verbose = verbose
    )
    pe
}


# Build chunk boundaries that respect duplicate-name groups: never split
# a run of raw geneDT rows that share the same name_to_row between two
# chunks. Returns a list of c(g_lo, g_hi) integer pairs covering 1..n.
# `target_size` is the desired chunk size in raw geneDT rows.
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
