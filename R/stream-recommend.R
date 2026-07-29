#' @include class-parquetExprStore.R
NULL

# stream-recommend ####
# Estimates an optimal streaming chunk size given the dataset shape, planned
# pipeline parameters, and currently free RAM.

#' @name sc_recommend_chunk
#' @title Recommend a streaming chunk size from free RAM
#' @description
#' Estimates the optimal `chunk_size` for streaming reads of a
#' [parquetExprStore-class] so that each pipeline pass stays within a
#' target fraction of currently free RAM. Call this once before opening
#' a dataset and pass the result to the store constructor.
#'
#' Prints a small comparison table across RAM budgets and returns the
#' recommended value invisibly.
#'
#' @param n_cells numeric. Total cell count in the dataset.
#' @param n_genes numeric. Total gene count.
#' @param density numeric. Approximate fill fraction of the matrix.
#' @param n_hvg integer. Planned number of HVG (default 2000).
#' @param k integer. Planned number of PCA dimensions (default 50).
#' @param ram_frac numeric. Fraction of free RAM to budget per chunk
#'   (default 0.25). Lower (0.15) on machines with heavy background
#'   load; higher (0.40) on a dedicated analysis server.
#' @param verbose logical. Print the comparison table and recommendation
#'   line. Default `TRUE`.
#' @return Recommended chunk size (integer), returned invisibly.
#' @examples
#' chunk <- sc_recommend_chunk(2e6, 5000, density = 0.30,
#'                              n_hvg = 1000, k = 50)
#' @export
sc_recommend_chunk <- function(
    n_cells,
    n_genes,
    density,
    n_hvg    = min(2000L, as.integer(n_genes)),
    k        = 50L,
    ram_frac = 0.25,
    verbose  = TRUE
) {
    n_cells <- as.numeric(n_cells)
    n_genes <- as.numeric(n_genes)
    n_hvg   <- as.numeric(n_hvg)
    k       <- as.numeric(k)
    p       <- 10                               # fixed oversampling

    free_bytes <- tryCatch(
        .sc_free_ram(),
        error = function(e) {
            if (isTRUE(verbose)) {
                message("[sc_recommend_chunk] RAM detection failed — ",
                        "assuming 4 GB free.")
            }
            4e9
        }
    )

    bytes_qc       <- n_genes * density * 12
    bytes_pca      <- n_hvg   * density * 12 + (k + p) * 8
    bytes_per_cell <- max(bytes_qc, bytes_pca)

    sketch_bytes <- n_cells * (k + p) * 8
    sketch_mb    <- round(sketch_bytes / 1e6)

    budget    <- free_bytes * ram_frac - sketch_bytes
    budget    <- max(budget, free_bytes * 0.05)
    chunk_rec <- as.integer(
        min(n_cells, max(10000, floor(budget / bytes_per_cell)))
    )

    if (isTRUE(verbose)) {
        fracs <- c(0.10, 0.15, 0.20, 0.25, 0.35, 0.50)
        # NSE bindings for data.table
        chunk_size <- recommended <- NULL
        tbl <- data.table::data.table(
            `RAM budget` = paste0(as.integer(fracs * 100), "%"),
            chunk_size   = vapply(fracs, function(f) {
                b <- free_bytes * f - sketch_bytes
                b <- max(b, free_bytes * 0.05)
                as.integer(min(n_cells, max(10000, floor(b / bytes_per_cell))))
            }, integer(1L))
        )
        tbl[, `n chunks/pass`   := as.integer(ceiling(n_cells / chunk_size))]
        tbl[, `chunk peak (MB)` := as.integer(round(chunk_size * bytes_per_cell / 1e6))]
        tbl[, recommended       := ifelse(chunk_size == chunk_rec, "<--", "")]

        cat(sprintf(
            "── sc_recommend_chunk %s\n",
            paste(rep("─", 47), collapse = "")
        ))
        fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)
        cat(sprintf("   Dataset  : %s cells × %s genes  (density %.2f)\n",
            fmt(n_cells), fmt(n_genes), density))
        cat(sprintf("   Pipeline : n_hvg = %s  |  k = %d\n",
            fmt(n_hvg), as.integer(k)))
        cat(sprintf("   Free RAM : %.1f GB detected\n", free_bytes / 1e9))
        cat(sprintf("   Sketch Y : %s MB  (persistent during PCA)\n\n",
            format(sketch_mb, big.mark = ",")))
        print(tbl, row.names = FALSE)
        cat(sprintf(
            "\n   ► Use chunk_size = %s  (%s MB peak per chunk)\n\n",
            format(chunk_rec, big.mark = ","),
            format(as.integer(round(chunk_rec * bytes_per_cell / 1e6)),
                   big.mark = ",")
        ))
    }
    invisible(chunk_rec)
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
