#!/usr/bin/env Rscript
# Driver: run bench/_runner.R against two trees and report the difference.
#
#   Rscript bench/regress.R                              # vs upstream/dev, synthetic
#   Rscript bench/regress.R --ref=HEAD~1 --reps=5
#   Rscript bench/regress.R --workers=1                  # serial, for comparison
#   Rscript bench/regress.R --report                     # re-print last run
#   GD_BENCH_H5=/path/cell_feature_matrix.h5 Rscript bench/regress.R
#
# Not a test: nothing asserts a correct answer. It reports time, peak memory and
# the ratio between two trees, and exits non-zero if anything crossed the
# threshold. Outside the built package (`^bench$` in .Rbuildignore); results are
# gitignored.
#
# Ratios, not baselines. Absolute times are not comparable across machines,
# arrow versions or background load, so a recorded-baselines file goes stale
# immediately. What is stable is the ratio between two trees measured back to
# back on one machine -- hence the ref is an argument and nothing is committed
# but the scripts.
#
# Results stream to bench/results/ as each case finishes, so an interrupted run
# keeps what it measured and `--report` can still read it.

suppressWarnings(suppressMessages(library(data.table)))

args <- commandArgs(trailingOnly = TRUE)
.arg <- function(nm, d) {
    h <- grep(paste0("^--", nm, "="), args, value = TRUE)
    if (length(h) == 0L) d else sub(paste0("^--", nm, "="), "", h[1L])
}
REF     <- .arg("ref", "upstream/dev")
REPS    <- as.integer(.arg("reps", "3"))
SLOWER  <- as.numeric(.arg("threshold", "1.20"))
CAP     <- as.numeric(.arg("cap", "900"))
WORKERS <- as.integer(.arg("workers", Sys.getenv("GD_PAR_WORKERS", "8")))
REPORT  <- any(args == "--report")
PKG     <- normalizePath(".")
H5SRC   <- Sys.getenv("GD_BENCH_H5", "")
RESULTS <- file.path(PKG, "bench", "results")
RUNNER  <- file.path(PKG, "bench", "_runner.R")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)

.run <- function(pkg, label) {
    out <- file.path(RESULTS, paste0(label, ".tsv"))
    unlink(out)
    cat(sprintf("\n-- %s  (%s)\n", label, pkg))
    system2("Rscript", c(shQuote(RUNNER), shQuote(pkg), shQuote(label), REPS,
            shQuote(H5SRC), shQuote(out), CAP, WORKERS))
    if (!file.exists(out)) stop("no results from ", label, call. = FALSE)
}

.read <- function(label) {
    f <- file.path(RESULTS, paste0(label, ".tsv"))
    if (!file.exists(f))
        stop("no results for ", label, " -- run without --report first", call. = FALSE)
    d <- fread(f, header = FALSE, sep = "\t",
               col.names = c("tree", "case", "sec", "peak_gb"))
    d[!duplicated(case, fromLast = TRUE)]
}

if (!REPORT) {
    wt <- file.path(tempdir(), paste0("gd-ref-", gsub("[^A-Za-z0-9]", "", REF)))
    if (dir.exists(wt)) system2("git", c("worktree", "remove", "--force", shQuote(wt)))
    if (system2("git", c("worktree", "add", "-f", shQuote(wt), shQuote(REF)),
                stdout = FALSE, stderr = FALSE) != 0L)
        stop("could not create a worktree at ref '", REF, "'", call. = FALSE)
    on.exit(system2("git", c("worktree", "remove", "--force", shQuote(wt)),
                    stdout = FALSE, stderr = FALSE), add = TRUE)
    cat(sprintf(paste0("\n  ref     : %s\n  now     : %s\n  data    : %s\n",
                       "  reps    : %d   cap %.0fs   workers %d\n"),
        REF, PKG, if (nzchar(H5SRC)) H5SRC else "synthetic 4000 x 40000 @ 0.12",
        REPS, CAP, WORKERS))
    .run(wt, "REF"); .run(PKG, "NOW")
}

ref <- .read("REF"); now <- .read("NOW")
r <- merge(ref[, .(case, ref = sec, ref_gb = peak_gb)],
           now[, .(case, now = sec, now_gb = peak_gb)], by = "case", sort = FALSE)
FLOOR <- 0.005      # timer resolution
MEM_FLOOR <- 0.05   # system-wide sampling cannot resolve below background noise
r[, ratio := ref / now]
r[, below := pmax(ref, now) < FLOOR]
r[, flag := ifelse(below, "~0",
            ifelse(ratio < 1 / SLOWER, "SLOWER",
            ifelse(ratio > SLOWER, "faster", "")))]
.gb <- function(v) if (is.na(v) || v < MEM_FLOOR) "-" else sprintf("%.2f", v)
cat(sprintf("\n%-34s %8s %8s %7s %7s %7s  %s\n",
    "case", "ref s", "now s", "x", "ref GB", "now GB", ""))
for (i in seq_len(nrow(r)))
    cat(sprintf("%-34s %8.3f %8.3f %7s %7s %7s  %s\n", r$case[i], r$ref[i], r$now[i],
        if (r$below[i]) "-" else sprintf("%.2f", r$ratio[i]),
        .gb(r$ref_gb[i]), .gb(r$now_gb[i]), r$flag[i]))
n_slow <- sum(r$flag == "SLOWER", na.rm = TRUE)
cat(sprintf("\n  %d slower than %.2fx, %d faster, %d below %.0f ms\n  results: %s\n\n",
    n_slow, SLOWER, sum(r$flag == "faster", na.rm = TRUE), sum(r$below),
    FLOOR * 1000, RESULTS))
if (n_slow > 0L) quit(status = 1L)
