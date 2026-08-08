#!/usr/bin/env Rscript
# Driver: run bench/_runner.R against two trees and report the difference.
#
#   Rscript bench/regress.R                        # every dataset, every case
#   Rscript bench/regress.R --data=synthetic       # one dataset
#   Rscript bench/regress.R --data=atera
#   Rscript bench/regress.R --cases=PCA            # one slice of cases (regex)
#   Rscript bench/regress.R --cases='featStats|cellStats' --reps=5
#   Rscript bench/regress.R --ref=HEAD~1 --workers=1
#   Rscript bench/regress.R --ref=A --now=B        # pin both ends (reproducible)
#   Rscript bench/regress.R --report               # re-print, no work
#
# `--ref` compares against the working tree. `--now` replaces the working tree
# with a second worktree, so a comparison stays reproducible after the branch
# moves on -- use it when citing numbers anywhere durable. The harness itself
# (_runner.R, _cases.R) always comes from the working tree either way, so both
# trees are measured by the same code and only the package differs.
#
# `atera` needs GD_BENCH_H5 pointing at a cell_feature_matrix.h5. Without it
# the default run does synthetic only and says so, so the default works
# anywhere.
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
NOWREF  <- .arg("now", "")   # "" = the working tree
DATA    <- .arg("data", "all")
CASESEL <- .arg("cases", "")
REPS    <- as.integer(.arg("reps", "3"))
SLOWER  <- as.numeric(.arg("threshold", "1.20"))
CAP     <- as.numeric(.arg("cap", "300"))   # matches the original harness
WORKERS <- as.integer(.arg("workers", Sys.getenv("GD_PAR_WORKERS", "8")))
REPORT  <- any(args == "--report")
PKG     <- normalizePath(".")
H5SRC   <- Sys.getenv("GD_BENCH_H5", "")
RESULTS <- file.path(PKG, "bench", "results")
RUNNER  <- file.path(PKG, "bench", "_runner.R")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)

# Datasets are separately runnable and all run by default. `atera` is skipped
# rather than erroring when GD_BENCH_H5 is unset, so `--data=all` is a sensible
# default on a machine that does not have the file.
ALL_DATA <- c("synthetic", "atera")
DATASETS <- if (identical(DATA, "all")) ALL_DATA else strsplit(DATA, ",")[[1L]]
if (!all(DATASETS %in% ALL_DATA))
    stop("--data must be one of: ", toString(ALL_DATA), ", all", call. = FALSE)
if ("atera" %in% DATASETS && !nzchar(H5SRC)) {
    if (identical(DATA, "all")) {
        message("  note: skipping atera -- GD_BENCH_H5 is not set")
        DATASETS <- setdiff(DATASETS, "atera")
    } else {
        stop("--data=atera needs GD_BENCH_H5 set to a cell_feature_matrix.h5",
             call. = FALSE)
    }
}
.src <- function(ds) if (identical(ds, "atera")) H5SRC else ""

.tsv <- function(label, ds) file.path(RESULTS, sprintf("%s-%s.tsv", label, ds))

.run <- function(pkg, label, ds) {
    out <- .tsv(label, ds)
    unlink(out)
    cat(sprintf("\n-- %s / %s  (%s)\n", ds, label, pkg))
    system2("Rscript", c(shQuote(RUNNER), shQuote(pkg), shQuote(label), REPS,
            shQuote(.src(ds)), shQuote(out), CAP, WORKERS, shQuote(CASESEL)))
    if (!file.exists(out)) stop("no results from ", label, "/", ds, call. = FALSE)
}

.read <- function(label, ds) {
    f <- .tsv(label, ds)
    if (!file.exists(f)) return(NULL)
    d <- fread(f, header = FALSE, sep = "\t",
               col.names = c("tree", "case", "sec", "peak_gb"))
    d[!duplicated(case, fromLast = TRUE)]
}

if (!REPORT) {
    # Both sides are just package paths, so either can be a worktree. With
    # `--now` unset the working tree is NOW, which is the development case;
    # supplying it pins both ends, which is what a bug report needs -- the
    # comparison then reproduces from any checkout state instead of silently
    # depending on where the working tree happens to sit.
    #
    # `prune` sweeps registrations whose directory is already gone -- a run
    # killed mid-flight leaves the metadata behind even though tempdir() went
    # with the session.
    created <- character(0)
    .wt_clean <- function() {
        for (p in created)
            system2("git", c("worktree", "remove", "--force", shQuote(p)),
                    stdout = FALSE, stderr = FALSE)
        system2("git", c("worktree", "prune"), stdout = FALSE, stderr = FALSE)
    }
    # Role is in the path so `--ref=X --now=X` does not collide on one dir.
    .wt_add <- function(ref, role) {
        p <- file.path(tempdir(),
                       sprintf("gd-%s-%s", role, gsub("[^A-Za-z0-9]", "", ref)))
        system2("git", c("worktree", "remove", "--force", shQuote(p)),
                stdout = FALSE, stderr = FALSE)
        system2("git", c("worktree", "prune"), stdout = FALSE, stderr = FALSE)
        if (system2("git", c("worktree", "add", "-f", shQuote(p), shQuote(ref)),
                    stdout = FALSE, stderr = FALSE) != 0L)
            stop("could not create a worktree at ref '", ref, "'", call. = FALSE)
        created <<- c(created, p)
        p
    }
    .wt_clean()
    # `finally`, not on.exit(): this is top-level script code, so there is no
    # function frame whose exit on.exit() could fire on. It silently registered
    # nothing and leaked one worktree per run.
    tryCatch({
        wt_ref <- .wt_add(REF, "ref")
        wt_now <- if (nzchar(NOWREF)) .wt_add(NOWREF, "now") else PKG
        cat(sprintf(paste0("\n  ref     : %s\n  now     : %s\n  data    : %s\n",
                           "  cases   : %s\n  reps    : %d   cap %.0fs   workers %d\n"),
            REF, if (nzchar(NOWREF)) NOWREF else PKG, toString(DATASETS),
            if (nzchar(CASESEL)) CASESEL else "all", REPS, CAP, WORKERS))
        for (ds in DATASETS) { .run(wt_ref, "REF", ds); .run(wt_now, "NOW", ds) }
    }, finally = .wt_clean())
}

FLOOR <- 0.005      # timer resolution
MEM_FLOOR <- 0.05   # system-wide sampling cannot resolve below background noise
.gb <- function(v) if (is.na(v) || v < MEM_FLOOR) "-" else sprintf("%.2f", v)

n_slow_total <- 0L
shown <- character(0)
for (ds in if (REPORT) ALL_DATA else DATASETS) {
    ref <- .read("REF", ds); now <- .read("NOW", ds)
    if (is.null(ref) || is.null(now)) next
    shown <- c(shown, ds)
    r <- merge(ref[, .(case, ref = sec, ref_gb = peak_gb)],
               now[, .(case, now = sec, now_gb = peak_gb)], by = "case", sort = FALSE)
    r[, ratio := ref / now]
    r[, below := pmax(ref, now) < FLOOR]
    r[, flag := ifelse(below, "~0",
                ifelse(ratio < 1 / SLOWER, "SLOWER",
                ifelse(ratio > SLOWER, "faster", "")))]
    cat(sprintf("\n=== %s %s\n", ds, strrep("=", max(0L, 60L - nchar(ds))))) 
    cat(sprintf("%-34s %8s %8s %7s %7s %7s  %s\n",
        "case", "ref s", "now s", "x", "ref GB", "now GB", ""))
    for (i in seq_len(nrow(r)))
        cat(sprintf("%-34s %8.3f %8.3f %7s %7s %7s  %s\n", r$case[i],
            r$ref[i], r$now[i],
            if (r$below[i]) "-" else sprintf("%.2f", r$ratio[i]),
            .gb(r$ref_gb[i]), .gb(r$now_gb[i]), r$flag[i]))
    ns <- sum(r$flag == "SLOWER", na.rm = TRUE)
    n_slow_total <- n_slow_total + ns
    cat(sprintf("  %d slower than %.2fx, %d faster\n",
        ns, SLOWER, sum(r$flag == "faster", na.rm = TRUE)))
}
if (length(shown) == 0L)
    stop("no results to report -- run without --report first", call. = FALSE)
cat(sprintf("\n  %d regression(s) across %s\n  results: %s\n\n",
    n_slow_total, toString(shown), RESULTS))
if (n_slow_total > 0L) quit(status = 1L)
