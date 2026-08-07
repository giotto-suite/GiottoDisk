#!/usr/bin/env Rscript
# Regression benchmark: compare this working tree against a git ref.
#
# Not a test. Nothing here asserts; it prints a table of ratios and flags what
# got slower. Lives outside the built package (`^bench$` in .Rbuildignore).
#
#   Rscript bench/regress.R                       # vs upstream/dev, synthetic
#   Rscript bench/regress.R --ref=HEAD~1
#   Rscript bench/regress.R --ref=e120c82 --reps=5
#   GD_BENCH_STORE=/path/to/store.parquet Rscript bench/regress.R
#
# Design notes, each of which is load-bearing:
#
# * PUBLIC API ONLY. The script runs against an older tree, which may not have
#   the internals -- or even the exported names -- that the current one does.
#   An earlier ad-hoc version died on `Giotto::reduceParam` not existing.
#
# * MEASURE VERBS, NOT INTERNALS, AT EXPLICIT CHAIN STATES. The one real
#   regression this harness was built for (22x on featStats) would have been
#   invisible to a micro-benchmark of the accumulator: the accumulator was
#   fine, and the verb was slow because the store's op placement changed which
#   execution path ran. So every case names the chain state it measures.
#
# * RATIOS, NOT BASELINES. Absolute times are not comparable across machines,
#   arrow versions, or background load. What is stable is the ratio between two
#   trees measured back to back on one machine, so there is no baselines file
#   to go stale.
#
# * ONE SCRIPT, TWO TREES. This file is always the runner, even when measuring
#   an old ref -- each side gets the same measurement code, only the package
#   under test differs. Never run the ref's own copy of this script.

suppressWarnings(suppressMessages({
    library(Matrix); library(data.table)
}))

# ---- args -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
.arg <- function(nm, default) {
    hit <- grep(paste0("^--", nm, "="), args, value = TRUE)
    if (length(hit) == 0L) return(default)
    sub(paste0("^--", nm, "="), "", hit[1L])
}
REF      <- .arg("ref", "upstream/dev")
REPS     <- as.integer(.arg("reps", "3"))
SLOWER   <- as.numeric(.arg("threshold", "1.20"))  # flag ratios past this
PKG      <- normalizePath(".")
STORE    <- Sys.getenv("GD_BENCH_STORE", "")

# ---- the cases --------------------------------------------------------------
# Each entry: label, and an expression evaluated in an env holding `pe` (raw
# store), `pn` (normalized: multiply + log), `u` (union, normalized), `sub`
# (small two-axis subset of pn), and `hvg` (feature id subset).
#
# Chain state is in the label because it is the axis that matters most: the
# same verb on a raw store and a normalized one can take entirely different
# execution paths.
# Expressions are STRINGS, not `quote()`d language objects: the list is
# serialized into the runner via deparse, and deparsing a language object yields
# its source text, so a `quote(FS(pe))` would come back as a bare `FS(pe)` and
# be evaluated when the runner builds the list -- before any timing starts.
CASES <- list(
    c("storeWrite (ingest)",             "WRITE()",                              "2"),
    c("processData libraryNorm  [raw]",  "LN(pe)",                               ""),
    c("processData logNorm      [norm]", "LG(pe)",                               ""),
    c("storeRead dgcmatrix      [norm]", "storeRead(sub, output = 'dgcmatrix')",  ""),
    c("analyzeData featStats    [raw]",  "FS(pe)",                               ""),
    c("analyzeData featStats    [norm]", "FS(pn)",                               ""),
    c("analyzeData cellStats    [norm]", "CS(pn)",                               ""),
    c("analyzeData cov_loess    [norm]", "HVF(pn)",                              ""),
    c("filterData               [norm]", "FILT(pn)",                             ""),
    c("reduceData random PCA    [norm]", "PCA(pn)",                              "2"),
    c("reduceData gram PCA      [norm]", "PCAG(pn)",                             "2"),
    c("analyzeData featStats  [union]",  "FS(u)",                                ""),
    c("analyzeData cov_loess  [union]",  "HVF(u)",                               ""),
    # End to end, from a raw store, with HVGs actually chosen by HVF rather
    # than faked. Catches a slowdown spread too thinly across steps to trip any
    # single threshold, and is the only case exercising the HVF -> PCA handoff.
    c("PIPELINE filter->norm->HVF->PCA",  "PIPELINE()",                           "1")
)

# ---- harness ----------------------------------------------------------------
# Written to a temp file and sourced inside a fresh R session per tree, so the
# two packages never share a namespace.
RUNNER <- tempfile(fileext = ".R")
writeLines(sprintf('
a <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(a[1], quiet = TRUE))
suppressWarnings(suppressMessages({library(Matrix); library(data.table)}))
LBL <- a[2]; REPS <- as.integer(a[3]); STORE <- a[4]

G <- 4000L; NC <- 40000L; DENS <- 0.12
set.seed(42)
if (nzchar(STORE)) {
    pe <- parquetExprStore(path = STORE)
    if (length(pe@cell_ids) == 0L)
        stop("GD_BENCH_STORE needs a store with cell_ids/feat_ids; ",
             "point at one saved with saveRDS instead, or unset it.")
    WRITE <- function() invisible(NULL)   # not meaningful on a supplied store
} else {
    m <- Matrix::rsparsematrix(G, NC, density = DENS,
         rand.x = function(n) as.double(rpois(n, 5L) + 1L))
    rownames(m) <- paste0("g", seq_len(G)); colnames(m) <- paste0("c", seq_len(NC))
    WRITE <- function()
        storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    pe <- WRITE()
}

LN  <- function(p) GiottoClass::processData(p, Giotto::normParam("library", scalefactor = 6e3))
LG  <- function(p) GiottoClass::processData(p, Giotto::normParam("log", base = 2, offset = 1))
FS  <- function(p) suppressWarnings(GiottoClass::analyzeData(p, Giotto::analyzeParam("feat_stats")))
CS  <- function(p) suppressWarnings(GiottoClass::analyzeData(p, Giotto::analyzeParam("cell_stats")))
HVF <- function(p) suppressWarnings(GiottoClass::analyzeData(p, Giotto::analyzeParam("cov_loess")))
FILT <- function(p) GiottoClass::filterData(p, Giotto::filterParam(
    expression_threshold = 1, feat_det_in_min_cells = 5, min_det_feats_per_cell = 5))
hvg <- rownames(pe)[seq_len(min(1000L, nrow(pe)))]
# `gramEigenPcaParam` is exported by GiottoDisk itself rather than routed
# through `Giotto::pcaParam()` (see R/pca-param.R) -- still public API, just a
# different package.
.pca <- function(p, feats, gram = FALSE) {
    prm <- if (gram) {
        gramEigenPcaParam(ncp = 20, feats_to_use = feats,
            center = TRUE, scale = FALSE)
    } else {
        Giotto::pcaParam("random", ncp = 20, feats_to_use = feats,
            center = TRUE, scale = FALSE, set_seed = TRUE, seed_number = 42L)
    }
    suppressWarnings(GiottoClass::reduceData(p, prm))
}
PCA  <- function(p) .pca(p, hvg)
PCAG <- function(p) .pca(p, hvg, gram = TRUE)

# Whole workflow from raw, HVGs taken from HVF output. `reps = 1` in the case
# list -- this is the slowest entry and its job is a trend line, not precision.
# `filterData` returns the surviving ids (masks), not a narrowed store -- the
# caller applies them. That is the shape a gobject-level filterGiotto wraps.
PIPELINE <- function() {
    m <- FILT(pe)
    keep <- pe[m$feats_keep, m$cells_keep]
    p  <- LG(LN(keep))
    st <- HVF(p)
    data.table::setorder(st, -cov_diff)
    .pca(p, utils::head(st$feats, 1000L))
}

pn  <- LG(LN(pe))
sub <- pn[seq_len(min(80L, nrow(pn))), seq_len(min(80L, ncol(pn)))]
u   <- local({
    half <- seq_len(ncol(pe) %%/%% 2L)
    p2 <- if (nzchar(STORE)) pe[, half] else
        storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
                   {mm <- m[, half]; colnames(mm) <- paste0("z", seq_along(half)); mm})
    LG(LN(unionParquetExprStore(list(pe, p2))))
})

CASES <- %s
for (cs in CASES) {
    reps <- if (nzchar(cs[3L])) as.integer(cs[3L]) else REPS
    ex <- parse(text = cs[2L])
    t <- vapply(seq_len(reps),
                function(i) system.time(eval(ex))[["elapsed"]], numeric(1L))
    cat(sprintf("%%s\\t%%s\\t%%.4f\\n", LBL, cs[1L], median(t)))
}
', paste(deparse(CASES), collapse = "\n")), RUNNER)

.run <- function(pkg, label) {
    out <- suppressWarnings(system2("Rscript",
        c(shQuote(RUNNER), shQuote(pkg), shQuote(label), REPS, shQuote(STORE)),
        stdout = TRUE, stderr = FALSE))
    rows <- grep(paste0("^", label, "\t"), out, value = TRUE)
    if (length(rows) == 0L) {
        cat(paste(utils::tail(out, 20L), collapse = "\n"), "\n")
        stop("no timings from ", label, " -- see output above", call. = FALSE)
    }
    d <- data.table::fread(text = paste(rows, collapse = "\n"),
                           header = FALSE, sep = "\t")
    data.table::setnames(d, c("tree", "case", "sec"))
    d[]
}

# ---- run both sides ---------------------------------------------------------
wt <- file.path(tempdir(), paste0("gd-ref-", gsub("[^A-Za-z0-9]", "", REF)))
if (dir.exists(wt)) system2("git", c("worktree", "remove", "--force", shQuote(wt)))
st <- system2("git", c("worktree", "add", "-f", shQuote(wt), shQuote(REF)),
              stdout = FALSE, stderr = FALSE)
if (st != 0L) stop("could not create a worktree at ref '", REF, "'", call. = FALSE)
on.exit(system2("git", c("worktree", "remove", "--force", shQuote(wt)),
                stdout = FALSE, stderr = FALSE), add = TRUE)

cat(sprintf("\n  ref  : %s\n  now  : %s\n  data : %s\n  reps : %d\n\n",
    REF, PKG, if (nzchar(STORE)) STORE else "synthetic 4000 x 40000 @ 0.12",
    REPS))

ref <- .run(wt,  "REF")
now <- .run(PKG, "NOW")

# ---- report -----------------------------------------------------------------
r <- merge(ref[, .(case, ref = sec)], now[, .(case, now = sec)],
           by = "case", sort = FALSE)
# Below the timer's useful resolution, a ratio is noise -- and some cases are
# meant to sit there. `logNorm` appends a record and does no data pass, so ~0 is
# the correct answer and a regression would show up as it becoming measurable.
FLOOR <- 0.005
r[, ratio := ref / now]                      # >1 faster, <1 slower
r[, below := pmax(ref, now) < FLOOR]
r[, flag := ifelse(below, "~0",
            ifelse(ratio < 1 / SLOWER, "SLOWER",
            ifelse(ratio > SLOWER, "faster", "")))]
cat(sprintf("%-34s %8s %8s %8s  %s\n", "case", "ref", "now", "x", ""))
for (i in seq_len(nrow(r))) {
    cat(sprintf("%-34s %8.3f %8.3f %8s  %s\n", r$case[i], r$ref[i], r$now[i],
        if (r$below[i]) "-" else sprintf("%.2f", r$ratio[i]), r$flag[i]))
}
n_slow <- sum(r$flag == "SLOWER", na.rm = TRUE)
cat(sprintf("\n  %d slower than %.2fx, %d faster, %d below %.0f ms\n\n",
            n_slow, SLOWER, sum(r$flag == "faster", na.rm = TRUE),
            sum(r$below), FLOOR * 1000))
if (n_slow > 0L) quit(status = 1L)
