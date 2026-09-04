# Measured side of the benchmark. Invoked once per tree by bench/regress.R:
#
#   Rscript bench/_runner.R <pkg> <label> <reps> <h5|""> <out.tsv> <cap_s> <workers> [case_regex]
#
# A real file rather than a template written at runtime: the previous version
# built this with sprintf + writeLines, which meant %% escaping throughout, no
# syntax checking, and no way to run it standalone when something broke.
#
# Everything here goes through the public API. Each side runs a different tree,
# which may not have the internals -- or the exported names -- the other does.

a <- commandArgs(trailingOnly = TRUE)
PKG <- a[1]; LBL <- a[2]; REPS <- as.integer(a[3]); H5SRC <- a[4]
OUT <- a[5]; CAP <- as.numeric(a[6]); WORKERS <- as.integer(a[7])
CASESEL <- if (length(a) >= 8L) a[8] else ""

suppressMessages(pkgload::load_all(PKG, quiet = TRUE))
suppressWarnings(suppressMessages({library(Matrix); library(data.table)}))
source(file.path(dirname(sub("^--file=", "", grep("^--file=",
    commandArgs(FALSE), value = TRUE)[1L])), "_cases.R"))

# ---- parallelism ------------------------------------------------------------
# Not optional. `.par_workers()` resolves from `options("giottodisk.par_workers")`
# or a non-uniprocess future plan, and defaults to 1 -- so a bench that sets
# neither measures the serial path. The h5 ingest comment records fork at 8
# workers as 9.0 s / +16.1 GB against 17.8 s via mirai, i.e. roughly 2x, so
# leaving this unset benchmarks a configuration nobody runs.
#
# A plan is set as well as the option, since consumers reached through
# lapply_flex read the plan rather than the option.
options(giotto.warn_sequential = FALSE)
options(giottodisk.par_workers = WORKERS)
if (WORKERS > 1L && requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.mirai", quietly = TRUE)) {
    future::plan(future.mirai::mirai_multisession, workers = WORKERS)
    # Warm the daemons and attach packages through the framework hook, not
    # library() in the worker body -- that hides the dependency and bypasses
    # R CMD check. Cost is paid once, outside any timed case.
    invisible(future.apply::future_lapply(seq_len(WORKERS), function(i) NULL,
        future.packages = c("arrow", "data.table", "Matrix", "hdf5r"),
        future.seed = NULL))
} else {
    if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
}
on.exit({
    if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
    if (requireNamespace("mirai", quietly = TRUE)) try(mirai::daemons(0), silent = TRUE)
}, add = TRUE)

# ---- peak footprint ---------------------------------------------------------
# Samples SYSTEM free memory and takes the minimum. arrow allocates in C++, so
# gc()-based accounting misses most of what these verbs do.
.free_cmd <- if (Sys.info()[["sysname"]] == "Darwin") {
    paste("vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}",
          "/Pages speculative/{s=$3}/Pages purgeable/{p=$3}",
          "END{gsub(/[^0-9]/,\"\",f);gsub(/[^0-9]/,\"\",i);",
          "gsub(/[^0-9]/,\"\",s);gsub(/[^0-9]/,\"\",p);print (f+i+s+p)*16384}'")
} else if (file.exists("/proc/meminfo")) {
    "awk '/MemAvailable/{print $2*1024}' /proc/meminfo"
} else NA_character_

.timed <- function(ex, reps, env) {
    secs <- rep(NA_real_, reps); peaks <- rep(NA_real_, reps)
    for (i in seq_len(reps)) {
        invisible(gc(full = TRUE, reset = TRUE))
        tf <- tempfile(); sn <- tempfile(); file.create(sn)
        base <- NA_real_
        if (!is.na(.free_cmd)) {
            base <- as.numeric(system(.free_cmd, intern = TRUE))
            system(sprintf("while [ -f %s ]; do %s >> %s; sleep 0.2; done",
                           sn, .free_cmd, tf),
                   wait = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
        }
        t0 <- proc.time()[["elapsed"]]
        ok <- tryCatch({
            setTimeLimit(elapsed = CAP, transient = TRUE)
            eval(ex, env); setTimeLimit(); TRUE
        }, error = function(e) { setTimeLimit(); message("  ! ", conditionMessage(e)); FALSE })
        el <- proc.time()[["elapsed"]] - t0
        unlink(sn); Sys.sleep(0.3)
        peaks[i] <- tryCatch((base - min(as.numeric(readLines(tf)))) / 1024^3,
                             warning = function(w) NA_real_,
                             error = function(e) NA_real_)
        unlink(tf)
        if (!ok) break
        secs[i] <- el
    }
    c(sec = suppressWarnings(median(secs, na.rm = TRUE)),
      peak_gb = suppressWarnings(max(peaks, na.rm = TRUE)))
}

.emit <- function(case, r) {
    cat(sprintf("%s\t%s\t%.4f\t%.3f\n", LBL, case, r[["sec"]], r[["peak_gb"]]),
        file = OUT, append = TRUE)
    cat(sprintf("  %-34s %8.3f s  %7.2f GB\n", case, r[["sec"]], r[["peak_gb"]]))
}

# ---- fixtures ---------------------------------------------------------------
E <- new.env(parent = globalenv())
E$h5mode <- nzchar(H5SRC)
local({
    set.seed(42)
    if (E$h5mode) {
        # Real data starts where a real workflow starts: at ingest. A store saved
        # from an older tree will not deserialize against a changed class, and
        # the Input path yields real ids and makes storeWrite a measured verb.
        #
        # `batch_cells` matters: the 250,000 default exceeds this dataset's cell
        # count, so the whole matrix would arrive in one hyperslab. 25,000 is
        # what was tuned for it -- see adr/0007 on read shape.
        E$INGEST <- function() {
            inp <- tenxH5Input(H5SRC, feature_id_col = 2L, batch_cells = 25000L)
            pe_all <- storeWrite(
                parquetExprStore(path = tempfile(fileext = ".parquet")), inp)
            h <- hdf5r::H5File$new(H5SRC, mode = "r")
            root <- names(h)[1L]
            ftype <- as.character(h[[paste0(root, "/features/feature_type")]][])
            h$close_all()
            pe_all[which(ftype == "Gene Expression"), ]
        }
    } else {
        G <- 4000L; NC <- 40000L
        m <- Matrix::rsparsematrix(G, NC, density = 0.12,
             rand.x = function(n) as.double(rpois(n, 5L) + 1L))
        rownames(m) <- paste0("g", seq_len(G))
        colnames(m) <- paste0("c", seq_len(NC))
        E$m <- m
        E$INGEST <- function()
            storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")), m)
    }
})

# scalefactor 1e4 is the original's `target_sum`; base exp(1) reproduces its
# plain `log1p()` (Giotto's house default is base 2, which is NOT what the
# reference profile was recorded with).
E$LN  <- function(p) GiottoClass::processData(p, Giotto::normParam("library", scalefactor = 1e4))
E$LG  <- function(p) GiottoClass::processData(p, Giotto::normParam("log", base = exp(1), offset = 1))
E$FS  <- function(p) suppressWarnings(GiottoClass::analyzeData(p, Giotto::analyzeParam("feat_stats")))
E$CS  <- function(p) suppressWarnings(GiottoClass::analyzeData(p, Giotto::analyzeParam("cell_stats")))
E$HVF <- function(p) suppressWarnings(GiottoClass::analyzeData(p, Giotto::analyzeParam("cov_loess")))
E$FILT <- function(p) GiottoClass::filterData(p, Giotto::filterParam(
    expression_threshold = 1, feat_det_in_min_cells = 5, min_det_feats_per_cell = 5))

# Ingest is the fixture AND the first timed case -- one pass, not three. At
# Atera scale a spare ingest is minutes and tens of GB of temp parquet.
# Ingest always runs -- it is the fixture -- but is only *reported* when it
# passes the case filter.
ing <- BENCH_CASES[[1L]]
r_ing <- .timed(quote(pe <- INGEST()), 1L, E)
if (!nzchar(CASESEL) || grepl(CASESEL, ing[[1L]])) .emit(ing[[1L]], r_ing)
cat(sprintf("  store: %s feats x %s cells\n",
    format(nrow(E$pe), big.mark = ","), format(ncol(E$pe), big.mark = ",")))

E$pn  <- E$LG(E$LN(E$pe))
E$sub <- E$pn[seq_len(min(80L, nrow(E$pn))), seq_len(min(80L, ncol(E$pn)))]
# Parameters aligned with the original Atera harness (NCP 30, HVG_N 2000) so
# the overlapping cases are comparable to its recorded profile.csv rather than
# being a fourth incompatible configuration.
E$NCP <- 30L; E$HVG_N <- 2000L

# Grouping for the windowed accumulator. Every other featStats case is
# ungrouped, and ungrouped never windows -- so without this the harness cannot
# see `.pe_accum_acero_windowed()` at all, which is the path a grouping puts an
# O(nonzeros) join in front of (adr/0011).
#
# Named by cell id, not positional: unnamed groups match by order and warn,
# and cell metadata is not guaranteed to share the store's order.
#
# 12 groups is a plausible cluster count and it sets the aggregate's width
# (features x groups). Derived outside any timed case.
E$N_GRP <- 12L
E$groups <- stats::setNames(
    factor(rep(paste0("clus", seq_len(E$N_GRP)), length.out = ncol(E$pn))),
    colnames(E$pn))
E$FSG <- function(p) suppressWarnings(GiottoClass::analyzeData(
    p, Giotto::analyzeParam("feat_stats"), groups = E$groups))

# HVGs come from HVF output, not `rownames(pe)[1:n]`. Feature selection picks
# the DENSEST rows on purpose (adr/0008), so an arbitrary prefix reads a
# sparser band than any real workflow would and makes PCA look faster than it
# is. Derived once here, outside any timed case.
E$hvg <- local({
    st <- E$HVF(E$pn)
    data.table::setorder(st, -cov_diff)
    utils::head(st$feats, min(E$HVG_N, nrow(E$pn)))
})
E$.pca <- function(p, feats, gram = FALSE) {
    # gramEigenPcaParam is exported by GiottoDisk itself rather than routed
    # through Giotto::pcaParam() -- see R/pca-param.R. Still public API.
    prm <- if (gram) {
        gramEigenPcaParam(ncp = E$NCP, feats_to_use = feats,
            center = TRUE, scale = FALSE)
    } else {
        Giotto::pcaParam("random", ncp = E$NCP, feats_to_use = feats,
            center = TRUE, scale = FALSE, set_seed = TRUE, seed_number = 42L)
    }
    suppressWarnings(GiottoClass::reduceData(p, prm))
}
E$PCA  <- function(p) E$.pca(p, E$hvg)
E$PCAG <- function(p) E$.pca(p, E$hvg, gram = TRUE)
# `HVG_N` and gram, not 1000 features and randomized: the original ran
# sc_pca(method = "gram") over n_top = 2000. The standalone PCA cases already
# use both, so hardcoding otherwise here made the end-to-end case the one
# configuration that matched nothing.
E$PIPELINE <- function() {
    mk <- E$FILT(E$pe)
    keep <- E$pe[mk$feats_keep, mk$cells_keep]
    p <- E$LG(E$LN(keep)); st <- E$HVF(p)
    data.table::setorder(st, -cov_diff)
    E$.pca(p, utils::head(st$feats, min(E$HVG_N, nrow(p))), gram = TRUE)
}
if (!E$h5mode) E$u <- local({
    half <- seq_len(ncol(E$pe) %/% 2L)
    p2 <- storeWrite(parquetExprStore(path = tempfile(fileext = ".parquet")),
                     {mm <- E$m[, half]; colnames(mm) <- paste0("z", seq_along(half)); mm})
    E$LG(E$LN(unionParquetExprStore(list(E$pe, p2))))
})

for (cs in BENCH_CASES[-1L]) {
    if (E$h5mode && !isTRUE(cs[[4L]])) next
    if (nzchar(CASESEL) && !grepl(CASESEL, cs[[1L]])) next
    reps <- if (is.na(cs[[3L]])) REPS else cs[[3L]]
    .emit(cs[[1L]], .timed(cs[[2L]], reps, E))
}
