# Case list, sourced by both the driver (for nothing but documentation) and the
# runner (which executes them). Sourced rather than interpolated so the
# expressions stay real code -- an earlier version serialized them through
# `deparse()`, which strips `quote()` and silently evaluated every case at
# list-construction time, before any timing began.
#
# Columns: label, expression, reps ("" = default), real-data ("no" = skip).
#
# Chain state is in the label because it is the axis that matters most: the same
# verb on a raw store and a normalized one can take entirely different execution
# paths. The 22x featStats regression this harness exists for showed up only in
# the [norm] case; [raw] was at parity.
#
# The union cases are skipped on real data -- they measure plan composition
# across substores, which the synthetic run already covers, and building a
# second substore from Atera means writing ~150M rows for no extra signal.
BENCH_CASES <- list(
    list("storeWrite (ingest)",             quote(INGEST()),                        1L, TRUE),
    list("processData libraryNorm  [raw]",  quote(LN(pe)),                          NA, TRUE),
    list("processData logNorm      [raw]",  quote(LG(pe)),                          NA, TRUE),
    list("storeRead dgcmatrix      [norm]", quote(storeRead(sub, output = "dgcmatrix")), NA, TRUE),
    list("analyzeData featStats    [raw]",  quote(FS(pe)),                          NA, TRUE),
    list("analyzeData cellStats    [raw]",  quote(CS(pe)),                          NA, TRUE),
    list("analyzeData featStats    [norm]", quote(FS(pn)),                          NA, TRUE),
    list("analyzeData cellStats    [norm]", quote(CS(pn)),                          NA, TRUE),
    list("analyzeData cov_loess    [norm]", quote(HVF(pn)),                         NA, TRUE),
    list("filterData               [norm]", quote(FILT(pn)),                        NA, TRUE),
    list("reduceData random PCA    [norm]", quote(PCA(pn)),                         2L, TRUE),
    list("reduceData gram PCA      [norm]", quote(PCAG(pn)),                        2L, TRUE),
    list("analyzeData featStats  [union]",  quote(FS(u)),                           NA, FALSE),
    list("analyzeData cov_loess  [union]",  quote(HVF(u)),                          NA, FALSE),
    list("PIPELINE filter->norm->HVF->PCA", quote(PIPELINE()),                      1L, TRUE)
)
