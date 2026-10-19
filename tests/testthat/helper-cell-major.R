# Shared check for the cell-major layout every expression store must have
# (AGENTS.md, "window the CELL axis"; adr/0017): each part file holds one
# range of cells, and the ranges do not overlap, so a cell window lowers to a
# `row_id` range that prunes the files and row groups outside it.

# Per part file under `dir`, the row_id range it holds, ordered by range start.
.part_ranges <- function(dir) {
    parts <- list.files(dir, pattern = "[.]parquet$", recursive = TRUE,
                        full.names = TRUE)
    r <- do.call(rbind, lapply(parts, function(f) {
        x <- range(arrow::read_parquet(f, col_select = "row_id")$row_id)
        data.frame(file = basename(f), lo = x[1L], hi = x[2L])
    }))
    r[order(r$lo, r$hi), , drop = FALSE]
}

# `strict = FALSE` allows two neighbouring parts to share their boundary cell,
# as a reader that cuts batches by record count rather than by cell does.
.expect_cell_major <- function(dir, strict = TRUE) {
    pr <- .part_ranges(dir)
    if (nrow(pr) > 1L) {
        nxt <- pr$lo[-1L]; prv <- pr$hi[-nrow(pr)]
        if (strict) expect_true(all(nxt > prv)) else expect_true(all(nxt >= prv))
    }
    invisible(pr)
}
