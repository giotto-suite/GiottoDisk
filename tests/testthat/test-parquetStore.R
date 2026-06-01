test_that("parquetStore: unwritten store does not exist", {
    ps <- parquetStore()
    expect_false(storeExists(ps))
})

test_that("parquetStore: write and existence", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_true(storeExists(ps))
})

test_that("parquetStore: dimensions after write", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_equal(nrow(ps), 32)
    expect_equal(ncol(ps), ncol(mtcars))
    expect_equal(dim(ps), c(32, ncol(mtcars)))
})

test_that("parquetStore: special cols hidden from colnames", {
    ps <- parquetStore() |> storeWrite(mtcars)
    cn <- colnames(ps)
    expect_false(any(c("row_index", "source_id") %in% cn))
    expect_true(all(colnames(mtcars) %in% cn))
})

test_that("parquetStore: storeRead query output is lazy arrow object", {
    ps <- parquetStore() |> storeWrite(mtcars)
    q <- storeRead(ps)
    expect_true(inherits(q, "ArrowObject"))
})

test_that("parquetStore: storeRead tibble roundtrip", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- storeRead(ps, output = "tibble")
    expect_s3_class(tbl, "data.frame")
    expect_equal(nrow(tbl), 32)
    expect_equal(sort(colnames(tbl)), sort(colnames(mtcars)))
})

test_that("parquetStore: subset op is recorded lazily", {
    ps <- parquetStore() |> storeWrite(mtcars)
    ps2 <- subset(ps, cyl == 4)
    expect_length(ps2@ops, 1L)
    # original store unaffected
    expect_length(ps@ops, 0L)
})

test_that("parquetStore: subset filters correctly at read time", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- subset(ps, cyl == 4) |> storeRead(output = "tibble")
    expect_true(all(tbl$cyl == 4))
    expect_equal(nrow(tbl), sum(mtcars$cyl == 4))
})

test_that("parquetStore: subset with local variable is inlined", {
    ps <- parquetStore() |> storeWrite(mtcars)
    target_cyl <- 6
    tbl <- subset(ps, cyl == target_cyl) |> storeRead(output = "tibble")
    expect_true(all(tbl$cyl == 6))
})

test_that("parquetStore: head limits rows", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- head(ps, 5) |> storeRead(output = "tibble")
    expect_lte(nrow(tbl), 5)
})

test_that("parquetStore: [,j] column selection", {
    ps <- parquetStore() |> storeWrite(mtcars)
    ps2 <- ps[, c("mpg", "cyl")]
    tbl <- storeRead(ps2, output = "tibble")
    expect_equal(sort(colnames(tbl)), sort(c("mpg", "cyl")))
})

test_that("parquetStore: nrow returns numeric not integer", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_type(nrow(ps), "double")
})

# [-join nomatch contract ####
# Two small stores sharing an `id` key. `target` has every id; `y` has a subset.
# Used by the join contract tests below.
.make_join_stores <- function() {
    target <- parquetStore() |>
        storeWrite(data.frame(id = letters[1:5], x = 1:5,
            stringsAsFactors = FALSE))
    y <- parquetStore() |>
        storeWrite(data.frame(id = letters[c(1, 3, 5)], y = c(10, 30, 50),
            stringsAsFactors = FALSE))
    list(target = target, y = y)
}

test_that("[-join: default queues op with nomatch='left'", {
    s <- .make_join_stores()
    j <- s$target[s$y, on = "id"]
    expect_length(j@ops, 1L)
    expect_equal(j@ops[[1L]]$type, "join")
    expect_equal(j@ops[[1L]]$nomatch, "left")
})

test_that("[-join: nomatch=NULL queues op with nomatch='inner'", {
    s <- .make_join_stores()
    j <- s$target[s$y, on = "id", nomatch = NULL]
    expect_equal(j@ops[[1L]]$nomatch, "inner")
})

test_that("[-join: nomatch=NA queues op with nomatch='left'", {
    s <- .make_join_stores()
    j <- s$target[s$y, on = "id", nomatch = NA]
    expect_equal(j@ops[[1L]]$nomatch, "left")
})

test_that("[-join: nomatch other than NULL/NA errors", {
    s <- .make_join_stores()
    expect_error(s$target[s$y, on = "id", nomatch = 0],
        "must be `NULL` \\(inner\\) or `NA` \\(left\\)")
    expect_error(s$target[s$y, on = "id", nomatch = "left"],
        "must be `NULL` \\(inner\\) or `NA` \\(left\\)")
})

test_that("[-join: missing `on=` errors", {
    s <- .make_join_stores()
    expect_error(s$target[s$y], "on")
})

test_that("[-join: queueing does not mutate the original target", {
    s <- .make_join_stores()
    j <- s$target[s$y, on = "id"]
    expect_length(s$target@ops, 0L)
    expect_length(j@ops, 1L)
})

test_that("[-join: default (left) preserves x rows with NA fill on miss", {
    s <- .make_join_stores()
    tbl <- s$target[s$y, on = "id"] |>
        storeRead(output = "tibble") |>
        dplyr::arrange(id)
    expect_equal(nrow(tbl), 5L)
    expect_equal(tbl$id, letters[1:5])
    expect_equal(tbl$y, c(10, NA, 30, NA, 50))
})

test_that("[-join: nomatch=NULL (inner) drops unmatched x rows", {
    s <- .make_join_stores()
    tbl <- s$target[s$y, on = "id", nomatch = NULL] |>
        storeRead(output = "tibble") |>
        dplyr::arrange(id)
    expect_equal(nrow(tbl), 3L)
    expect_equal(tbl$id, letters[c(1, 3, 5)])
    expect_equal(tbl$y, c(10, 30, 50))
})

test_that("[-join: y's special cols are dropped from result except join keys", {
    s <- .make_join_stores()
    tbl <- s$target[s$y, on = "id", nomatch = NULL] |>
        storeRead(output = "tibble", omit_internals = FALSE)
    # row_index / source_id from y should not leak; only the join key + value col
    expect_false("source_id.y" %in% colnames(tbl))
    expect_true("id" %in% colnames(tbl))
    expect_true("y" %in% colnames(tbl))
})


# col-select + filter order-of-ops composition ####
# `[, j]` narrows @fields independently of the @ops queue. Filter ops in @ops
# must still be able to reference cols that `[, j]` dropped from the visible
# schema. The lazy_fields layer widens the projection at storeRead time so
# filter exprs resolve; .pbase_storeread_processing narrows back to user
# fields after ops run.

test_that("col-select + subset on different col: [, j] first, tibble", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- ps[, "mpg"] |> subset(cyl == 4) |> storeRead(output = "tibble")
    expect_setequal(names(tbl), "mpg")
    expect_equal(nrow(tbl), sum(mtcars$cyl == 4))
})

test_that("col-select + subset on different col: subset first, tibble", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- subset(ps, cyl == 4)[, "mpg"] |> storeRead(output = "tibble")
    expect_setequal(names(tbl), "mpg")
    expect_equal(nrow(tbl), sum(mtcars$cyl == 4))
})

test_that("col-select + multi-col filter (>1 dropped col)", {
    ps <- parquetStore() |> storeWrite(mtcars)
    tbl <- ps[, "mpg"] |> subset(cyl == 4 & hp > 50) |>
        storeRead(output = "tibble")
    expect_setequal(names(tbl), "mpg")
    expect_equal(nrow(tbl), sum(mtcars$cyl == 4 & mtcars$hp > 50))
})

test_that("programmatic subset (quote=FALSE) resolves dropped cols", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expr <- quote(cyl == 4)
    tbl <- ps[, "mpg"] |> subset(expr, quote = FALSE) |>
        storeRead(output = "tibble")
    expect_setequal(names(tbl), "mpg")
    expect_equal(nrow(tbl), sum(mtcars$cyl == 4))
})

test_that("col-select + subset composes through query output", {
    ps <- parquetStore() |> storeWrite(mtcars)
    q <- ps[, "mpg"] |> subset(cyl == 4) |> storeRead(output = "query")
    df <- dplyr::collect(q)
    expect_setequal(names(df), "mpg")
    expect_equal(nrow(df), sum(mtcars$cyl == 4))
})

test_that("col-select + subset composes with head (no row_index warning)", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_warning(
        tbl <- ps[, "mpg"] |> subset(cyl == 4) |> head(2) |>
            storeRead(output = "tibble"),
        NA  # no warning expected
    )
    expect_lte(nrow(tbl), 2L)
    expect_setequal(names(tbl), "mpg")
})

test_that("subset() can reference any disk col regardless of @fields", {
    # `.inline_local_vars` uses effective_schema, not colnames -- so a
    # [, j] narrow does not break subset()'s free-var resolution.
    ps <- parquetStore() |> storeWrite(mtcars)
    narrowed <- ps[, "mpg"]
    # `cyl` is not in colnames(narrowed) but IS in effective_schema
    expect_no_error(subset(narrowed, cyl == 4))
})


# post-join composition: y-side cols visible via effective_schema ####
# `[, j]` and `subset()` consult `.pstore_effective_schema` which walks
# @ops and pulls y-side cols from queued join ops. Lets users narrow to,
# filter on, and introspect cols that came from y after a join.

.make_xy_for_join_compose <- function() {
    x <- parquetStore() |> storeWrite(data.frame(
        cell_ID = letters[1:8],
        gene_count = c(50, 200, 30, 180, 90, 220, 10, 150)
    ))
    y <- parquetStore() |> storeWrite(data.frame(
        cell_ID = letters[1:8],
        qc_pass = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE),
        n_reads = c(100, 500, 50, 400, 250, 800, 20, 350)
    ))
    list(
        x = x, y = y,
        joined = subset(x, gene_count > 50)[
            subset(y, qc_pass), on = "cell_ID", nomatch = NULL
        ]
    )
}

test_that("colnames() reflects post-join effective schema", {
    s <- .make_xy_for_join_compose()
    expect_setequal(colnames(s$joined),
        c("cell_ID", "gene_count", "qc_pass", "n_reads"))
})

test_that("post-join NSE subset on y-side col", {
    s <- .make_xy_for_join_compose()
    res <- s$joined |> subset(n_reads > 300) |> storeRead(output = "tibble")
    expect_setequal(res$cell_ID, c("b", "d", "f", "h"))
    expect_true(all(res$n_reads > 300))
})

test_that("[-narrow to a y-side col after join", {
    s <- .make_xy_for_join_compose()
    res <- s$joined[, "n_reads"] |> storeRead(output = "tibble")
    expect_setequal(names(res), "n_reads")
    expect_equal(nrow(res), 5L)  # inner intersection: b, d, e, f, h
})

test_that("[-narrow to mixed x-side + y-side cols after join", {
    s <- .make_xy_for_join_compose()
    res <- s$joined[, c("cell_ID", "n_reads")] |>
        storeRead(output = "tibble")
    expect_setequal(names(res), c("cell_ID", "n_reads"))
})

test_that("chained joins: cols from a 2nd-level y are visible to subset+[", {
    s <- .make_xy_for_join_compose()
    z <- parquetStore() |> storeWrite(data.frame(
        cell_ID = c("b", "f"), label = c("hi", "lo"),
        stringsAsFactors = FALSE
    ))
    chained <- s$joined[z, on = "cell_ID", nomatch = NULL]
    expect_true("label" %in% colnames(chained))
    # subset on 2nd-level y col
    res <- chained |> subset(label == "hi") |>
        storeRead(output = "tibble")
    expect_equal(res$cell_ID, "b")
    # [-narrow to 2nd-level y col
    res2 <- chained[, "label"] |> storeRead(output = "tibble")
    expect_setequal(names(res2), "label")
})

test_that("non-joined store: colnames unchanged (no regression)", {
    ps <- parquetStore() |> storeWrite(mtcars)
    expect_setequal(colnames(ps), colnames(mtcars))
})


# join op persistence + frozen-snapshot invariants ####
# These lock the contract that the coordinator depends on: a queued join is
# a value-copy snapshot of its inner store at queue time. Surviving saveRDS
# is the persistence pathway after qs's deprecation, and immunity to later
# mutations of the original y is what lets "[-as-view" be safe to assemble
# at getter-call time.

test_that("[-join: chained join saveRDS roundtrip preserves ops + result", {
    x <- parquetStore() |> storeWrite(data.frame(
        cell_ID = letters[1:6], v = 1:6
    ))
    y <- parquetStore() |> storeWrite(data.frame(
        cell_ID = letters[1:6], w = c(10, 20, 30, 40, 50, 60)
    ))
    z <- parquetStore() |> storeWrite(data.frame(
        cell_ID = letters[c(2, 4, 6)], t = c("a", "b", "c")
    ))
    chained <- subset(x, v > 1)[y, on = "cell_ID", nomatch = NULL][
        z, on = "cell_ID", nomatch = NULL
    ]
    expect_length(chained@ops, 3L)
    expect_equal(
        vapply(chained@ops, function(o) o$type, character(1L)),
        c("filter", "join", "join")
    )
    pre <- storeRead(chained, output = "tibble")

    tmp <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp), add = TRUE)
    saveRDS(chained, tmp)
    rt <- readRDS(tmp)
    expect_equal(
        vapply(rt@ops, function(o) o$type, character(1L)),
        c("filter", "join", "join")
    )
    post <- storeRead(rt, output = "tibble")
    expect_equal(pre, post)
})

test_that("[-join: op$y is a frozen snapshot -- later mutations to y don't leak", {
    # If anyone refactors parquetStore to reference semantics, this catches
    # the regression: queuing a join must value-copy y, not alias it.
    x <- parquetStore() |> storeWrite(data.frame(id = letters[1:5], v = 1:5))
    y <- parquetStore() |> storeWrite(data.frame(id = letters[c(1, 3, 5)], w = c(10, 30, 50)))
    j <- x[y, on = "id"]
    expect_length(j@ops[[1L]]$y@ops, 0L)

    # Mutate the original y -- queue more ops on the source store
    y_mut <- subset(y, w > 20)
    expect_length(y_mut@ops, 1L)
    # The queued join's op$y must NOT pick up the mutation
    expect_length(j@ops[[1L]]$y@ops, 0L)
    # Re-running storeRead reflects the original y, not the mutated version
    res <- storeRead(j, output = "tibble")
    expect_equal(sort(stats::na.omit(res$w)), c(10, 30, 50))
})

test_that("[-join: multi-key join on c('k1', 'k2')", {
    # `op_referenced_cols` adds both keys to the upstream projection; effective
    # schema walks join op cols regardless of key arity.
    x <- parquetStore() |> storeWrite(data.frame(
        k1 = c("a", "a", "b", "b"),
        k2 = c(1L, 2L, 1L, 2L),
        v = c(10, 20, 30, 40)
    ))
    y <- parquetStore() |> storeWrite(data.frame(
        k1 = c("a", "b"),
        k2 = c(2L, 1L),
        w = c(100, 200)
    ))
    j <- x[y, on = c("k1", "k2"), nomatch = NULL]
    res <- storeRead(j, output = "tibble")
    expect_equal(nrow(res), 2L)
    expect_setequal(paste(res$k1, res$k2), c("a 2", "b 1"))
    expect_setequal(res$w, c(100, 200))
})
