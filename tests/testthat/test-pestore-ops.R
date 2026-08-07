# Tests for the op-chain edit helpers on parquetExprStore:
#   * .pe_push_op phase routing + the monotonic rule
#   * .pe_demote_ops cascade, selection by index or type, order preservation

.ops_store <- function(ops = list(), post_ops = list()) {
    pe <- parquetExprStore(path = tempfile(fileext = ".parquet"))
    pe@ops <- ops
    pe@post_ops <- post_ops
    pe
}

.op_types <- function(x) vapply(x, function(o) o$type, character(1L))

test_that(".pe_push_op routes by phase and rejects lazy after post", {
    pe <- .ops_store()
    pe <- .pe_push_op(pe, list(type = "a"), phase = "lazy")
    expect_equal(.op_types(pe@ops), "a")
    expect_length(pe@post_ops, 0L)

    pe <- .pe_push_op(pe, list(type = "b"), phase = "post")
    expect_equal(.op_types(pe@post_ops), "b")

    expect_error(.pe_push_op(pe, list(type = "c"), phase = "lazy"),
                 "cannot queue a lazy op after a post op")
})

test_that(".pe_demote_ops moves the selected op and all after it", {
    pe <- .ops_store(ops = list(list(type = "a"), list(type = "b"),
                                list(type = "c")))
    out <- .pe_demote_ops(pe, 2L)

    expect_equal(.op_types(out@ops), "a")
    expect_equal(.op_types(out@post_ops), c("b", "c"))
})

test_that(".pe_demote_ops selects by op type", {
    pe <- .ops_store(ops = list(list(type = "a"), list(type = "norm_libsize"),
                                list(type = "log")))
    out <- .pe_demote_ops(pe, "norm_libsize")

    expect_equal(.op_types(out@ops), "a")
    expect_equal(.op_types(out@post_ops), c("norm_libsize", "log"))
    # payload travels with the record, not just the type tag
    pe@ops[[2]]$scalef <- data.table::data.table(x = 1)
    expect_equal(.pe_demote_ops(pe, "norm_libsize")@post_ops[[1]]$scalef,
                 data.table::data.table(x = 1))
})

test_that(".pe_demote_ops lands the moved block ahead of existing post ops", {
    pe <- .ops_store(ops = list(list(type = "a"), list(type = "b")),
                     post_ops = list(list(type = "r1"), list(type = "r2")))
    out <- .pe_demote_ops(pe, 1L)

    expect_length(out@ops, 0L)
    expect_equal(.op_types(out@post_ops), c("a", "b", "r1", "r2"))
})

test_that(".pe_demote_ops is stackable and reaches a fully demoted chain", {
    pe <- .ops_store(ops = list(list(type = "a"), list(type = "b"),
                                list(type = "c")))
    out <- .pe_demote_ops(.pe_demote_ops(pe, 3L), 1L)

    expect_length(out@ops, 0L)
    expect_equal(.op_types(out@post_ops), c("a", "b", "c"))
    # nothing left to demote
    expect_error(.pe_demote_ops(out, 1L), "must select one of the 0 ops")
})

test_that(".pe_demote_ops rejects a selection that is not on @ops", {
    pe <- .ops_store(ops = list(list(type = "a")),
                     post_ops = list(list(type = "log")))

    expect_error(.pe_demote_ops(pe, "log"), "no op of type 'log' on @ops")
    expect_error(.pe_demote_ops(pe, 2L), "must select one of the 1 ops")
    expect_error(.pe_demote_ops(pe, 0L), "must select one of the 1 ops")
})
