#' @include pkg_imports.R
NULL

# Custom chaining and lazy eval permits delay of spatial operations until
# ready to collect, even for operations that would normally immediately force
# eval such as sf::st_buffer(). Having this level of control, allows proper
# handling of how the materialization is performed, ie through chunking or
# tiled operations


setMethod("%>%", c("parquetStore", "function"), \(lhs, rhs) {
    op <- new("spatialOperation",
        fn = y,
        params = list(),
        name = deparse(substitute(y))
    )

    new("computeResult",
        type = "lazy",
        operations = list(op),
        source = x
    )
})

setMethod("%>%", c("computeResult", "function"), \(lhs, rhs) {
    op <- new("spatialOperation",
        fn = y,
        params = list(),
        name = deparse(substitute(y))
    )

    x@operations <- c(x@operations, list(op))
    x
})

setMethod("|>", c("parquetStore", "function"), \(lhs, rhs) {
    op <- new("spatialOperation",
        fn = y,
        params = list(),
        name = deparse(substitute(y))
    )

    new("computeResult",
        type = "lazy",
        operations = list(op),
        source = x
    )
})

setMethod("|>", c("computeResult", "function"), \(lhs, rhs) {
    op <- new("spatialOperation",
        fn = y,
        params = list(),
        name = deparse(substitute(y))
    )

    x@operations <- c(x@operations, list(op))
    x
})
