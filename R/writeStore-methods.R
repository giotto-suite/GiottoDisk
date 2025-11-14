

setMethod("writeStore", "parquetStore", function(x, data, path, ...) {
    checkmate::assert_character(path, len = 1L)
    schema <- eval(x@schema)
    arrow::write_dataset(data,
        path = path,
        format = "parquet",
        schema = schema,
        ...
    )
    invisible(TRUE)
})


setMethod("writeStore", "tileBatch", function(x, data, ...) {
    path <- .tilepath(x@file, x@idx)
    .write_parquet(x, data, path, ...)
    invisible(TRUE)
})

# helpers ####

.write_parquet <- \(x, data, path, ...) {
    arrow::write_dataset(
        dataset = data,
        path = path,
        format = "parquet",
        schema = eval(x@schema),
        ...
    )
}

.tilepath <- \(dir, idx, ...) {
    file.path(dir, sprintf("03d%", as.integer(idx)), "tile")
}

