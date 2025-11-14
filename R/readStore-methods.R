
setMethod("readStore", "parquetStore", function(x, i = NULL, j = NULL, tile = NULL, ...) {
    schema <- eval(x@schema)
    ds <- arrow::open_dataset(x@src, schema = schema)

    if (!is.null(tile)) {
        extent <- ext(do.call(`[`, c(list(x@tiles), tile)))
        ds <- ds |>
            dplyr::filter(
                x >= extent[1], x <= extent[2],
                y >= extent[3], y <= extent[4]
            )
    }

    ds |> dplyr::collect()
})

