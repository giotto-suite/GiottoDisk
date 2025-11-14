

setMethod("process", "dataStore", \(x, FUN, lazy = TRUE, ntiles = NULL, ...) {

})

setMethod("process", "nativeBatch", \(x, FUN, ...) {

})

setMethod("process", "selectBatch", \(x, FUN, ...) {

})

setMethod("process", "tileBatch", \(x, FUN, ...) {
    b <- readStore(x) # batch
    a <- list(b, ...) # argslist
    res <- do.call(FUN, a)

    writeStore(x, res)
})
