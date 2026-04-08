# docs ####


# definitions ####

setGeneric("storeRead", function(store, ...) standardGeneric("storeRead"))
setGeneric("storeWrite", function(store, data, ...) standardGeneric("storeWrite"))
setGeneric("storeExists", function(x, ...) standardGeneric("storeExists"))
setGeneric("specialCols", function(store) standardGeneric("specialCols"))
setGeneric("sourceWrite", function(src, data, ...) standardGeneric("sourceWrite"))
setGeneric("sourcePrune", function(src, ...) standardGeneric("sourcePrune"))
setGeneric("snapshotSave", function(src, x, ...) standardGeneric("snapshotSave"))
setGeneric("snapshotLoad", function(src, ...) standardGeneric("snapshotLoad"))
setGeneric("snapshotDelete", function(src, name, ...) standardGeneric("snapshotDelete"))
setGeneric("rowSample", function(x, size, ...) standardGeneric("rowSample"))