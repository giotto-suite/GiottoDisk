# docs ####


# definitions ####

setGeneric("storeRead", function(store, ...) standardGeneric("storeRead"))
setGeneric("storeWrite", function(store, data, ...) standardGeneric("storeWrite"))
setGeneric("sourceWrite", function(src, data, ...) standardGeneric("sourceWrite"))
setGeneric("sourcePrune", function(src, ...) standardGeneric("sourcePrune"))
