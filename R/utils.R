
# Helper function to convert a SpatExtent object to a simple numeric vector
# and strip the names
.ext_to_num_vec <- function(x) {
    out <- x[]
    names(out) <- NULL
    out
}
