# get number of rows in data
.dplyr_nrow <- function(data) {
    data %>%
        dplyr::summarize(n = n()) %>%
        dplyr::pull(n, as_vector = TRUE)
}

# determine if a function has no content
.is_empty_fun <- function(f) {
    body_code <- deparse(body(f))
    return(identical(body_code, "NULL") ||
               identical(body_code, "{\n    NULL\n}"))
}

#' @name .dt_set_row_index
#' @title data.table Utility: Add row index
#' @description
#' Internal that adds a sequential row index column to a data.frame using
#' data.table's by-reference modification. This modifies the input object
#' in-place for memory efficiency with large in-memory datasets.
#'
#' This is the in-memory counterpart to [.arrow_add_row_index()] for use with
#' materialized data.frames rather than lazy arrow queries.
#' @param x `data.frame` or `data.table`. Will be converted to data.table if
#'   not already.
#' @param offset `integer`. Value to add to the row index. Default is `1`,
#'   meaning the first row will be assigned index `1`. Use higher offsets when
#'   processing data chunks that continue from previous batches.
#' @param col `character`. Name of the column to create for the row index.
#'   Default is `"row_index"`.
#' @return The modified data.table with the new index column added. Note that
#'   the input `x` is also modified by reference.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Add row_index starting from 1
#' df <- data.frame(x = 1:100, y = letters[1:100])
#' .dt_set_row_index(df)
#'
#' # Add row_index starting from 1000 (for tile/chunk continuation)
#' .dt_set_row_index(df, offset = 1000)
#'
#' # Custom column name
#' .dt_set_row_index(df, col = "id")
#' }
.dt_set_row_index <- function(x, offset = 1, col = "row_index") {
    checkmate::assert_data_frame(x)
    x <- data.table::setDT(x)
    x[, (col) := .I + offset]
    x
}

# make a unique identifier across nodes/process
# randomness is enforced despite seed setting using withr-style
# seed from the time and a worker-specific count
.make_uid <- function(n = 8L, include_node = NULL, include_pid = NULL) {
    include_node <- include_node %||% 
        getOption("giottodisk.uid_include_node", FALSE)
    include_pid <- include_pid %||% 
        getOption("giottodisk.uid_include_pid", TRUE)
    count <- getOption("giottodisk.uid_count", 1L)
    on.exit({
        options("giottodisk.uid_count" = count + 1L)
    })
    GiottoUtils::gwith_seed(seed = Sys.time() + count, {
        sampleset <- sample(c(letters, LETTERS, as.character(0:9)), 
            size = n,
            replace = TRUE
        )
        rand <- paste(sampleset, collapse = "")
    })

    parts <- rand
    if (include_pid) {
        pid <- Sys.getpid()
        parts <- c(pid, parts)
    }
    if (include_node) {
        host <- substr(Sys.info()[["nodename"]], 1, 8)
        parts <- c(host, parts)
    }
    paste0(parts, collapse = "_")
}

# generate a character timestamp
.timestamp <- function() {
    as.character(format(Sys.time()))
}

# generate hash
.hash <- function(x) {
    digest::digest(x)
}
