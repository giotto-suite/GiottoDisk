
#' @name .arrow_add_row_index
#' @title arrow Utility: Add row index
#' @description
#' Internal that adds a sequential row index column to arrow data processed
#' batch-by-batch. The index counter is maintained across batches using a
#' closure, ensuring continuous sequential numbering even when data is
#' processed in chunks.
#'
#' This is particularly useful when writing large datasets to parquet where
#' row indices need to be assigned without loading the entire dataset into
#' memory.
#' @param data An object coercible to [arrow::RecordBatchReader]. Common inputs
#'   include arrow queries (`arrow_dplyr_query`), `FileSystemDataset`, arrow
#'   `Table`, or `data.frame`/`tibble`.
#' @param col `character`. Name of the column to create for the row index.
#'   Default is `"row_index"`.
#' @param offset `integer`. Starting value for the row index. Default is `0`,
#'   meaning the first row will be assigned index `1`. Use non-zero offsets
#'   when appending to existing indexed data or processing tiles.
#' @return An arrow [RecordBatchReader] with the new index column added.
#'   The indices continue sequentially across all batches.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Add row_index starting from 1
#' ds <- arrow::open_dataset("path/to/data")
#' indexed_ds <- .arrow_add_row_index(ds)
#'
#' # Add row_index starting from 1000 (for appending)
#' indexed_ds <- .arrow_add_row_index(ds, offset = 1000)
#'
#' # Custom column name
#' indexed_ds <- .arrow_add_row_index(ds, col = "id")
#' }
.arrow_add_row_index <- function(data, col = "row_index", offset = 0) {
    counter <- offset
    # function to apply to each batch
    add_idx <- function(batch) {
        n_rows <- dim(batch)[1L]
        idx <- seq(from = counter + 1L, length.out = n_rows)
        # update counter for next batch
        counter <<- counter + n_rows
        batch[[col]] <- idx
        batch
    }

    data <- .arrow_map_batches(data, add_idx, .lazy = TRUE)
}

#' @name .arrow_map_batches
#' @title arrow Utility: Map batches
#' @description
#' Internal that applies a transformation function `FUN` to an arrow dataset
#' batch by batch rather than all at once, returning a new lazy
#' `RecordBatchReader` with the transformed data. This can be used to tell
#' arrow how to perform a data transformation that is not natively supported
#' in a memory-friendly way by performing it in batches. This implementation
#' is based on `arrow::map_batches()`. This may be superceded in arrow's own API
#' in the future. The code runs in R instead of arrow C++ code
#'
#' **Key behaviors:**
#' * Lazy evaluation: Processes batches on-demand (unless `.lazy = FALSE`)
#' * Schema inference: If schema not provided, processes first batch to infer
#'   output schema, then uses that schema for remaining batches
#' * Streaming: Returns a reader that can be consumed incrementally without
#'   loading entire dataset into memory
#' * Preserves batch structure: Each batch is transformed independently
#' @param X An object coercible to [arrow::RecordBatchReader]. Common inputs
#'   include arrow queries (`arrow_dplyr_query`), `FileSystemDataset`, arrow
#'   `Table`, or `data.frame`/`tibble`.
#' @param FUN `function`. A function to apply. The first parameter should
#' accept a batch subset of the full dataset from  `X`.
#' @param \dots additional params to pass to `FUN`
#' @param .schema arrow `Schema` (optional). If not explicitly provided, the
#' schema will be inferred from the first batch.
#' @param .lazy `logical`. If `TRUE` (default), returns lazy reader. If `FALSE`,
#'   materializes all batches immediately.
#' @return An arrow [RecordBatchReader] with transformed batches
#' @examples
#' \dontrun{
#' # Add row indices to batches
#' ds <- arrow::open_dataset("path/to/data")
#' .arrow_map_batches(ds, function(batch) {
#'     batch$row_id <- seq_len(nrow(batch))
#'     batch
#' })
#' }
#' @keywords internal
.arrow_map_batches <- function(X, FUN, ..., .schema = NULL, .lazy = TRUE) {
    reader <- arrow::as_record_batch_reader(X)
    dots <- list(...)
    if (is.null(.schema)) {
        batch <- reader$read_next_batch()
        if (is.null(batch)) {
            stop("Can't infer schema from a RecordBatchReader",
                 " with zero batches\n")
        }
        first_result <- arrow::as_record_batch(
            do.call(FUN, c(list(batch), dots))
        )
        .schema <- first_result$schema
        fun <- function() {
            if (!is.null(first_result)) {
                result <- first_result
                first_result <<- NULL
                result
            }
            else {
                batch <- reader$read_next_batch()
                if (is.null(batch)) {
                    NULL
                }
                else {
                    arrow::as_record_batch(
                        x = do.call(FUN, c(list(batch), dots)),
                        schema = .schema
                    )
                }
            }
        }
    }
    else {
        fun <- function() {
            batch <- reader$read_next_batch()
            if (is.null(batch)) {
                return(NULL)
            }
            arrow::as_record_batch(
                x = do.call(FUN, c(list(batch), dots)),
                schema = .schema
            )
        }
    }
    reader_out <- arrow::as_record_batch_reader(fun, schema = .schema)
    if (!.lazy) {
        reader_out <- arrow::RecordBatchReader$create(
            batches = reader_out$batches(),
            schema = .schema)
    }
    reader_out
}

 .arrow_sample_max_rows <- function(atab, nmax) {                        
      nmax <- as.integer(nmax)                                            
                                                                          
      total <- atab |>
          dplyr::count() |>
          dplyr::collect() |>
          dplyr::pull(n)

      if (total <= nmax) return(atab)

      k <- as.integer(ceiling(total / nmax))
      dplyr::filter(atab, (row_index - 1L) %% k == 0L)
  }
