# atera-zarr-v6 cell_feature_matrix (AnnData, CSR by cell) as an
# exprInput.
#
# `X` is already cell-major with the 10x h5 `indptr`, so storeRead() cuts
# budget-sized cell ranges the way the h5 input does: one range read of
# `indices` / `data` per range, no reorder pass and no rescans (cf. the
# feature-major xenium-zarr-v1 matrix in tenxZarrInput).

# Bytes held per stored value while a range is decoded and written
# (measured on the full Atera kidney matrix).
.zarr_csr_bytes_per_record <- 224

# AnnData categorical column -> character (codes index categories; -1 is
# missing).
.zarr_categorical <- function(src, prefix) {
    codes <- as.integer(.zarr_array(src, paste0(prefix, "/codes")))
    cats <- as.character(.zarr_array(src, paste0(prefix, "/categories")))
    out <- rep(NA_character_, length(codes))
    ok <- codes >= 0L
    out[ok] <- cats[codes[ok] + 1L]
    out
}

#' @title 10x Atera zarr CSR expression input
#' @name tenxZarrCsrInput-class
#' @description
#' [tenxZarrInput()] for an Atera v1 `cell_feature_matrix.zarr.zip`
#' (AnnData, CSR by cell). Construction reads cell barcodes, the feature
#' catalog and `indptr`; [storeRead()] streams cell windows.
#' @slot indptr numeric. Per-cell offsets into `X/indices` / `X/data`.
#' @keywords internal
setClass("tenxZarrCsrInput",
    contains = "tenxZarrInput",
    slots = list(indptr = "numeric")
)

.tenx_zarr_csr_input <- function(zarr_path, feature_id_col = 2L) {
    src <- .zarr_open(zarr_path)
    on.exit(.zarr_close(src), add = TRUE)
    xattr <- .zarr_attrs(src, "X")
    if (!identical(xattr$`encoding-type`, "csr_matrix")) {
        stop("[tenxZarrInput] expected a CSR (cell-major) matrix in X, ",
            "found '", xattr$`encoding-type`, "'", call. = FALSE)
    }
    indptr <- .zarr_array(src, "X/indptr")
    feat_name <- .zarr_array(src, "var/feature_name")
    feat_id <- if (.zarr_exists(src, "var/feature_id")) {
        .zarr_array(src, "var/feature_id")
    } else {
        feat_name
    }
    new("tenxZarrCsrInput",
        path = normalizePath(zarr_path),
        cell_ids = .zarr_array(src, "obs/barcode"),
        feat_ids = .disambiguate_feat_ids(
            if (feature_id_col == 1L) feat_id else feat_name
        ),
        n_cells = length(indptr) - 1L,
        n_genes = length(feat_name),
        feat_types = .zarr_categorical(src, "var/feature_type"),
        nnz = indptr[length(indptr)],
        indptr = indptr
    )
}

setMethod("storeRead", signature("tenxZarrCsrInput"), function(store, ...) {
    src <- .zarr_open(store@path)
    ip <- store@indptr
    ranges <- .ingest_cell_ranges(diff(ip), .ingest_batch_rows(
        .zarr_csr_bytes_per_record,
        pinned = getOption("giottodisk.h5_batch_rows")))
    k <- 1L
    next_batch <- function() {
        while (k <= length(ranges)) {
            c0 <- ranges[[k]][1L]
            c1 <- ranges[[k]][2L]
            k <<- k + 1L
            lo <- ip[c0] + 1
            hi <- ip[c1 + 1L]
            if (hi < lo) next
            out <- data.table::data.table(
                row_id = rep.int(c0:c1, diff(ip[c0:(c1 + 1L)])),
                col_id = .zarr_array(src, "X/indices",
                    range = c(lo, hi), u4_as_integer = TRUE) + 1L,
                value = as.double(.zarr_array(src, "X/data",
                    range = c(lo, hi)))
            )
            # AnnData does not require sorted indices within a row
            if (any(diff(out$col_id) < 0L & diff(out$row_id) == 0L)) {
                data.table::setorder(out, row_id, col_id)
            }
            return(out)
        }
        NULL
    }
    list(
        next_batch = next_batch,
        close = function() .zarr_close(src),
        cell_ids = function() store@cell_ids,
        feat_ids = function() store@feat_ids,
        n_cells = function() store@n_cells,
        n_genes = function() store@n_genes
    )
})
