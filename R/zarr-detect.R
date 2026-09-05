# Versioned layout detection for Xenium/Atera zarr output directories.
#
# ALL layout assumptions of the converter (grids/0 transcript tiles,
# polygon_sets/{0,1}, CSC-by-feature cell_features, a-p encoded cell ids)
# are keyed on the descriptor this file returns. When 10x ships the final
# Atera format (or zarr v3 / sharded stores, which Rarr cannot read),
# detection is where the change lands — parsers consume the descriptor.

# The layout the current parsers implement (Xenium Onboard Analysis
# v1-v4 zarr; Atera preview data ships the same layout).
.zarr_layout_version <- "xenium-zarr-v1"

# Is `path` a zarr source (zipped archive or unzipped tree)?
.is_zarr_path <- function(path) {
    grepl("\\.zarr(\\.zip)?$", path, ignore.case = TRUE)
}

# Locate one product's zarr in `dir`: prefers `<name>.zarr.zip`, falls
# back to an unzipped `<name>.zarr` directory. NULL when absent.
.zarr_find_product <- function(dir, name) {
    zipped <- file.path(dir, paste0(name, ".zarr.zip"))
    if (file.exists(zipped)) return(zipped)
    unzipped <- file.path(dir, paste0(name, ".zarr"))
    if (dir.exists(unzipped)) return(unzipped)
    NULL
}

# Marker entries that identify the layout inside each product.
.zarr_product_markers <- list(
    cells = c("cell_id", "cell_summary", "polygon_sets/0"),
    transcripts = "grids/0",
    cell_feature_matrix = c(
        "cell_features/indptr", "cell_features/indices",
        "cell_features/data"
    )
)

# Validate one product archive against its markers. A zarr v3 store
# (zarr.json instead of .zarray/.zgroup) gets a targeted error.
.zarr_validate_product <- function(path, name) {
    src <- .zarr_open(path)
    on.exit(.zarr_close(src), add = TRUE)
    if (.zarr_exists(src, "zarr.json")) {
        stop("[detectZarrLayout] ", path, " is a zarr v3 store; only ",
            "zarr v2 (", .zarr_layout_version, ") is supported. ",
            "Rarr cannot read zarr v3 — please report the format so ",
            "support can be added.", call. = FALSE)
    }
    markers <- .zarr_product_markers[[name]]
    present <- vapply(markers, function(m) .zarr_exists(src, m), logical(1L))
    if (!all(present)) {
        stop("[detectZarrLayout] ", path, " does not match the expected ",
            .zarr_layout_version, " layout: missing ",
            paste(markers[!present], collapse = ", "),
            ". The export may use a newer 10x layout — please report it.",
            call. = FALSE)
    }
    invisible(NULL)
}

.zarr_layout_detect <- function(dir, validate = TRUE) {
    checkmate::assert_directory_exists(dir)
    products <- list(
        cells = .zarr_find_product(dir, "cells"),
        transcripts = .zarr_find_product(dir, "transcripts"),
        cell_feature_matrix = .zarr_find_product(dir, "cell_feature_matrix"),
        analysis = .zarr_find_product(dir, "analysis")
    )
    found <- !vapply(products, is.null, logical(1L))
    if (!any(found)) {
        return(list(
            layout_version = NA_character_,
            products = products,
            gene_panel_json = NULL,
            experiment_info = NULL
        ))
    }
    if (isTRUE(validate)) {
        for (name in intersect(names(.zarr_product_markers),
            names(products)[found])) {
            .zarr_validate_product(products[[name]], name)
        }
    }
    gp <- file.path(dir, "gene_panel.json")
    xe <- list.files(dir, pattern = "\\.xenium$", full.names = TRUE)
    list(
        layout_version = .zarr_layout_version,
        products = products,
        gene_panel_json = if (file.exists(gp)) gp else NULL,
        experiment_info = if (length(xe)) xe[[1L]] else NULL
    )
}

#' @name detectZarrLayout
#' @title Detect the zarr layout of a Xenium/Atera output directory
#' @description
#' Inspect an output directory for Xenium/Atera `.zarr.zip` archives (or
#' unzipped `.zarr` trees) and return a layout descriptor naming the
#' detected layout version and per-product paths. Unsupported layouts
#' (zarr v3 stores, or archives missing the expected internal structure)
#' error with a message naming what was found.
#' @param dir Xenium/Atera output directory
#' @returns list with `layout_version` (`NA` when the directory holds no
#'   zarr products), `products` (named paths for `cells`, `transcripts`,
#'   `cell_feature_matrix`, `analysis`; `NULL` when absent),
#'   `gene_panel_json`, and `experiment_info`.
#' @examples
#' \dontrun{
#' detectZarrLayout("/path/to/xenium_output")
#' }
#' @export
detectZarrLayout <- function(dir) {
    .zarr_layout_detect(dir, validate = TRUE)
}
