# Versioned layout detection for Xenium/Atera zarr output directories.
#
# Each product archive is detected on its own against the layout specs
# below; parsers dispatch on the returned layout id. When 10x ships a new
# format (or zarr v3 / sharded stores, which Rarr cannot read), detection
# is where the change lands.
#
#   xenium-zarr-v1  Xenium Onboard Analysis v1-v4 (grids/0 transcript
#                   tiles, CSC-by-feature cell_features)
#   atera-zarr-v6   Atera v1 (grid/<i,j> transcript tiles with per-row
#                   cell_id, AnnData CSR cell_feature_matrix)
# The cells archive keeps one layout across both.

.zarr_layout_version <- "xenium-zarr-v1"

# Marker entries per product, per layout, tried in order.
.zarr_layout_specs <- list(
    cells = list(
        "xenium-zarr-v1" = c("cell_id", "cell_summary", "polygon_sets/0")
    ),
    transcripts = list(
        "xenium-zarr-v1" = "grids/0",
        "atera-zarr-v6" = c("grid", "codeword_category")
    ),
    cell_feature_matrix = list(
        "xenium-zarr-v1" = c(
            "cell_features/indptr", "cell_features/indices",
            "cell_features/data"
        ),
        "atera-zarr-v6" = c(
            "X/indptr", "X/indices", "X/data", "obs/barcode",
            "var/feature_name"
        )
    )
)

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

# Layout id of one product archive. Errors on a zarr v3 store or when no
# spec matches.
.zarr_detect_product <- function(path, name) {
    specs <- .zarr_layout_specs[[name]]
    if (is.null(specs)) stop("[zarr] unknown product: ", name, call. = FALSE)
    src <- .zarr_open(path)
    on.exit(.zarr_close(src), add = TRUE)
    if (.zarr_exists(src, "zarr.json")) {
        stop("[detectZarrLayout] ", path, " is a zarr v3 store; only ",
            "zarr v2 is supported. Rarr cannot read zarr v3 -- please ",
            "report the format so support can be added.", call. = FALSE)
    }
    for (id in names(specs)) {
        present <- vapply(specs[[id]], function(m) .zarr_exists(src, m),
            logical(1L))
        if (all(present)) return(id)
    }
    stop("[detectZarrLayout] ", path, " does not match a known ", name,
        " layout (", paste(names(specs), collapse = ", "), "). The export ",
        "may use a newer 10x layout -- please report it.", call. = FALSE)
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
    product_layouts <- list()
    if (isTRUE(validate)) {
        for (name in intersect(names(.zarr_layout_specs),
            names(products)[found])) {
            product_layouts[[name]] <- .zarr_detect_product(
                products[[name]], name
            )
        }
    }
    ids <- unique(unlist(product_layouts))
    gp <- file.path(dir, "gene_panel.json")
    xe <- list.files(dir, pattern = "\\.xenium$", full.names = TRUE)
    list(
        layout_version = if ("atera-zarr-v6" %in% ids) "atera-zarr-v6"
            else .zarr_layout_version,
        products = products,
        product_layouts = product_layouts,
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
#'   `product_layouts` (the layout id detected for each product),
#'   `gene_panel_json`, and `experiment_info`.
#' @examples
#' \dontrun{
#' detectZarrLayout("/path/to/xenium_output")
#' }
#' @export
detectZarrLayout <- function(dir) {
    .zarr_layout_detect(dir, validate = TRUE)
}
