# Technology-agnostic pieces shared by the disk readers.
#
# Readers of binned-grid platforms (VisiumHD, Stereo-seq) relate one unit to
# the next: a 2 um bin sits in one 8 um bin and at most one cell. That
# hierarchy is recorded here from a plain table, so no reader carries its own
# copy.



## unit hierarchy ####

# Record which parent unit each unit belongs to, as cell metadata.
#
# `mapping` is a plain table: an `id_col` naming units of `spat_unit`, plus
# one column per candidate parent unit, named by that unit's spat_unit
# (e.g. `bin008`, `cell`). A parent column is kept only where the relation is
# a function -- every unit maps to at most one parent -- so a column always
# answers "which X is this in". Multivalued pairs (an 8 um bin spanning two
# cells) are dropped with a message. Logical flag columns are carried as-is.
#
# Parent columns are ordinary metadata: grouped statistics, subsetting and
# plotting use them with no dedicated verb.
.add_parent_units <- function(gobject, spat_unit, mapping,
                              id_col = "cell_ID", feat_type = NULL,
                              verbose = NULL) {
    mapping <- data.table::as.data.table(mapping)
    checkmate::assert_subset(id_col, names(mapping))
    ids <- GiottoClass::spatIDs(gobject, spat_unit = spat_unit)
    mapping <- mapping[mapping[[id_col]] %in% ids]

    cand <- setdiff(names(mapping), id_col)
    one_row_each <- !anyDuplicated(mapping[[id_col]])
    keep <- vapply(cand, function(p) {
        if (one_row_each) return(TRUE)
        v <- mapping[[p]]
        if (is.logical(v)) return(TRUE)
        pairs <- unique(mapping[!is.na(v), c(id_col, p), with = FALSE])
        !anyDuplicated(pairs[[id_col]])
    }, logical(1L))
    if (any(!keep)) {
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "[add_parent_units] %s: not one-to-one with %s; skipped",
            toString(cand[!keep]), spat_unit))
    }
    cols <- cand[keep]
    if (!length(cols)) return(gobject)

    if (one_row_each) {
        agg <- mapping[, c(id_col, cols), with = FALSE]
        data.table::setnames(agg, id_col, "cell_ID")
        return(GiottoClass::addCellMetadata(gobject,
            spat_unit = spat_unit, feat_type = feat_type,
            new_metadata = agg, by_column = TRUE, column_cell_ID = "cell_ID"))
    }

    # one row per unit; for flags, TRUE if any member row is TRUE
    flag_cols <- cols[vapply(cols, function(p) is.logical(mapping[[p]]),
                             logical(1L))]
    id_cols <- setdiff(cols, flag_cols)
    agg <- mapping[, c(
        lapply(.SD[, id_cols, with = FALSE], function(v) {
            v <- v[!is.na(v)]
            if (length(v)) v[[1L]] else NA_character_
        }),
        lapply(.SD[, flag_cols, with = FALSE], any)
    ), by = id_col, .SDcols = cols]
    data.table::setnames(agg, id_col, "cell_ID")

    GiottoClass::addCellMetadata(
        gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        new_metadata = agg,
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )
}
