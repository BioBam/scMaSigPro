#' Cell metadata of the Sparse slot
#'
#' `cSparse()` returns the cell metadata (`colData`) of the `Sparse` slot of
#' a [ScMaSigPro-class] object; `cSparse<-` replaces it.
#'
#' @param object A [ScMaSigPro-class] object.
#' @param value A data frame with one row per cell, named after the cells of
#'   the `Sparse` slot and in their order, as `cSparse()` returns it. Without
#'   row names, the rows are taken in the order of the cells.
#'
#' @details
#' The cells are the columns of the counts, and [sc.squeeze()] finds them by
#' name, so a `value` with other cells or another order is an error.
#' `cSparse(object, value)` sets `value` as `cSparse(object) <- value` does.
#'
#' @return `cSparse(object)`: a data frame with one row per cell.
#' `cSparse<-` and `cSparse(object, value)`: the updated object, returned
#' invisibly by the latter.
#'
#' @family accessors
#'
#' @importFrom SummarizedExperiment assay assay<- assayNames colData<-
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' cells <- cSparse(scmp.ob)
#' head(cells[, c("Pseudotime", "Path")])
#'
#' # Add a column to the cell metadata
#' cells$early <- cells$Pseudotime < 50
#' cSparse(scmp.ob) <- cells
#'
#' @export
setGeneric("cSparse", function(object, value = "missing") standardGeneric("cSparse"))

#' @rdname cSparse
#' @export
setGeneric("cSparse<-", function(object, value) standardGeneric("cSparse<-"))

#' @rdname cSparse
#' @export
setMethod("cSparse", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(as.data.frame(colData(object@Sparse))) # Getter
    } else {
        object@Sparse <- .set_csparse(object@Sparse, value) # Setter
        validObject(object)
        return(invisible(object))
    }
})

#' @rdname cSparse
#' @export
setReplaceMethod("cSparse", "ScMaSigPro", function(object, value) {
    object@Sparse <- .set_csparse(object@Sparse, value)
    validObject(object)
    return(object)
})
###############################################################################
#' Bin metadata of the Dense slot
#'
#' `cDense()` returns the bin metadata (`colData`) of the `Dense` slot of a
#' [ScMaSigPro-class] object, one row per pseudotime bin; `cDense<-`
#' replaces it.
#'
#' @param object A [ScMaSigPro-class] object.
#' @param value A data frame of bin metadata, one row per bin.
#'
#' @details
#' The bins are the columns of the pseudo-bulk counts. A `value` with one
#' row per bin replaces the metadata and keeps the counts; its rows may be
#' renamed, and [pb_counts()] then recomputes the counts. A `value` with
#' another number of rows, such as bins rebuilt by hand, cannot keep the
#' counts: the `Dense` slot then holds the metadata without counts, a
#' message says so, and [pb_counts()] must be run before the next step.
#' `cDense(object, value)` sets `value` as `cDense(object) <- value` does.
#'
#' @return `cDense(object)`: a data frame with one row per bin. `cDense<-`
#' and `cDense(object, value)`: the updated object, returned invisibly by the
#' latter.
#'
#' @family accessors
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' bins <- cDense(scmp.ob)
#' head(bins[, c("Path", "scmp_binned_pseudotime", "scmp_bin_size")])
#'
#' # Add a column to the bin metadata
#' bins$early <- bins$scmp_binned_pseudotime <= 4
#' cDense(scmp.ob) <- bins
#'
#' @export
setGeneric("cDense", function(object, value = "missing") standardGeneric("cDense"))

#' @rdname cDense
#' @export
setGeneric("cDense<-", function(object, value) standardGeneric("cDense<-"))

#' @rdname cDense
#' @export
setMethod("cDense", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(as.data.frame(colData(object@Dense)))
    } else {
        object@Dense <- .set_cdense(object@Dense, value)
        validObject(object)
        return(invisible(object))
    }
})

#' @rdname cDense
#' @export
setReplaceMethod("cDense", "ScMaSigPro", function(object, value) {
    object@Dense <- .set_cdense(object@Dense, value)
    validObject(object)
    return(object)
})
###############################################################################
#' Cell counts of the Sparse slot
#'
#' `eSparse()` returns an assay of the `Sparse` slot of a [ScMaSigPro-class]
#' object, genes by cells, as a dense matrix; `eSparse<-` sets one.
#'
#' @param object A [ScMaSigPro-class] object.
#' @param i The name of the assay. With the default, `"counts"`, `eSparse()`
#'   returns the first assay when no assay has that name.
#' @param ... Not used.
#' @param value A matrix, a `dgCMatrix` or a data frame of numbers, with one
#'   row per gene and one column per cell of the `Sparse` slot. Its row and
#'   column names, when present, must be the gene and cell names of the slot
#'   in the same order; missing names are filled in.
#'
#' @details
#' `eSparse<-` stores `value` as a `dgCMatrix`. It gives a message when it
#' creates a new assay and a warning when it replaces one.
#'
#' @return `eSparse()`: a numeric matrix. `eSparse<-`: the updated object.
#'
#' @family accessors
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' counts <- eSparse(scmp.ob)
#' counts[1:3, 1:4]
#'
#' # Add a second assay
#' eSparse(scmp.ob, "logcounts") <- log1p(counts)
#'
#' @export
setGeneric("eSparse", function(object, i = "counts", ...) standardGeneric("eSparse"))

#' @rdname eSparse
#' @export
setGeneric("eSparse<-", function(object, i = "counts", value) standardGeneric("eSparse<-"))

#' @rdname eSparse
#' @export
setMethod("eSparse", signature(object = "ScMaSigPro"), function(object, i = "counts", ...) {
    assays <- assayNames(object@Sparse)
    if (i %in% assays) {
        return(as.matrix(assay(object@Sparse, i, withDimnames = FALSE)))
    }
    if (i == "counts") {
        # Without a "counts" assay the default returns the first assay
        if (length(assays) > 0) {
            return(as.matrix(assay(object@Sparse, 1L, withDimnames = FALSE)))
        }
        stop("No assays found in the object.")
    }
    stop("Assay '", i, "' not found in the object.")
})

#' @rdname eSparse
#' @export
setMethod("eSparse<-", signature(object = "ScMaSigPro", i = "character", value = "ANY"), function(object, i, value) {
    value <- .check_assay(object@Sparse, value, "Sparse", "cell")

    # Check if the specified assay exists; if not, create a new entry
    if (!i %in% assayNames(object@Sparse)) {
        message("Creating a new assay as '", i, "'.")
    } else {
        warning("Overwriting assay '", i, "'.")
    }
    assay(object@Sparse, i, withDimnames = FALSE) <- value
    validObject(object)

    # Return the modified object
    return(object)
})

#' @rdname eSparse
#' @export
setMethod(
    "eSparse<-",
    signature(object = "ScMaSigPro", i = "missing", value = "ANY"),
    function(object, i, value) {
        eSparse(object, "counts") <- value
        object
    }
)

##############################################################################
#' Pseudo-bulk counts of the Dense slot
#'
#' `eDense()` returns an assay of the `Dense` slot of a [ScMaSigPro-class]
#' object, genes by bins, as a dense matrix; `eDense<-` sets one.
#'
#' @param object A [ScMaSigPro-class] object.
#' @param i The name of the assay. With the default, `"bulk.counts"`,
#'   `eDense()` returns the first assay when no assay has that name.
#' @param ... Not used.
#' @param value A matrix, a `dgCMatrix` or a data frame of numbers, with one
#'   row per gene and one column per bin of the `Dense` slot. Its row and
#'   column names, when present, must be the gene and bin names of the slot
#'   in the same order; missing names are filled in.
#'
#' @details
#' `eDense<-` stores `value` as a `dgCMatrix`. It gives a message when it
#' creates a new assay and a warning when it replaces one.
#'
#' @return `eDense()`: a numeric matrix. `eDense<-`: the updated object.
#'
#' @family accessors
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' bulk <- eDense(scmp.ob)
#' bulk[1:3, 1:4]
#'
#' # Add a second assay
#' eDense(scmp.ob, "logcounts") <- log1p(bulk)
#'
#' @export
setGeneric("eDense", function(object, i = "bulk.counts", ...) standardGeneric("eDense"))

#' @rdname eDense
#' @export
setGeneric("eDense<-", function(object, i = "bulk.counts", value) standardGeneric("eDense<-"))

#' @rdname eDense
#' @export
setMethod("eDense", signature(object = "ScMaSigPro"), function(object, i = "bulk.counts", ...) {
    assays <- assayNames(object@Dense)
    if (i %in% assays) {
        return(as.matrix(assay(object@Dense, i, withDimnames = FALSE)))
    }
    if (i == "bulk.counts") {
        # Without a "bulk.counts" assay the default returns the first assay
        if (length(assays) > 0) {
            return(as.matrix(assay(object@Dense, 1L, withDimnames = FALSE)))
        }
        stop("No assays found in the object.")
    }
    stop("Assay '", i, "' not found in the object.")
})

#' @rdname eDense
#' @export
setMethod("eDense<-", signature(object = "ScMaSigPro", i = "character", value = "ANY"), function(object, i, value) {
    value <- .check_assay(object@Dense, value, "Dense", "bin")

    # Check if the specified assay exists; if not, create a new entry
    if (!i %in% assayNames(object@Dense)) {
        message("Creating a new assay as '", i, "'.")
    } else {
        warning("Overwriting assay '", i, "'.")
    }
    assay(object@Dense, i, withDimnames = FALSE) <- value
    validObject(object)

    # Return the modified object
    return(object)
})

#' @rdname eDense
#' @export
setMethod(
    "eDense<-",
    signature(object = "ScMaSigPro", i = "missing", value = "ANY"),
    function(object, i, value) {
        eDense(object, "bulk.counts") <- value
        object
    }
)

##############################################################################
#' Path assignment matrix of the design
#'
#' `pathAssign()` returns the path assignment matrix that [sc.set.poly()]
#' stores in the `Design` slot of a [ScMaSigPro-class] object: one row per
#' bin, with the binned pseudotime, a replicate index and a 0/1 column per
#' path. `pathAssign<-` replaces it.
#'
#' @param object A [ScMaSigPro-class] object.
#' @param value A matrix with the rows of the design matrix ([predictors()]),
#'   one per bin and with the same names.
#'
#' @details
#' `pathAssign(object, value)` sets `value` as `pathAssign(object) <- value`
#' does.
#'
#' @return `pathAssign(object)`: a numeric matrix. `pathAssign<-` and
#' `pathAssign(object, value)`: the updated object, returned invisibly by the
#' latter.
#'
#' @family accessors
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' head(pathAssign(scmp.ob))
#'
#' @export
setGeneric("pathAssign", function(object, value = "missing") standardGeneric("pathAssign"))

#' @rdname pathAssign
#' @export
setGeneric("pathAssign<-", function(object, value) standardGeneric("pathAssign<-"))

#' @rdname pathAssign
#' @export
setMethod("pathAssign", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(object@Design@assignment_matrix) # Getter
    } else {
        object@Design@assignment_matrix <- as.matrix(value) # Setter
        validObject(object)
        return(invisible(object))
    }
})

#' @rdname pathAssign
#' @export
setReplaceMethod("pathAssign", "ScMaSigPro", function(object, value) {
    object@Design@assignment_matrix <- as.matrix(value)
    validObject(object)
    return(object)
})
##############################################################################
#' Design matrix of the polynomial model
#'
#' `predictors()` returns the design matrix that [sc.set.poly()] stores in
#' the `Design` slot of a [ScMaSigPro-class] object: one row per bin and one
#' column per term of the full model. `predictors<-` replaces it.
#'
#' @param object A [ScMaSigPro-class] object.
#' @param value A matrix with the rows of the path assignment matrix
#'   ([pathAssign()]), one per bin and with the same names, and one column
#'   per term of the design.
#'
#' @details
#' `predictors(object, value)` sets `value` as `predictors(object) <- value`
#' does.
#'
#' @return `predictors(object)`: a numeric matrix. `predictors<-` and
#' `predictors(object, value)`: the updated object, returned invisibly by the
#' latter.
#'
#' @family accessors
#'
#' @examples
#' data("scmp.ob", package = "scMaSigPro")
#' head(predictors(scmp.ob))
#'
#' @export
setGeneric("predictors", function(object, value = "missing") standardGeneric("predictors"))

#' @rdname predictors
#' @export
setGeneric("predictors<-", function(object, value) standardGeneric("predictors<-"))

#' @rdname predictors
#' @export
setMethod("predictors", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(object@Design@predictor_matrix) # Getter
    } else {
        object@Design@predictor_matrix <- as.matrix(value) # Setter
        validObject(object)
        return(invisible(object))
    }
})

#' @rdname predictors
#' @export
setReplaceMethod("predictors", "ScMaSigPro", function(object, value) {
    object@Design@predictor_matrix <- as.matrix(value)
    validObject(object)
    return(object)
})
##############################################################################
#' Checks behind the setters of the Sparse and Dense slots
#'
#' `.set_csparse()` replaces the cell metadata of the `Sparse` slot. The
#' cells are the columns of the counts, so `value` needs one row per cell,
#' named after the cells in their order; a table without rownames is taken
#' in the order of the cells and keeps their names.
#'
#' `.set_cdense()` replaces the bin metadata of the `Dense` slot. Bins may be
#' renamed, or lose their rownames, before [pb_counts()] recomputes the
#' counts. A table with another number of bins cannot keep the old counts:
#' the slot then holds the bin metadata only, as [sc.squeeze()] builds it
#' before calling [pb_counts()].
#'
#' `.check_assay()` returns `value` ready to be stored as an assay of `sce`:
#' it must have the dimensions of `sce`, and its row and column names, when
#' present, must be those of `sce` in the same order. Missing names are
#' filled in, and the result is a `dgCMatrix`.
#'
#' @param sce The `Sparse` or `Dense` slot.
#' @param value The new table or matrix.
#' @param slot_name `"Sparse"` or `"Dense"`, for the messages.
#' @param unit What a column of `sce` is, `"cell"` or `"bin"`.
#'
#' @return `.set_csparse()` and `.set_cdense()`: the updated
#'   `SingleCellExperiment`. `.check_assay()`: a `dgCMatrix`.
#'
#' @name set_slot
#' @noRd
NULL

.set_csparse <- function(sce, value) {
    value <- DataFrame(value)
    .check(nrow(value) == ncol(sce),
        msg = paste0(
            "'value' must have one row per cell of the 'Sparse' slot: ",
            nrow(value), " rows for ", ncol(sce), " cells."
        )
    )
    if (is.null(rownames(value))) {
        rownames(value) <- colnames(sce)
    }
    .check(identical(rownames(value), colnames(sce)),
        msg = paste(
            "The rownames of 'value' must be the cell names of the 'Sparse'",
            "slot, in the same order."
        )
    )
    colData(sce) <- value
    sce
}

.set_cdense <- function(sce, value) {
    value <- DataFrame(value)
    if (nrow(value) == ncol(sce)) {
        colData(sce) <- value
        return(sce)
    }
    message(
        "'value' has ", nrow(value), " bins instead of ", ncol(sce),
        ": the pseudo-bulk counts of the 'Dense' slot are removed; run ",
        "pb_counts() to compute them for the new bins."
    )
    SingleCellExperiment(colData = value)
}

.check_assay <- function(sce, value, slot_name, unit) {
    if (is.data.frame(value)) {
        value <- as.matrix(value)
    }
    .check(identical(dim(value), dim(sce)),
        msg = paste0(
            "'value' must have one row per gene and one column per ", unit,
            " of the '", slot_name, "' slot: ", nrow(value), " x ",
            ncol(value), " given, ", nrow(sce), " x ", ncol(sce), " needed."
        )
    )
    if (is.null(rownames(value))) {
        rownames(value) <- rownames(sce)
    }
    if (is.null(colnames(value))) {
        colnames(value) <- colnames(sce)
    }
    .check(identical(rownames(value), rownames(sce)),
        msg = paste0(
            "The rownames of 'value' must be the gene names of the '",
            slot_name, "' slot, in the same order."
        )
    )
    .check(identical(colnames(value), colnames(sce)),
        msg = paste0(
            "The colnames of 'value' must be the ", unit, " names of the '",
            slot_name, "' slot, in the same order."
        )
    )
    if (!inherits(value, "dgCMatrix")) {
        value <- as(value, "dgCMatrix")
    }
    value
}
