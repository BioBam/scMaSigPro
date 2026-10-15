#' @title Get the data for Sparse Slot.
#'
#' @description
#' Get or set cell level metadata for the Sparse slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot of Sparse. (When Setting)
#'
#' @return `colData` for the `Sparse` Slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom SummarizedExperiment assay assay<- assayNames colData<-
#'
#' @export
setGeneric("cSparse", function(object, value = "missing") standardGeneric("cSparse"))

#' @title Set the data for Sparse Slot
#'
#' @description
#' Set cell level metadata for the Sparse slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot. (When Setting)
#'
#' @return  Modified `ScMaSigPro` object when setting.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("cSparse<-", function(object, value) standardGeneric("cSparse<-"))

#' @title Set or get the Sparse Column Data of an ScMaSigPro Object
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot. (When Setting)
#'
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#'
#' @details
#' When `value` is given, it is checked as by \code{\link{cSparse<-}}.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("cSparse", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(as.data.frame(colData(object@Sparse))) # Getter
    } else {
        object@Sparse <- .set_csparse(object@Sparse, value) # Setter
        return(invisible(object))
    }
})

#' @title Replacement method for cSparse
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot. (When Setting)
#'
#' @return Modified `ScMaSigPro` object.
#'
#' @details
#' `value` must have one row per cell, named after the cells of the `Sparse`
#' slot and in their order, as [cSparse()] returns it; a table without
#' rownames is taken in the order of the cells. Anything else is an error:
#' the cells are the columns of the count matrix, and [sc.squeeze()] finds
#' them by name.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setReplaceMethod("cSparse", "ScMaSigPro", function(object, value) {
    object@Sparse <- .set_csparse(object@Sparse, value)
    return(object)
})
###############################################################################
#' @title Get the data for Dense Slot.
#'
#' @description
#' Get or set cell level metadata for the Dense slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot of Dense. (When Setting)
#'
#' @return `colData` for the `Dense` Slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("cDense", function(object, value = "missing") standardGeneric("cDense"))

#' @title Set the data for Dense Slot
#'
#' @description
#' Set cell level metadata for the Dense slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot. (When Setting)
#'
#' @return  Modified `ScMaSigPro` object when setting.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("cDense<-", function(object, value) standardGeneric("cDense<-"))

#' @title Set or get the Dense Column Data of an ScMaSigPro Object
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot. (When Setting)
#'
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#'
#' @details
#' When `value` is given, it is set as by \code{\link{cDense<-}}.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("cDense", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(as.data.frame(colData(object@Dense)))
    } else {
        object@Dense <- .set_cdense(object@Dense, value)
        return(invisible(object))
    }
})

#' @title Replacement method for cDense
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The new value for the `colData` slot. (When Setting)
#'
#' @return Modified `ScMaSigPro` object.
#'
#' @details
#' The bins are the columns of the pseudo-bulk counts. A table with one row
#' per bin replaces the bin metadata; rows may be renamed, and [pb_counts()]
#' then recomputes the counts for the new bins. A table with another number
#' of rows, such as bins rebuilt by hand, cannot keep the old counts: the
#' `Dense` slot then holds the new bin metadata without counts, a message
#' says so, and [pb_counts()] must be run before the next step.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setReplaceMethod("cDense", "ScMaSigPro", function(object, value) {
    object@Dense <- .set_cdense(object@Dense, value)
    return(object)
})
###############################################################################
#' @title Get Expression Counts from Sparse Slot
#'
#' @description
#' `eSparse` is a generic function for retrieving expression counts stored in the
#' Sparse slot of an \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param i The name of the assay to retrieve expression counts for.
#' Default is "counts", which will return the default expression matrix.
#' @param ... Additional arguments (not used currently).
#'
#' @return An expression matrix for the specified assay from the Sparse slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eSparse", function(object, i = "counts", ...) standardGeneric("eSparse"))

#' @title Set Expression Counts in Sparse Slot
#'
#' @description
#' `eSparse<-` is a generic function for setting expression counts in the Sparse slot
#' of an \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param i The name of the assay for which to set expression counts.
#' Default is "counts", which targets the default expression counts slot.
#' @param value The new expression matrix to be set for the specified assay.
#'
#' @return The modified \code{\link{ScMaSigPro}} object with the updated
#' expression data in the Sparse slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eSparse<-", function(object, i = "counts", value) standardGeneric("eSparse<-"))

#' @title Retrieve Expression Counts from Sparse Slot
#'
#' @description
#' Retrieves expression counts from the Sparse slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An \code{\link{ScMaSigPro}} object from which to retrieve expression data.
#' @param i The name of the assay to retrieve expression counts for.
#' The default value "counts" is a placeholder that results in the method returning
#' the first assay's data found in the Sparse slot.
#' @param ... Additional arguments (not used currently).
#'
#' @return An expression matrix for the specified (or first, by default) assay
#' from the Sparse slot.
#'
#' @seealso \code{\link{eSparse<-}} for setting expression counts in the Sparse slot.
#'
#' @export
setMethod("eSparse", signature(object = "ScMaSigPro"), function(object, i = "counts", ...) {
    if (i == "counts") {
        if (length(assayNames(object@Sparse)) > 0) {
            return(as.matrix(assay(object@Sparse, 1L, withDimnames = FALSE)))
        } else {
            stop("No assays found in the object.")
        }
    } else {
        # Return the specific assay's data matrix
        if (!i %in% assayNames(object@Sparse)) {
            stop("Assay '", i, "' not found in the object.")
        }
        return(as.matrix(assay(object@Sparse, i, withDimnames = FALSE)))
    }
})

#' @title Set Expression Counts in Sparse Slot
#'
#' @description
#' Sets or updates expression counts in the Sparse slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An \code{\link{ScMaSigPro}} object to be modified.
#' @param i The name of the assay for which to set or update expression counts.
#' Must be provided as a character string.
#' @param value The new expression data to set for the specified assay.
#'
#' @return The modified \code{\link{ScMaSigPro}} object with the updated or new expression data
#' in the Sparse slot.
#'
#' @details
#' `value` must have one row per gene and one column per cell of the
#' `Sparse` slot. Its row and column names, when present, must be the gene
#' and cell names of the slot in the same order; missing names are taken
#' from the slot. It is stored as a `dgCMatrix`.
#'
#' @seealso \code{\link{eSparse}} for retrieving expression counts from the Sparse slot.
#'
#' @export
setMethod("eSparse<-", signature(object = "ScMaSigPro", i = "character", value = "matrix"), function(object, i, value) {
    value <- .check_assay(object@Sparse, value, "Sparse", "cell")

    # Check if the specified assay exists; if not, create a new entry
    if (!i %in% assayNames(object@Sparse)) {
        message("Creating a new assay as '", i, "'.")
    } else {
        warning("Overwriting assay '", i, "'.")
    }
    assay(object@Sparse, i, withDimnames = FALSE) <- value

    # Return the modified object
    return(object)
})

##############################################################################
#' @title Retrieve Expression Counts from Dense Slot
#'
#' @description
#' `eDense` is a generic function for retrieving expression counts stored in the Dense slot
#' of an \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param i The name of the dataset to retrieve expression counts for.
#' Default is "bulk.counts", which will return the default expression matrix.
#' @param ... Additional arguments (not used currently).
#'
#' @return An expression matrix for the specified dataset from the Dense slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eDense", function(object, i = "bulk.counts", ...) standardGeneric("eDense"))

#' @title Set Expression Counts in Dense Slot
#'
#' @description
#' `eDense<-` is a generic function for updating or setting expression counts in the Dense slot
#' of an \code{\link{ScMaSigPro}} object.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param i The name of the dataset for which to set expression counts.
#' Default is "bulk.counts".
#' @param value The new expression matrix to be set for the specified dataset.
#'
#' @return The modified \code{\link{ScMaSigPro}} object with updated expression
#' data in the Dense slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eDense<-", function(object, i = "bulk.counts", value) standardGeneric("eDense<-"))

#' @title Retrieve Expression Data from Dense Slot
#'
#' @description
#' Retrieves expression data stored in the Dense slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An \code{\link{ScMaSigPro}} object from which to retrieve
#' the expression data.
#' @param i The name of the dataset to retrieve. The default value is "bulk.counts".
#' @param ... Additional arguments (not used currently).
#'
#' @return Returns an expression matrix for the specified dataset from the Dense slot.
#'
#' @seealso \code{\link{eDense<-}} for setting expression data in the Dense slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("eDense", signature(object = "ScMaSigPro"), function(object, i = "bulk.counts", ...) {
    if (i == "bulk.counts") {
        if (length(assayNames(object@Dense)) > 0) {
            return(as.matrix(assay(object@Dense, 1L, withDimnames = FALSE)))
        } else {
            stop("No assays found in the object.")
        }
    } else {
        # Return the specific assay's data matrix
        if (!i %in% assayNames(object@Dense)) {
            stop("Assay '", i, "' not found in the object.")
        }
        return(as.matrix(assay(object@Dense, i, withDimnames = FALSE)))
    }
})

#' @title Set or Update Expression Data in Dense Slot
#'
#' @description
#' `eDense<-` updates or sets new expression data within the Dense slot of an
#' \code{\link{ScMaSigPro}} object.
#'
#' @param object An \code{\link{ScMaSigPro}} object to be modified.
#' @param i The name of the dataset within the Dense slot to set or update.
#' @param value The expression data to set for the specified dataset.
#'
#' @return The modified \code{\link{ScMaSigPro}} object with the updated or new expression data
#' in the Dense slot.
#'
#' @details
#' `value` must have one row per gene and one column per bin of the `Dense`
#' slot. Its row and column names, when present, must be the gene and bin
#' names of the slot in the same order; missing names are taken from the
#' slot. It is stored as a `dgCMatrix`.
#'
#' @seealso \code{\link{eDense}} for retrieving expression data from the Dense slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("eDense<-", signature(object = "ScMaSigPro", i = "character", value = "matrix"), function(object, i, value) {
    value <- .check_assay(object@Dense, value, "Dense", "bin")

    # Check if the specified assay exists; if not, create a new entry
    if (!i %in% assayNames(object@Dense)) {
        message("Creating a new assay as '", i, "'.")
    } else {
        warning("Overwriting assay '", i, "'.")
    }
    assay(object@Dense, i, withDimnames = FALSE) <- value

    # Return the modified object
    return(object)
})

##############################################################################
#' @title Get or set the Branch Assignment Matrix
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the branch Assignment Matrix.
#'
#' @return Get the Assignment Matrix.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("pathAssign", function(object, value = "missing") standardGeneric("pathAssign"))

#' @title Replacement method for pathAssign.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the branch Assignment Matrix.
#'
#' @return Modified `ScMaSigPro` object with new Assignment Matrix.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("pathAssign<-", function(object, value) standardGeneric("pathAssign<-"))

#' @title Set or get the Assignment Matrix.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the branch Assignment Matrix.
#'
#' @return Modified `ScMaSigPro` object with new Assignment Matrix.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("pathAssign", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(object@Design@assignment_matrix) # Getter
    } else {
        object@Design@assignment_matrix <- as.matrix(value) # Setter
        return(invisible(object))
    }
})

#' @title Replacement method for pathAssign
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the branch Assignment Matrix.
#'
#' @return Modified `ScMaSigPro` object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setReplaceMethod("pathAssign", "ScMaSigPro", function(object, value) {
    object@Design@assignment_matrix <- as.matrix(value)
    return(object)
})
##############################################################################
#' @title Get or set the Predictor Matrix
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the Predictor Matrix.
#'
#' @return Get the Predictor Matrix.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("predictors", function(object, value = "missing") standardGeneric("predictors"))

#' @title Replacement method for predictors.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the Predictors Matrix.
#'
#' @return Modified `ScMaSigPro` object with new Predictor Matrix.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("predictors<-", function(object, value) standardGeneric("predictors<-"))

#' @title Set or get the Predictor Matrix.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the Predictors Matrix.
#'
#' @return Modified `ScMaSigPro` object with new Predictor Matrix.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("predictors", "ScMaSigPro", function(object, value) {
    if (identical(value, "missing")) {
        return(object@Design@predictor_matrix) # Getter
    } else {
        object@Design@predictor_matrix <- as.matrix(value) # Setter
        return(invisible(object))
    }
})

#' @title Replacement method for predictors
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value of the Predictors Matrix.
#'
#' @return Modified `ScMaSigPro` object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setReplaceMethod("predictors", "ScMaSigPro", function(object, value) {
    object@Design@predictor_matrix <- as.matrix(value)
    return(object)
})
##############################################################################
setMethod(
    "show",
    "ScMaSigPro",
    function(object) {
        .ScMaSigPro_show(object)
    }
)
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
