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
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("cSparse", "ScMaSigPro", function(object, value) {
  if (identical(value, "missing")) {
    return(as.data.frame(object@Sparse@colData)) # Getter
  } else {
    object@Sparse@colData <- DataFrame(value) # Setter
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
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setReplaceMethod("cSparse", "ScMaSigPro", function(object, value) {
  object@Sparse@colData <- DataFrame(value)
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
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("cDense", "ScMaSigPro", function(object, value) {
  if (identical(value, "missing")) {
    return(as.data.frame(object@Dense@colData))
  } else {
    object@Dense@colData <- DataFrame(value)
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
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setReplaceMethod("cDense", "ScMaSigPro", function(object, value) {
  object@Dense@colData <- DataFrame(value)
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
    assayNames <- names(object@Sparse@assays@data@listData)
    if (length(assayNames) > 0) {
      return(as.matrix(object@Sparse@assays@data@listData[[assayNames[1]]]))
    } else {
      stop("No assays found in the object.")
    }
  } else {
    # Return the specific assay's data matrix
    if (!i %in% names(object@Sparse@assays@data@listData)) {
      stop(paste0("Assay '", i, "' not found in the object."))
    }
    return(as.matrix(object@Sparse@assays@data@listData[[i]]))
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
#' @seealso \code{\link{eSparse}} for retrieving expression counts from the Sparse slot.
#'
#' @export
setMethod("eSparse<-", signature(object = "ScMaSigPro", i = "character", value = "matrix"), function(object, i, value) {
  if (!inherits(value, "dgCMatrix")) {
    value <- as(value, "dgCMatrix")
  }

  # Check if the specified assay exists; if not, create a new entry
  if (!i %in% names(object@Sparse@assays@data@listData)) {
    message(paste0("Creating a new assay as '", i, "'."))
    # Update the assay data
    object@Sparse@assays@data@listData[[i]] <- value
  } else if (i %in% names(object@Sparse@assays@data@listData)) {
    warning(paste0("Overwritting assay '", i, "'."))
    # Update the assay data
    object@Sparse@assays@data@listData[[i]] <- value
  }

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
    assayNames <- names(object@Dense@assays@data@listData)
    if (length(assayNames) > 0) {
      return(as.matrix(object@Dense@assays@data@listData[[assayNames[1]]]))
    } else {
      stop("No assays found in the object.")
    }
  } else {
    # Return the specific assay's data matrix
    if (!i %in% names(object@Dense@assays@data@listData)) {
      stop(paste0("Assay '", i, "' not found in the object."))
    }
    return(as.matrix(object@Dense@assays@data@listData[[i]]))
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
#' @seealso \code{\link{eDense}} for retrieving expression data from the Dense slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("eDense<-", signature(object = "ScMaSigPro", i = "character", value = "matrix"), function(object, i, value) {
  if (!inherits(value, "dgCMatrix")) {
    value <- as(value, "dgCMatrix")
  }

  # Check if the specified assay exists; if not, create a new entry
  if (!i %in% names(object@Dense@assays@data@listData)) {
    message(paste0("Creating a new assay as '", i, "'."))
    # Update the assay data
    object@Dense@assays@data@listData[[i]] <- value
  } else if (i %in% names(object@Dense@assays@data@listData)) {
    warning(paste0("Overwritting assay '", i, "'."))
    # Update the assay data
    object@Dense@assays@data@listData[[i]] <- value
  }

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
###############################################################################
# Constructor
ScMaSigPro <- function(Sparse = new("SingleCellExperiment"),
                       Profile = new("VariableProfiles"),
                       Estimate = new("Estimates"),
                       Dense = new("SingleCellExperiment"),
                       Design = new("MatrixDesign"),
                       Parameters = new("ParameterConfig"),
                       Significant = new("Significant")) {
  new("ScMaSigPro",
    Sparse = Sparse,
    Profile = Profile,
    Estimate = Estimate,
    Dense = Dense,
    Design = Design,
    ParameterConfig = Parameters,
    Significant = Significant # new("Significant"),
  )
}
