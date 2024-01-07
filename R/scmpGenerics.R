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
#' @title Set value for expression counts.
#'
#' @description
#' `eSparse<-` is a generic function for setting expression counts for Sparse
#' in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set.
#' @param i Name of the assay.
#'
#' @return Modified `ScMaSigPro` object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eSparse<-", function(object, i, value) standardGeneric("eSparse<-"))

#' @title Get Value for expression counts.
#'
#' @description
#' `eSparse<-` is a generic function for getting expression counts for Sparse
#' in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set. (If Setting)
#'
#' @return The expression matrix for Sparse slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eSparse", function(object, value = "missing") standardGeneric("eSparse"))

#' @title Get Value for expression counts.
#'
#' @importFrom SummarizedExperiment assay assay<-
#'
#' @description
#' `eSparse` is a generic function for setting/getting expression counts for
#' Sparse in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set. (If Setting)
#'
#' @return The expression matrix for Sparse slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("eSparse", "ScMaSigPro", function(object, value = "missing") {
  if (identical(value, "missing")) {
    return(as.matrix(assay(object@Sparse)))
  } else {
    return(as.matrix(assay(object@Sparse, value)))
  }
})

#' @title Get Value for expression counts.
#'
#' @description
#' `eSparse<-` is a generic function for setting expression counts for
#' Sparse in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set. (If Setting)
#' @param i Name of the Assay.
#'
#' @return Modified `ScMaSigPro` object.
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("eSparse<-", "ScMaSigPro", function(object, i, value) {
  assay(object@Sparse, i) <- value
  return(invisible(object))
})
##############################################################################
#' @title Set value for expression counts.
#'
#' @description
#' `eDense<-` is a generic function for setting expression counts for Dense
#' in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set.
#' @param i Name of the assay.
#'
#' @return Modified `ScMaSigPro` object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eDense<-", function(object, i, value) standardGeneric("eDense<-"))

#' @title Get Value for expression counts.
#'
#' @description
#' `eDense<-` is a generic function for getting expression counts for Dense
#' in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set. (If Setting)
#'
#' @return The expression matrix for Dense slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setGeneric("eDense", function(object, value = "missing") standardGeneric("eDense"))

#' @title Get Value for expression counts.
#'
#' @importFrom SummarizedExperiment assay
#'
#' @description
#' `eDense` is a generic function for setting/getting expression counts for
#' Dense in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set. (If Setting)
#'
#' @return The expression matrix for Dense slot.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("eDense", "ScMaSigPro", function(object, value = "missing") {
  if (identical(value, "missing")) {
    return(assay(object@Dense))
  } else {
    return(assay(object@Dense, value))
  }
})

#' @title Get Value for expression counts.
#'
#' @description
#' `eDense<-` is a generic function for setting expression counts for
#' Sparse in an \code{\link{ScMaSigPro}}.
#'
#' @param object An object of class \code{\link{ScMaSigPro}}.
#' @param value The value to be set. (If Setting)
#' @param i Name of the Assay.
#'
#' @return Modified `ScMaSigPro` object.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @export
setMethod("eDense<-", "ScMaSigPro", function(object, i, value) {
  assay(object@Dense, i) <- value
  return(invisible(object))
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
