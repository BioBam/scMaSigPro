#' Get or set the sparse column data of a ScMaSigPro object
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot to set. Optional for getting.
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#' @export
setGeneric("cSparse", function(object, value = "missing") standardGeneric("cSparse"))

#' Replacement method for cSparse
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot.
#' @return Modified `ScMaSigPro` object.
#' @export
setGeneric("cSparse<-", function(object, value) standardGeneric("cSparse<-"))

#' Set or get the Sparse Column Data of an ScMaSigPro Object
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot. Optional for getting.
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#' @export
setMethod("cSparse", "ScMaSigPro", function(object, value) {
  if (identical(value, "missing")) {
    return(as.data.frame(object@sparse@colData)) # Getter
  } else {
    object@sparse@colData <- DataFrame(value) # Setter
    return(invisible(object))
  }
})

#' Replacement method for cSparse
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot.
#' @return Modified `ScMaSigPro` object.
#' @export
setReplaceMethod("cSparse", "ScMaSigPro", function(object, value) {
  object@sparse@colData <- DataFrame(value)
  return(object)
})
###############################################################################
#' Get or set the sparse column data of a ScMaSigPro object
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot to set. Optional for getting.
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#' @export
setGeneric("cDense", function(object, value = "missing") standardGeneric("cDense"))

#' Replacement method for cSparse
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot.
#' @return Modified `ScMaSigPro` object.
#' @export
setGeneric("cDense<-", function(object, value) standardGeneric("cDense<-"))

#' Set or get the Sparse Column Data of an ScMaSigPro Object
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot. Optional for getting.
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#' @export
setMethod("cDense", "ScMaSigPro", function(object, value) {
  if (identical(value, "missing")) {
    return(as.data.frame(object@dense@colData)) # Getter
  } else {
    object@dense@colData <- DataFrame(value) # Setter
    return(invisible(object))
  }
})

#' Replacement method for cSparse
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot.
#' @return Modified `ScMaSigPro` object.
#' @export
setReplaceMethod("cDense", "ScMaSigPro", function(object, value) {
  object@dense@colData <- DataFrame(value)
  return(object)
})
###############################################################################
#' Set eSparse value for an object
#'
#' @description
#' `eSparse<-` is a generic function for setting eSparse value in an object.
#' Currently, this functionality is not implemented for all object types.
#'
#' @param object The object to be modified.
#' @param value The value to be set for eSparse.
#' @param i object name
#' @return None
#' @export
setGeneric("eSparse<-", function(object, i, value) standardGeneric("eSparse<-"))

#' eSparse value of an object
#'
#' @description
#' `eSparse` is a generic function that can act as both a getter and a setter for the eSparse value of an object.
#' The setter functionality is not implemented for all object types.
#'
#' @param object The object to be accessed or modified.
#' @param value The value to be set for eSparse, if it's a setter.
#' @return The eSparse value of the object if it's a getter.
#' @export
setGeneric("eSparse", function(object, value = "missing") standardGeneric("eSparse"))

#' Get eSparse value for ScMaSigPro objects
#'
#' @description
#' Method to get the eSparse value from ScMaSigPro class objects.
#'
#' @param object The ScMaSigPro object.
#' @param value Dummy parameter, not used.
#' @return The eSparse value from the ScMaSigPro object.
#'
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("eSparse", "ScMaSigPro", function(object, value = "missing") {
  if (identical(value, "missing")) {
    return(assay(object@sparse)) # Getter: Replace 'assay' with the appropriate getter function for your object
  } else {
    return(assay(object@sparse, value)) # Assuming you don't want a setter for this method
  }
})

#' Set eSparse value for ScMaSigPro objects (Not Implemented)
#'
#' @description
#' Method to set the eSparse value for ScMaSigPro class objects.
#' This method is currently not implemented.
#'
#' @param object The ScMaSigPro object.
#' @param value The value to set.
#' @param i object name
#' @return None
#' @export
setMethod("eSparse<-", "ScMaSigPro", function(object, i, value) {
  assay(object@sparse, i) <- value
  return(invisible(object))
})
##############################################################################
#' Set eDense value for an object
#'
#' @description
#' `eDense<-` is a generic function for setting eDense value in an object.
#' Currently, this functionality is not implemented for all object types.
#'
#' @param object The object to be modified.
#' @param value The value to be set for eSparse.
#' @param i object name
#' @return None
#' @export
setGeneric("eDense<-", function(object, i, value) standardGeneric("eDense<-"))

#' eDense value of an object
#'
#' @description
#' `eDense` is a generic function that can act as both a getter and a setter for the eSparse value of an object.
#' The setter functionality is not implemented for all object types.
#'
#' @param object The object to be accessed or modified.
#' @param value The value to be set for eSparse, if it's a setter.
#' @return The eSparse value of the object if it's a getter.
#' @export
setGeneric("eDense", function(object, value = "missing") standardGeneric("eDense"))

#' Get eDense value for ScMaSigPro objects
#'
#' @description
#' Method to get the eDense value from ScMaSigPro class objects.
#'
#' @param object The ScMaSigPro object.
#' @param value Dummy parameter, not used.
#' @return The eDense value from the ScMaSigPro object.
#'
#' @importFrom SummarizedExperiment assay `assay<-`
#' @export
setMethod("eDense", "ScMaSigPro", function(object, value = "missing") {
  if (identical(value, "missing")) {
    return(assay(object@dense)) # Getter: Replace 'assay' with the appropriate getter function for your object
  } else {
    return(assay(object@dense, value)) # Assuming you don't want a setter for this method
  }
})

#' Set eDense value for ScMaSigPro objects (Not Implemented)
#'
#' @description
#' Method to set the eDense value for ScMaSigPro class objects.
#' This method is currently not implemented.
#'
#' @param object The ScMaSigPro object.
#' @param value The value to set.
#' @param i object name
#' @return None
#' @export
setMethod("eDense<-", "ScMaSigPro", function(object, i, value) {
  assay(object@dense, i) <- value
  return(invisible(object))
})
##############################################################################
#' Get or set the branch allocation matrix
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot to set. Optional for getting.
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#' @export
setGeneric("bAlloc", function(object, value = "missing") standardGeneric("bAlloc"))

#' Replacement method for bAlloc
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot.
#' @return Modified `ScMaSigPro` object.
#' @export
setGeneric("bAlloc<-", function(object, value) standardGeneric("bAlloc<-"))

#' Set or get the Sparse Column Data of an ScMaSigPro Object
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot. Optional for getting.
#' @return `colData` when getting, modified `ScMaSigPro` object when setting.
#' @export
setMethod("bAlloc", "ScMaSigPro", function(object, value) {
  if (identical(value, "missing")) {
    return(object@design@assignment_matrix) # Getter
  } else {
    object@design@assignment_matrix <- as.matrix(value) # Setter
    return(invisible(object))
  }
})

#' Replacement method for bAlloc
#' @param object An object of class `ScMaSigPro`.
#' @param value The new value for the `colData` slot.
#' @return Modified `ScMaSigPro` object.
#' @export
setReplaceMethod("bAlloc", "ScMaSigPro", function(object, value) {
  object@design@assignment_matrix <- as.matrix(value)
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
ScMaSigPro <- function(sparse = new("SingleCellExperiment"),
                       profile = new("VariableProfiles"),
                       estimate = new("Estimates"),
                       dense = new("SingleCellExperiment"),
                       design = new("MatrixDesign"),
                       param = new("ParameterConfig"),
                       sig.genes = new("Significant")) {
  new("ScMaSigPro",
    sparse = sparse,
    profile = profile,
    estimate = estimate,
    dense = dense,
    design = design,
    ParameterConfig = param,
    sig.genes = sig.genes # new("Significant"),
  )
}
