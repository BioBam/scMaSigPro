#' Get or set the sparse column data of a scmp object
#'
#' This function can be used to either retrieve or set the `colData` 
#' from the `sparse` slot of a `scmp` object.
#'
#' @param object An object of class `scmp`.
#' @param value The new value for the `colData` slot to set. 
#'              This parameter is optional for getting the data.
#' @return Returns the `colData` slot when getting. 
#'         Returns the modified `scmp` object when setting.
#' @export
cSparse <- function(object, value = "missing") {
    standardGeneric("cSparse")
}

#' Replacement method for cSparse
#'
#' This method allows replacing the `colData` of the `sparse` slot
#' of a `scmp` object.
#'
#' @param object An object of class `scmp`.
#' @param value The new value to set in the `colData` slot.
#' @return Returns the modified `scmp` object.
#' @export
`cSparse<-` <- function(object, value) {
    standardGeneric("cSparse<-")
}

##############################################################################
setMethod(
    "show",
    "scmp",
    function(object) {
        .scmp_show(object)
    }
)
#' Set the Sparse Column Data of an scmp Object
#' 
#' @param object An object of class `scmp`.
#' @param value The new value for the `colData` slot to set.
#' @return Returns the `colData` slot when getting, and the modified `scmp` object when setting.
#' @export
setMethod("cSparse", "scmp", function(object, value) {
    if (identical(value, "missing")) {
        return(object@sparse@colData)  # Getter
    } else {
        object@sparse@colData <- value  # Setter
        return(invisible(object))
    }
})

#' Set the Sparse Column Data of an scmp Object
#'
#' This method allows setting the `colData` of the `sparse` slot of a `scmp` object.
#'
#' @param object An object of class `scmp`.
#' @param value The new value to set in the `colData` slot.
#' @return Returns the modified `scmp` object.
#' @export
setReplaceMethod("cSparse", "scmp", function(object, value) {
    object@sparse@colData <- value
    return(object)
})



###############################################################################
scmp <- function(sparse = new("SingleCellExperiment"), 
                 profile = new("sigProfileClass"),
                 estimate = new("estimateClass"),
                 dense = new("SingleCellExperiment"),
                 design = new("designClass"),
                 param = new("paramClass"),
                 sig.genes = new("sigClass")) {
    new("scmp",
        sparse = sparse,
        profile = profile,
        estimate = estimate,
        dense = dense,
        design = design,
        paramClass = param,
        sig.genes = sig.genes # new("sigClass"),
    )
}

