#' scMaSigProClass
#'
#' A class to represent the ScMaSigPro analysis results and associated data.
#' Inherits from \code{SingleCellExperiment}.
#'
#' @name scMaSigProClass
#' @aliases scMaSigProClass
#' @exportClass scMaSigProClass
setClass(
    "scMaSigProClass",
    representation(
        SingleCellExperiment = "SingleCellExperiment",
        scPVectorClass = "scPVectorClass",
        scTFitClass = "scTFitClass"
    ),
    validity = function(object) {
        # Check SingleCellExperiment slot
        if (!validObject(object@SingleCellExperiment)) {
            stop("SingleCellExperiment slot is not a valid SingleCellExperiment object.")
        }
        
        # Check scPVectorClass slot
        if (!validObject(object@scPVectorClass)) {
            stop("scPVectorClass slot is not a valid scPVectorClass object.")
        }
        
        # Check scTFitClass slot
        if (!validObject(object@scTFitClass)) {
            stop("scTFitClass slot is not a valid scTFitClass object.")
        }
    },
    prototype = list(
        scPVectorClass = new("scPVectorClass"),  # Assuming you've defined scPVectorClass with its prototype
        scTFitClass = new("scTFitClass")  # Assuming you've defined scTFitClass with its prototype
    )
)

scMaSigProClass <- function(SingleCellExperiment,  # Remove default SingleCellExperiment
                            scPVectorClass = new("scPVectorClass"), 
                            scTFitClass = new("scTFitClass")) {
    new("scMaSigProClass", 
        SingleCellExperiment = SingleCellExperiment, 
        scPVectorClass = scPVectorClass, 
        scTFitClass = scTFitClass)
}

setMethod(
    "show",
    "scMaSigProClass",
    function(object) {
        .smsp_show(object)
    }
)