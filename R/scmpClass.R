#' scMaSigProClass
#'
#' A class to represent the ScMaSigPro analysis results and associated data.
#' Inherits from \code{SingleCellExperiment}.
#'
#' @name scMaSigProClass
#' @aliases scMaSigProClass
#' @exportClass scMaSigProClass
#' 
#' @importFrom methods is new
#' 
setClass(
  "scMaSigProClass",
  representation(
    sce = "SingleCellExperiment",
    scPVector = "scPVectorClass",
    scTFit = "scTFitClass",
    compress.sce = "SingleCellExperiment",
    edesign = "edesignClass",
    siggenes = "sigClass",
    addParams = "addParamClass"
  ),
  validity = function(object) {
    # Check sce slot
    if (!validObject(object@sce)) {
      stop("sce slot is not a valid SingleCellExperiment object.")
    }

    # Check scPVectorClass slot
    if (!validObject(object@scPVector)) {
      stop("scPVector slot is not a valid scPVectorClass object.")
    }

    # Check scTFitClass slot
    if (!validObject(object@scTFit)) {
      stop("scTFitClass slot is not a valid scTFitClass object.")
    }

    # Check scTFitClass slot
    if (!validObject(object@compress.sce)) {
      stop("compress.sce slot is not a valid SingleCellExperiment object.")
    }
    if (!validObject(object@edesign)) {
      stop("edesign slot is not a valid edesignClass object.")
    }
    if (!validObject(object@siggenes)) {
      stop("siggenes slot is not a valid sigClass object.")
    }
    if (!validObject(object@addParams)) {
      stop("addParams slot is not a valid addParamClass object.")
    }
  },
  prototype = list(
    scPVector = new("scPVectorClass"), # Assuming you've defined scPVectorClass with its prototype
    scTFit = new("scTFitClass"), # Assuming you've defined scTFitClass with its prototype
    siggenes = new("sigClass"), # Assuming you've defined scTFitClass with its prototype
    addParams = new("addParamClass") # Assuming you've defined scTFitClass with its prototype
  )
)

scMaSigProClass <- function(sce = new("SingleCellExperiment"), # Remove default sce
                            scPVector = new("scPVectorClass"),
                            scTFit = new("scTFitClass"),
                            compress.sce = new("SingleCellExperiment"),
                            edesign = new("edesignClass"),
                            siggenes = new("sigClass"),
                            addParams = new("addParamClass")) {
  new("scMaSigProClass",
    sce = sce,
    scPVector = scPVector,
    scTFit = scTFit,
    compress.sce = compress.sce,
    edesign = edesign,
    siggenes = siggenes,
    addParamClass = addParamClass
  )
}

setMethod(
  "show",
  "scMaSigProClass",
  function(object) {
    .smsp_show(object)
  }
)
