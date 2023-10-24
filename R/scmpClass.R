#' scMaSigProClass
#'
#' A class to represent the ScMaSigPro analysis results and associated data.
#' Inherits from \code{SingleCellExperiment}.
#'
#' @slot sce Object of Class SingleCellExperiment. See \pkg{SingleCellExperiment} for more details.
#' @slot scPVector Object of Class scPVectorClass See \pkg{scPVectorClass} for more details.
#' @slot scTFit Object of Class scTFitClass. See \code{\link{scTFitClass}} for more details.
#' @slot compress.sce
#' @slot edesign Object of Class edesignClass. See \code{\link{edesignClass}} for more details.
#' @slot siggenes
#' @slot addParams Object of Class addParamClass. See \code{\link{addParamClass}} for more details.
#'
#' @name scMaSigProClass
#' @aliases scMaSigProClass-class
#' @rdname scMaSigProClass-class
#' @exportClass scMaSigProClass
#' @importFrom methods is new
#' @keywords classes


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

    # Check compress.sce slot
    if (!validObject(object@compress.sce)) {
      stop("compress.sce slot is not a valid SingleCellExperiment object.")
    }

    # Check edesignClass slot
    if (!validObject(object@edesign)) {
      stop("edesign slot is not a valid edesignClass object.")
    }

    # Check sigClass slot
    if (!validObject(object@siggenes)) {
      stop("siggenes slot is not a valid sigClass object.")
    }

    # Check addParamClass slot
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
    addParamClass = addParams
  )
}

setMethod(
  "show",
  "scMaSigProClass",
  function(object) {
    .smsp_show(object)
  }
)
