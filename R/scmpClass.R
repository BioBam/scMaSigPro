#' scMaSigProClass
#'
#' A class to represent the ScMaSigPro analysis results and associated data.
#' Inherits from \code{SingleCellExperiment}.
#'
#' @slot sce Object of Class SingleCellExperiment. See \pkg{SingleCellExperiment} for more details.
#' @slot scPVector Object of Class scPVectorClass See \pkg{scPVectorClass} for more details.
#' @slot scTFit Object of Class scTFitClass. See \code{\link{scTFitClass}} for more details.
#' @slot compress.sce ABC
#' @slot edesign Object of Class edesignClass. See \code{\link{edesignClass}} for more details.
#' @slot addParams Object of Class addParamClass. See \code{\link{addParamClass}} for more details.
#' @slot sig.genes ABC
#' @slot distribution The distribution function to be used in the glm model.
#'
#' @name scMaSigProClass
#' @aliases scMaSigProClass-class
#' @rdname scMaSigProClass-class
#' @exportClass scMaSigProClass
#' @importFrom methods is new as
#' @keywords classes


setClass(
    "sigClass",
    representation(
        summary = "ANY",
        sig.genes = "list"
    ),
    validity = function(object) {
        if (!is.list(object@sig.genes)) {
            stop("sig.genes slot must be a list")
        }
    },
    prototype = list(
        summary = list(),
        sig.genes = list()
    )
)

setClass(
  "scMaSigProClass",
  representation(
    sce = "SingleCellExperiment",
    scPVector = "scPVectorClass",
    scTFit = "scTFitClass",
    compress.sce = "SingleCellExperiment",
    edesign = "edesignClass",
    addParams = "addParamClass",
    sig.genes = "sigClass",
    distribution = "ANY"
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

    # Check addParamClass slot
    if (!validObject(object@addParams)) {
      stop("addParams slot is not a valid addParamClass object.")
    }
      # Check addParamClass slot
      if (!validObject(object@addParams)) {
          stop("'sig.genes' slot is not a valid addParamClass object.")
      }
  },
  prototype = list(
    scPVector = new("scPVectorClass"), # Assuming you've defined scPVectorClass with its prototype
    scTFit = new("scTFitClass"), # Assuming you've defined scTFitClass with its prototype
    addParams = new("addParamClass"), # Assuming you've defined scTFitClass with its prototype
    sig.genes = new("sigClass"),
    distribution = MASS::negative.binomial(theta = 1)
  )
)

scMaSigProClass <- function(sce = new("SingleCellExperiment"), # Remove default sce
                            scPVector = new("scPVectorClass"),
                            scTFit = new("scTFitClass"),
                            compress.sce = new("SingleCellExperiment"),
                            edesign = new("edesignClass"),
                            addParams = new("addParamClass"),
                            sig.genes = new("sigClass"),
                            distribution = MASS::negative.binomial(theta = 1)) {
  new("scMaSigProClass",
    sce = sce,
    scPVector = scPVector,
    scTFit = scTFit,
    compress.sce = compress.sce,
    edesign = edesign,
    addParamClass = addParams,
    sig.genes = new("sigClass"),
    distribution = distribution
  )
}

setMethod(
  "show",
  "scMaSigProClass",
  function(object) {
    .smsp_show(object)
  }
)
