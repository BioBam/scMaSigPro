# scmpClass with all the super-classes
# Description


###############################################################################
#' Class "paramClass"
#'
#' A Class holding information about the parameters used.
#'
#' @slot bin_pseudotime_colname A character string representing the column name for binned pseudotime values.
#' @slot path_prefix A character string representing the prefix for path labeling.
#' @slot root_label A character string representing the label for the tree root.
#' @slot pseudotime_colname A character string representing the column name for pseudotime values.
#' @slot bin_method A character string representing the algorithm used for binning.
#' @slot path_colname A character string representing the column name for path values.
#' @slot bin_colname A character string representing the column name for bin values.
#' @slot bin_size_colname A character string representing the column name for bin sizes.
#' @slot bin_members_colname A character string representing the column name for bin members.
#' @slot annotation_col A character string representing the column name for cell type annotation.
#' @slot g Integer. Number of genes taken in the regression fit.
#' @slot Q Numeric. Significance Level.
#' @slot min.obs Numeric. Minimum value to estimate the model (degree+1) x Groups + 1. (Default = 6).
#' @slot MT.adjust A character string specifying the Pvalue correction method used.
#' @slot epsilon Numeric. Convergence tolerance.
#' @slot step.method A character string specifying the imputed step method for stepwise regression.
#' @slot useWeights A logical value specifying whether to use weights during model fitting.
#' @slot offset A logical value specifying whether to use offset during model fitting.
#' @slot useInverseWeights A logical value specifying whether to take inverse of the weights.
#' @slot logOffset A logical value specifying whether to take the logarithm of the offsets during model fitting.
#' @slot logWeights A logical value specifying whether to take the logarithm of the weights during model fitting.
#' @slot max_it Integer. Maximum number of iterations to fit the model.
#' @slot poly_degree Integer with the polynomial degree to fit the regression. 1
#' specifies a linear regression, 2 a quadratic regression, etc.
#'
#' @name paramClass
#' @aliases paramClass-class
#' @rdname paramClass-class
#' @exportClass paramClass
#' @importFrom methods is new
#' @keywords classes

setClass(
  "paramClass",
  representation(
    bin_pseudotime_colname = "character",
    path_prefix = "character",
    root_label = "character",
    pseudotime_colname = "character",
    step.method = "character",
    bin_method = "character",
    path_colname = "character",
    bin_colname = "character",
    bin_size_colname = "character",
    bin_members_colname = "character",
    annotation_col = "character",
    g = "integer",
    Q = "numeric",
    min.obs = "numeric",
    MT.adjust = "character",
    epsilon = "numeric",
    useWeights = "logical",
    offset = "logical",
    useInverseWeights = "logical",
    logOffset = "logical",
    logWeights = "logical",
    max_it = "integer",
    poly_degree = "integer"
  ),
  validity = function(object) {
    errors <- character(0)

    # Check if any character slots are empty or not of type character
    char_slots <- c(
      "bin_pseudotime_colname", "path_prefix", "root_label",
      "pseudotime_colname", "bin_method",
      "path_colname", "bin_colname", "bin_size_colname",
      "bin_members_colname", "MT.adjust", "step.method",
      "annotation_col"
    )

    for (slot_name in char_slots) {
      slot_value <- slot(object, slot_name) # Corrected line
      if (length(slot_value) == 0 || !is.character(slot_value)) { # Corrected line
        errors <- c(errors, paste(slot_name, "must not be empty and should be of type character."))
      }
    }

    # Check for slot g
    if (!is.integer(object@g)) {
      stop("Slot 'g' must be an integer.")
    }


    # Check for slot Q
    if (!is.numeric(object@Q)) {
      stop("Slot 'Q' must be numeric.")
    }

    if (!is.numeric(object@epsilon)) {
      stop("Slot 'epsilon' must be numeric.")
    }

    # Check for slot min.obs
    if (!is.numeric(object@min.obs)) {
      stop("Slot 'min.obs' must be an integer.")
    }

    # Check if any of the character slots have multiple values
    for (slot_name in char_slots) {
      slot_value <- slot(object, slot_name) # Corrected line
      if (length(slot_value) > 1) { # Corrected line
        errors <- c(errors, paste(slot_name, "should not contain multiple values."))
      }
    }

    if (length(errors) == 0) TRUE else errors
  },
  prototype = list(
    bin_pseudotime_colname = "scmp_binned_pseudotime",
    path_prefix = "Path",
    root_label = "root",
    pseudotime_colname = "Pseudotime",
    path_colname = "Path",
    bin_method = "Sturges",
    bin_colname = "scmp_bin",
    g = 0L,
    Q = 0.05,
    min.obs = 1,
    bin_size_colname = "scmp_bin_size",
    bin_members_colname = "scmp_bin_members",
    MT.adjust = "BH",
    step.method = "backward",
    epsilon = 1e-8,
    useWeights = TRUE,
    offset = TRUE,
    useInverseWeights = FALSE,
    logOffset = FALSE,
    max_it = 100L,
    logWeights = FALSE,
    annotation_col = "cell_type",
    poly_degree = 2L
  )
)

###############################################################################
#' Class "edesignClass"
#'
#' An S4 class to represent an edesign object with associated data.
#' This class contains three slots: dis, groups.vector, and edesign.
#'
#' @slot dis A data frame containing the design matrix of dummies for fitting the GLM.
#' @slot groups.vector A character vector specifying the experimental group to which
#' each variable belongs to.
#' @slot edesign A data frame describing the experimental design. Rows must contain
#' cells and columns experiment descriptors. The matrix must be binarized.
#'
#' @name edesignClass
#' @aliases edesignClass-class
#' @rdname edesignClass-class
#' @exportClass edesignClass
#' @importFrom methods is new
#' @keywords classes
#'
#' @section Validity:
#'   Valid objects must have:
#'   \itemize{
#'     \item{dis}{A valid data frame.}
#'     \item{groups.vector}{A valid character vector.}
#'     \item{edesign}{A valid data frame.}
#'   }
#'

setClass(
    "edesignClass",
    representation(
        dis = "matrix",
        groups.vector = "character",
        edesign = "matrix"
    ),
    prototype = list(
        dis = matrix(NA, nrow = 0, ncol = 0),
        groups.vector = character(),
        edesign = matrix(NA, nrow = 0, ncol = 0)
    ),
    validity = function(object) {
        if (!validObject(object@dis)) {
            stop("dis slot is not a valid data frame.")
        }
        if (!validObject(object@groups.vector)) {
            stop("groups.vector slot is not a valid character vector.")
        }
        if (!validObject(object@edesign)) {
            stop("poly_degree slot is not a valid matrix")
        }
        TRUE
    }
)
###############################################################################


#' scMaSigProClass
#'
#' A class to represent the ScMaSigPro analysis results and associated data.
#' Inherits from \code{SingleCellExperiment}.
#'
#' @slot sparse Object of Class SingleCellExperiment. See \pkg{SingleCellExperiment} for more details.
#' @slot scPVector Object of Class scPVectorClass See \pkg{scPVectorClass} for more details.
#' @slot scTFit Object of Class scTFitClass. See \code{\link{scTFitClass}} for more details.
#' @slot dense  Object of Class SingleCellExperiment. See \pkg{SingleCellExperiment} for more details.
#' @slot edesign Object of Class edesignClass. See \code{\link{edesignClass}} for more details.
#' @slot param Object of Class paramClass. See \code{\link{paramClass}} for more details.
#' @slot sig.genes ABC
#' @slot distribution The distribution function to be used in the glm model.
#'
#' @name scMaSigProClass
#' @aliases scMaSigProClass-class
#' @rdname scMaSigProClass-class
#' @exportClass scMaSigProClass
#' @importFrom methods is new as
#' @keywords classes

# SigClass
setClass(
  "sigClass",
  representation(
    sig.genes = "list",
    feature.clusters = "list"
  ),
  validity = function(object) {
    if (!is.list(object@sig.genes)) {
      stop("sig.genes slot must be a list")
    }
  },
  prototype = list(
    sig.genes = list(),
    feature.clusters = list()
  )
)

setClass(
  "scMaSigProClass",
  representation(
    sparse = "SingleCellExperiment",
    scPVector = "scPVectorClass",
    scTFit = "scTFitClass",
    dense = "SingleCellExperiment",
    edesign = "edesignClass",
    param = "paramClass",
    sig.genes = "sigClass",
    distribution = "ANY"
  ),
  validity = function(object) {
    # Check sparse slot
    if (!validObject(object@sparse)) {
      stop("sparse slot is not a valid SingleCellExperiment object.")
    }

    # Check scPVectorClass slot
    if (!validObject(object@scPVector)) {
      stop("scPVector slot is not a valid scPVectorClass object.")
    }

    # Check scTFitClass slot
    if (!validObject(object@scTFit)) {
      stop("scTFitClass slot is not a valid scTFitClass object.")
    }

    # Check dense slot
    if (!validObject(object@dense)) {
      stop("dense slot is not a valid SingleCellExperiment object.")
    }

    # Check edesignClass slot
    if (!validObject(object@edesign)) {
      stop("edesign slot is not a valid edesignClass object.")
    }

    # Check paramClass slot
    if (!validObject(object@param)) {
      stop("param slot is not a valid paramClass object.")
    }
    # Check paramClass slot
    if (!validObject(object@param)) {
      stop("'sig.genes' slot is not a valid paramClass object.")
    }
  },
  prototype = list(
    scPVector = new("scPVectorClass"), # Assuming you've defined scPVectorClass with its prototype
    scTFit = new("scTFitClass"), # Assuming you've defined scTFitClass with its prototype
    param = new("paramClass"), # Assuming you've defined scTFitClass with its prototype
    sig.genes = new("sigClass"),
    distribution = negative.binomial(theta = 10)
  )
)

scMaSigProClass <- function(sparse = new("SingleCellExperiment"), # Remove default sparse
                            scPVector = new("scPVectorClass"),
                            scTFit = new("scTFitClass"),
                            dense = new("SingleCellExperiment"),
                            edesign = new("edesignClass"),
                            param = new("paramClass"),
                            sig.genes = new("sigClass"),
                            distribution = negative.binomial(theta = 10)) {
  new("scMaSigProClass",
    sparse = sparse,
    scPVector = scPVector,
    scTFit = scTFit,
    dense = dense,
    edesign = edesign,
    paramClass = param,
    sig.genes = sig.genes, # new("sigClass"),
    distribution = distribution
  )
}

setMethod(
  "show",
  "scMaSigProClass",
  function(object) {
    .scmp_show(object)
  }
)
