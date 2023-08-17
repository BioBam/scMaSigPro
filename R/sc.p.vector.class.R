#' A class to store regression fit results for time series gene expression experiments
#'
#' The scPVectorClass is designed to hold the results of a regression fit for each gene in a time series gene expression experiment.
#' It contains slots to store various results, such as expression values for significant genes, computed p-values, adjusted p-values,
#' number of input genes, number of genes taken in the regression fit, and more. This class is useful for analyzing time course
#' microarray experiments and identifying significant differential expression profiles.
#'
#' @slot SELEC Matrix containing the expression values for significant genes
#' @slot sc.p.vector Matrix containing the computed p-values
#' @slot p.adjusted Vector of FDR-adjusted p-values
#' @slot G Total number of input genes
#' @slot g Number of genes taken in the regression fit
#' @slot FDR P-value at FDR \code{Q} control when Benjamini & Hochberg (BH) correction is used
#' @slot i Number of significant genes
#' @slot dis Design matrix used in the regression fit
#' @slot dat Matrix of expression value data used in the regression fit
#' @slot min.obs Minimum value to estimate the model (degree+1) x Groups + 1. Default is 6.
#' @slot Q Significance level (default is 0.05)
#' @slot groups.vector List containing groups information
#' @slot edesign Experimental design data frame
#' @slot family Distribution function to be used in the glm model. If NULL, the family will be \code{negative.binomial(theta)} when \code{counts = TRUE} or \code{gaussian()} when \code{counts = FALSE}.
#'
#' @export
#' @keywords regression
#'
#' @examples
#' #### GENERATE TIME COURSE DATA
#' ## generates n random gene expression profiles of a data set with
#' ## one control plus 3 treatments, 3 time points and r replicates per time point.
#'
#' # ... (add example usage of the function here)
#'
#' @author Priyansh Srivastava <spriyansh29@gmail.com>
#' @seealso \code{\link{T.fit}}, \code{\link{lm}}
#' @importFrom stats family gaussian poisson negative.binomial
#' @importFrom utils data combn

# Define the scPVectorClass with the following slots:
setClass("scPVectorClass",
  slots = c(
    SELEC = "matrix", # Matrix containing the expression values for significant genes
    sc.p.vector = "matrix", # Matrix containing the computed p-values
    p.adjusted = "numeric", # Vector of FDR-adjusted p-values
    G = "integer", # Total number of input genes
    g = "integer", # Number of genes taken in the regression fit
    FDR = "numeric", # P-value at FDR Q control when Benjamini & Hochberg (BH) correction is used
    i = "integer", # Number of significant genes
    dis = "data.frame", # Design matrix used in the regression fit
    dat = "matrix", # Matrix of expression value data used in the regression fit
    min.obs = "numeric", # Minimum value to estimate the model (degree+1) x Groups + 1
    Q = "numeric", # Significance level (default is 0.05)
    groups.vector = "character", # List containing groups information
    edesign = "matrix", # Experimental design data frame
    family = "ANY" # Distribution function to be used in the glm model
  ),
  validity = function(object) {
    # Check for slot SELEC
    if (!is.matrix(object@SELEC)) {
      stop("Slot 'SELEC' must be a matrix.")
    }

    # Check for slot sc.p.vector
    if (!is.matrix(object@sc.p.vector)) {
      stop("Slot 'sc.p.vector' must be a matrix.")
    }

    # Check for slot p.adjusted
    if (!is.numeric(object@p.adjusted)) {
      stop("Slot 'p.adjusted' must be numeric.")
    }

    # Check for slot G
    if (!is.integer(object@G)) {
      stop("Slot 'G' must be an integer.")
    }

    # Check for slot g
    if (!is.integer(object@g)) {
      stop("Slot 'g' must be an integer.")
    }

    # Check for slot FDR
    if (!is.numeric(object@FDR)) {
      stop("Slot 'FDR' must be numeric.")
    }

    # Check for slot i
    if (!is.integer(object@i)) {
      stop("Slot 'i' must be an integer.")
    }

    # Check for slot dis
    if (!is.data.frame(object@dis)) {
      stop("Slot 'dis' must be a data frame.")
    }

    # Check for slot dat
    if (!is.matrix(object@dat)) {
      stop("Slot 'dat' must be a matrix.")
    }

    # Check for slot min.obs
    if (!is.numeric(object@min.obs)) {
      stop("Slot 'min.obs' must be an integer.")
    }

    # Check for slot Q
    if (!is.numeric(object@Q)) {
      stop("Slot 'Q' must be numeric.")
    }

    # Check for slot groups.vector
    if (!is.character(object@groups.vector)) {
      stop("Slot 'groups.vector' must be a list.")
    }

    # Check for slot edesign
    if (!is.matrix(object@edesign)) {
      stop("Slot 'edesign' must be a data frame.")
    }

    # # Check for slot family
    # if (!inherits(object@family, "ANY")) {
    #   stop("Slot 'family' must inherit from class 'family'.")
    # }
    # Return TRUE if all checks pass
    TRUE
  },
  prototype = list(
    SELEC = matrix(NA, nrow = 0, ncol = 0), # Empty matrix for SELEC
    sc.p.vector = matrix(NA, nrow = 0, ncol = 0), # Empty matrix for sc.p.vector
    p.adjusted = numeric(0), # Empty numeric vector for p.adjusted
    G = 0L, # Default G value is 0
    g = 0L, # Default g value is 0
    FDR = 0, # Default FDR value is 0
    i = 0L, # Default i value is 0
    dis = data.frame(), # Empty data frame for dis
    dat = matrix(NA, nrow = 0, ncol = 0), # Empty matrix for dat
    min.obs = 6, # Default min.obs value is 0
    Q = 0.05, # Default Q value is 0
    groups.vector = character(), # Empty list for groups.vector
    edesign = matrix(NA, nrow = 0, ncol = 0), # Empty data frame for edesign
    family = gaussian() # Default family value is NULL
  )
)

setGeneric("family", function(object) standardGeneric("family"))

setMethod("family", "scPVectorClass", function(object) {
  if (isS4(object)) {
    return(object@family)
  } else if (is(object, "list")) {
    return(object$family)
  } else {
    stop("Object must be of class 'scPVectorClass' (S4) or list (S3).")
  }
})
