#' A class to store regression fit results for time series gene expression experiments
#'
#' The scPVectorClass is designed to hold the results of a regression fit for each gene in a time series gene expression experiment.
#' It contains slots to store various results, such as expression values for significant genes, computed p-values, adjusted p-values,
#' number of input genes, number of genes taken in the regression fit, and more. This class is useful for analyzing time course
#' microarray experiments and identifying significant differential expression profiles.
#'
#' @slot SELEC dgCMatrix containing the expression values for significant genes.
#' @slot p.vector Numeric vector containing the computed p-values.
#' @slot p.adjusted Numeric vector of FDR-adjusted p-values.
#' @slot FDR P-value at FDR \code{Q} control when Benjamini & Hochberg (BH) correction is used.
#' @slot dis Data frame containing the matrix used in the regression fit.
#' @slot groups.vector Character list containing groups information.
#' @slot family Distribution function to be used in the glm model. If NULL, the
#' family will be \code{negative.binomial(theta)} when \code{counts = TRUE} or
#' \code{gaussian()} when \code{counts = FALSE}.
#'
#' @name scPVectorClass
#' @aliases scPVectorClass-class
#' @rdname scPVectorClass-class
#' @exportClass scPVectorClass
#' @keywords classes regression
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
#' @importFrom stats family gaussian poisson
#' @importFrom utils data combn
#' @importFrom MASS negative.binomial

# Define the scPVectorClass with the following slots:
setClass("scPVectorClass",
  slots = c(
    SELEC = "dgCMatrix", # Matrix containing the expression values for significant genes
    p.vector = "numeric", # Matrix containing the computed p-values
    p.adjusted = "numeric", # Vector of FDR-adjusted p-values
    FDR = "numeric", # P-value at FDR Q control when Benjamini & Hochberg (BH) correction is used
    dis = "data.frame", # Design matrix used in the regression fit
    groups.vector = "character", # List containing groups information
    family = "ANY" # Distribution function to be used in the glm model
  ),
  validity = function(object) {
    # Check for slot SELEC
    if (!validObject(object@SELEC)) {
      stop("Slot 'SELEC' must be a dgCMatrix")
    }

    # Check for slot sc.p.vector
    if (!is.numeric(object@p.vector)) {
      stop("Slot 'sc.p.vector' must be a numeric")
    }

    # Check for slot p.adjusted
    if (!is.numeric(object@p.adjusted)) {
      stop("Slot 'p.adjusted' must be numeric.")
    }

    # Check for slot FDR
    if (!is.numeric(object@FDR)) {
      stop("Slot 'FDR' must be numeric.")
    }

    # Check for slot dis
    if (!is.data.frame(object@dis)) {
      stop("Slot 'dis' must be a data frame.")
    }

    # Check for slot groups.vector
    if (!is.character(object@groups.vector)) {
      stop("Slot 'groups.vector' must be a list.")
    }

    # # Check for slot family
    # if (!inherits(object@family, "ANY")) {
    #   stop("Slot 'family' must inherit from class 'family'.")
    # }
    # Return TRUE if all checks pass
    TRUE
  },
  prototype = list(
    SELEC = as(matrix(NA, nrow = 0, ncol = 0), "dgCMatrix"), # Empty matrix for SELEC
    p.vector = numeric(0), # Empty matrix for sc.p.vector
    p.adjusted = numeric(0), # Empty numeric vector for p.adjusted
    FDR = 0, # Default FDR value is 0
    dis = data.frame(), # Empty data frame for dis
    groups.vector = character() # Empty list for groups.vector
  )
)
