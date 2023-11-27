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
#' @slot poly_degree Integer with the polynomial degree to fit the regression. 1 
#' specifies a linear regression, 2 a quadratic regression, etc.  
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
    edesign = "matrix",
    poly_degree = "integer"
  ),
  prototype = list(
    dis = matrix(NA, nrow = 0, ncol = 0),
    groups.vector = character(),
    edesign = matrix(NA, nrow = 0, ncol = 0),
    poly_degree = as.integer(2)
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
    if (!is.integer(object@poly_degree)) {
      stop("edesign slot is not a valid integer")
    }
    TRUE
  }
)
