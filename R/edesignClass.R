#' edesignClass
#'
#' An S4 class to represent an edesign object with associated data.
#' This class contains three slots: dis, groups.vector, and edesign.
#'
#' @slot dis A data frame.
#' @slot groups.vector A character vector.
#' @slot edesign A data frame.
#' @slot poly_degree Polynomial degree
#'
#' @section Validity:
#'   Valid objects must have:
#'   \itemize{
#'     \item{dis}{A valid data frame.}
#'     \item{groups.vector}{A valid character vector.}
#'     \item{edesign}{A valid data frame.}
#'   }
#'
#' @exportClass edesignClass
setClass(
  "edesignClass",
  representation(
    dis = "matrix",
    groups.vector = "character",
    edesign = "data.frame",
    poly_degree = "integer"
  ),
  prototype = list(
    dis = matrix(NA, nrow = 0, ncol = 0),
    groups.vector = character(),
    edesign = data.frame(),
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
      stop("poly_degree slot is not a valid data frame.")
    }
    if (!is.integer(object@poly_degree)) {
      stop("edesign slot is not a valid integer")
    }
    TRUE
  }
)
