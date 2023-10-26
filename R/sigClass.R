#' @title Class "sigClass"
#'
#' @description A class for Storing significant genes.
#'
#' @slot summary A dataframe listing significant genes for the variables given 
#' by the function parameters.
#' @slot sig.genes A list with detailed information on the significant genes found 
#' for the variables given by the function parameters. Each element of the list is also a list containing:
#'     \itemize{
#'       \item{\code{sig.profiles}:}{Expression values of significant genes.}
#'       \item{\code{coefficients}:}{Regression coefficients of the adjusted models.}
#'       \item{\code{group.coeffs}:}{Regression coefficients of the implicit models of each experimental group.}
#'       \item{\code{sig.pvalues}:}{P-values of the regression coefficients for significant genes.}
#'       \item{\code{g}:}{Number of genes.}
#'       \item{\code{...}:}{Arguments passed by previous functions.}
#'     }
#'
#' @name sigClass
#' @aliases sigClass-class
#' @rdname sigClass-class
#' @exportClass sigClass
#' @importFrom methods is new as
#' @keywords classes
#'
setClass(
  "sigClass",
  representation(
    summary = "data.frame",
    sig.genes = "list"
  ),
  validity = function(object) {
    if (!is.data.frame(object@summary)) {
      stop("summary slot must be a data.frame")
    }
    if (!is.list(object@sig.genes)) {
      stop("sig.genes slot must be a list")
    }
  },
  prototype = list(
    summary = data.frame(),
    sig.genes = list()
  )
)
