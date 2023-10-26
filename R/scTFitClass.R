#' @title Class "scTFitClass"
#'
#' @description A class for fitting a model.
#'
#' @slot sol A data frame for summary results of the stepwise regression. 
#' For each selected gene the following values are given:
#'   \itemize{
#'     \item{p-value}{of the regression ANOVA.}
#'     \item{R-squared}{of the model.}
#'     \item{p-value}{of the regression coefficients of the selected variables.}
#'   }
#' @slot coefficients A data frame containing regression coefficients for the adjusted models.
#' @slot group.coeffs A matrix with the coefficients of the implicit models of each experimental group.
#' @slot t.score A data frame containing tscores for each covariate in polynomial glm.
#' @slot variables A character vector containing the variables in the complete regression model.
#' @slot groups.vector A character vector containing the branching path.
#' @slot influ.info A matrix with genes containing influencial data.
#'
#' @name scTFitClass
#' @aliases scTFitClass-class
#' @rdname scTFitClass-class
#' @exportClass scTFitClass
#' @importFrom methods is new
#' @keywords classes


setClass(
  "scTFitClass",
  representation(
    sol = "data.frame",
    coefficients = "data.frame",
    group.coeffs = "matrix",
    t.score = "data.frame",
    variables = "character",
    groups.vector = "character",
    influ.info = "matrix"
  ),
  validity = function(object) {
    if (!is.data.frame(object@sol)) {
      stop("sol slot must be a data.frame")
    }
    if (!is.data.frame(object@coefficients)) {
      stop("coefficients slot must be a data.frame")
    }
    if (!is.matrix(object@group.coeffs)) {
      stop("group.coeffs slot must be a matrix.")
    }
    if (!is.data.frame(object@t.score)) {
      stop("t.score slot must be a data.frame")
    }
    if (!is.character(object@variables)) {
      stop("variables slot must be a character vector.")
    }
    if (!is.character(object@groups.vector)) {
      stop("groups.vector slot must be a character.")
    }
    if (!is.matrix(object@influ.info)) {
      stop("influ.info slot must be a matrix.")
    }
  },
  prototype = list(
    sol = data.frame(),
    coefficients = data.frame(),
    group.coeffs = matrix(0, nrow = 1, ncol = 1),
    t.score = data.frame(),
    variables = "not_selected",
    groups.vector = character(),
    influ.info = matrix(NA, nrow = 0, ncol = 0)
  )
)
