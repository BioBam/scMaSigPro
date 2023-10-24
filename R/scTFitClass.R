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
#' @slot sig.profiles A data frame with the expression values for the genes contained in sol.
#' @slot coefficients A data frame containing regression coefficients for the adjusted models.
#' @slot group.coeffs A matrix with the coefficients of the implicit models of each experimental group.
#' @slot t.score A data frame containing tscores for each covariate in polynomial glm.
#' @slot variables A character vector containing the variables in the complete regression model.
#' @slot G An integer representing the total number of input genes.
#' @slot g An numeric value representing the number of genes taken in the regression fit.
#' @slot dat A matrix with the input analysis data.
#' @slot dis A data frame with the regression design.
#' @slot step.method A character specifying the imputed step method for stepwise regression.
#' @slot groups.vector A character vector containing the branching path.
#' @slot edesign Experimental design in matrix format.
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
    sig.profiles = "data.frame",
    coefficients = "data.frame",
    group.coeffs = "matrix",
    t.score = "data.frame",
    variables = "character",
    G = "integer",
    g = "numeric",
    dat = "matrix",
    dis = "data.frame",
    step.method = "character",
    groups.vector = "character",
    edesign = "matrix",
    influ.info = "matrix"
  ),
  validity = function(object) {
    if (!is.data.frame(object@sol)) {
      stop("sol slot must be a data.frame")
    }
    if (!is.data.frame(object@sig.profiles)) {
      stop("sig.profiles slot must be a data.frame")
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
    if (!is.integer(object@G)) {
      stop("G slot must be an integer.")
    }
    if (!is.numeric(object@g)) {
      stop("g slot must be an integer.")
    }
    if (!is.matrix(object@dat)) {
      stop("dat slot must be a matrix.")
    }
    if (!is.data.frame(object@dis)) {
      stop("dis slot must be a data.frame")
    }
    if (!is.character(object@step.method)) {
      stop("step.method slot must be a character.")
    }
    if (!is.character(object@groups.vector)) {
      stop("groups.vector slot must be a character.")
    }
    if (!is.matrix(object@edesign)) {
      stop("edesign slot must be a matrix.")
    }
    if (!is.matrix(object@influ.info)) {
      stop("influ.info slot must be a matrix.")
    }
  },
  prototype = list(
    sol = data.frame(),
    sig.profiles = data.frame(),
    coefficients = data.frame(),
    group.coeffs = matrix(0, nrow = 1, ncol = 1),
    t.score = data.frame(),
    variables = "not_selected",
    G = integer(0),
    g = 0,
    dat = matrix(0, nrow = 1, ncol = 1),
    dis = data.frame(),
    step.method = "backward",
    groups.vector = character(),
    edesign = matrix(0, nrow = 1, ncol = 1),
    influ.info = matrix(NA, nrow = 0, ncol = 0)
  )
)
