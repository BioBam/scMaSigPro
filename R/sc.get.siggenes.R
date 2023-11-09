#' @title Extract significant genes for sets of variables in time series gene expression experiments
#'
#' @description
#' `sc.get.siggenes()` creates lists of significant genes for a set of variables
#' whose significance value has been computed with the \code{sc.T.fit} function.
#'
#' @param scmpObj Object of Class \code{\link{scMaSigProClass}} in which the
#' \code{sc.T.fit} has been run.
#' @param rsq Cut-off level at the R-squared value for the stepwise regression fit.
#' @param includeInflu description
#' Only genes with R-squared more than 'rsq' are selected. (Default = 0.7).
#' @param vars Variables for which to extract significant genes. There are 3 possible values:
#' \itemize{
#'   \item \code{"all"}: generates one single matrix or gene list with all significant genes.
#'   \item \code{"each"}: generates as many significant genes extractions as variables in the general regression model.
#'   \item \code{"groups"}: generates a significant genes extraction for each experimental group.
#' }
#' @param significant.intercept Experimental groups for which significant intercept coefficients are considered.
#'
#' @details
#' Refer to the function description for details on the arguments and their usage.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{summary}: A vector or matrix listing significant genes for the variables given by the function parameters.
#'   \item \code{sig.genes}: A list with detailed information on the significant genes found for the variables given by the function parameters.
#'   Each element of the list is also a list containing:
#'     \describe{
#'       \item{\code{sig.profiles}:}{Expression values of significant genes.}
#'       \item{\code{coefficients}:}{Regression coefficients of the adjusted models.}
#'       \item{\code{group.coeffs}:}{Regression coefficients of the implicit models of each experimental group.}
#'       \item{\code{sig.pvalues}:}{P-values of the regression coefficients for significant genes.}
#'       \item{\code{g}:}{Number of genes.}
#'       \item{\code{...}:}{Arguments passed by previous functions.}
#'     }
#'   }
#'
#' @references
#' Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. (2006).
#' maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
#' Bioinformatics, 22(9), 1096-1102. \url{https://doi.org/10.1093/bioinformatics/btl056}
#'
#' @author Ana Conesa and Maria Jose Nueda (mj.nueda@@ua.es)
#'
#' @examples
#' # Example usage of the function can be placed here.
#'
#' @keywords manip
#' @importFrom maSigPro get.siggenes
#'
#' @export
sc.get.siggenes <- function(scmpObj, rsq = 0.7,
                            vars = c("all", "each", "groups"),
                            significant.intercept = "dummy",
                            includeInflu = TRUE) {
  # Check Validity of the object
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigProClass'"
  )

  # Create a named tstep
  tstep <- list(
    dis = scmpObj@edesign@dis,
    edesign = scmpObj@edesign@edesign,
    groups.vector = scmpObj@scTFit@groups.vector,
    sol = showSol(scmpObj, return = TRUE, view = FALSE, includeInflu = includeInflu),
    coefficients = showCoeff(scmpObj, return = TRUE, view = FALSE, includeInflu = includeInflu),
    sig.profiles = showSigProf(scmpObj, return = TRUE, view = FALSE, includeInflu = includeInflu),
    group.coeffs = scmpObj@scTFit@group.coeffs
  )


  sig.genes.s3 <- get.siggenes(tstep,
    rsq = rsq, add.IDs = FALSE, IDs = NULL, matchID.col = 1,
    only.names = FALSE, vars = vars,
    significant.intercept = significant.intercept,
    groups.vector = NULL, trat.repl.spots = "none",
    r = 0.7
  )

  # Call the maSigPro get sig
  # Create Object
  siggenes.object <- new("sigClass",
    summary = sig.genes.s3$summary,
    sig.genes = sig.genes.s3$sig.genes
  )

  # Update the slot
  scmpObj@sig.genes <- siggenes.object

  # Return
  return(scmpObj)
}
