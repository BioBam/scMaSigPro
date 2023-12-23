#' @title Extract significant genes for sets of variables in time series gene expression experiments
#'
#' @description
#' `sc.filter()` creates lists of significant genes for a set of variables
#' whose significance value has been computed with the \code{sc.t.fit} function.
#'
#' @param scmpObj Object of Class \code{\link{ScMaSigPro}} in which the
#' \code{sc.t.fit} has been run.
#' @param rsq Cut-off level at the R-squared value for the stepwise regression fit.
#' @param includeInflu description
#' @param Q overall model significance
#' @param term.Q Term wise significance
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
sc.filter <- function(scmpObj,
                      rsq = 0.7,
                      Q = scmpObj@param@Q,
                      vars = c("all", "each", "groups"),
                      significant.intercept = "dummy",
                      term.Q = 0.05,
                      includeInflu = TRUE) {
  # Check Validity of the object
  assert_that(is(scmpObj, "ScMaSigPro"),
    msg = "Please provide object of class 'ScMaSigPro'"
  )

  assert_that(
    all(
      vars %in% c("all", "each", "groups")
    ),
    msg = "Invalid selction for 'vars'"
  )

  # Initiate empty list
  sig.list <- list()

  # Extract sol
  sol <- showSol(scmpObj, includeInflu = includeInflu)

  # Replace the NAs with corresponding values
  sol[is.na(sol$`p-value`), "p-value"] <- 1
  sol[is.na(sol$`R-squared`), "R-squared"] <- 0
  sol[, !(colnames(sol) %in% c("p-value", "R-squared"))][is.na(sol[, !(colnames(sol) %in% c("p-value", "R-squared"))])] <- 1

  # Filter based on RSQ and P-value
  sol <- sol[sol$`p-value` <= Q, ]
  sol <- sol[sol$`R-squared` >= rsq, ]

  # If`all`
  if (vars == "all") {
    # Return a list of gene names
    sig.list[["all"]] <- rownames(sol)
  } else if (vars == "each") {
    # Subset
    sol.sub <- sol[, !(colnames(sol) %in% c("p-value", "R-squared")), drop = FALSE]

    # Get the genes
    sig.list <- apply(sol.sub, 2, function(column) row.names(sol.sub)[which(column <= term.Q)])
    names(sig.list) <- stringr::str_remove(names(sig.list), "p.valor_")
  } else if (vars == "groups") {
    # Subset
    sol.sub <- sol[, !(colnames(sol) %in% c("p-value", "R-squared")), drop = FALSE]

    # Get group_vector, from T fit
    group_vector <- scmpObj@estimate@path

    # Based on the dummy, none and all
    if (significant.intercept == "all") {
      select_cols <- seq(from = 1, to = ncol(sol.sub))
    } else if (significant.intercept == "dummy") {
      select_cols <- seq(from = 2, to = ncol(sol.sub))
    } else if (significant.intercept == "none") {
      select_cols <- seq(from = 3, to = ncol(sol.sub))
    }

    # subset groupVector and sol
    group_vector <- group_vector[select_cols]
    sol.sub <- sol.sub[, select_cols, drop = FALSE]

    # Create Named list for traversal
    group_term_list <- as.list(as.vector(colnames(sol.sub)))
    names(group_term_list) <- group_vector

    # Available groups
    avail_groups <- unique(group_vector)

    # Get per path
    sig.list <- lapply(avail_groups, function(group_i) {
      # Get names per path
      group_term_list_sub <- group_term_list[names(group_term_list) %in% group_i]

      # get subset
      sol.sub.term <- sol.sub[, unlist(group_term_list_sub), drop = FALSE]

      # Get group wise genes
      sig.list.tmp <- row.names(sol.sub.term)[apply(sol.sub.term, 1, function(x) any(x <= term.Q))]

      # Return
      return(sig.list.tmp)
    })

    # Add names
    names(sig.list) <- avail_groups
  }

  # # Create a named tstep
  ## Donot remove this ##
  # tstep <- list(
  #   dis = scmpObj@design@predictor,
  #   edesign = scmpObj@design@alloc,
  #   groups.vector = scmpObj@estimate@path,
  #   sol = showSol(scmpObj, return = TRUE, view = FALSE, includeInflu = includeInflu),
  #   coefficients = showCoeff(scmpObj, return = TRUE, view = FALSE, includeInflu = includeInflu),
  #   sig.profiles = showSigProf(scmpObj, return = TRUE, view = FALSE, includeInflu = includeInflu),
  #   group.coeffs = scmpObj@estimate@group.coeffs
  # )

  # sig.genes.s3 <- get.siggenes(tstep,
  #   rsq = rsq, add.IDs = FALSE, IDs = NULL, matchID.col = 1,
  #   only.names = FALSE, vars = vars,
  #   significant.intercept = significant.intercept,
  #   groups.vector = NULL, trat.repl.spots = "none",
  #   r = 0.7
  # )

  # Create Object
  # siggenes.object <- new("Significant",
  #   summary = sig.genes.s3$summary,
  #   sig.genes = sig.genes.s3$sig.genes
  # )

  # Create Object
  siggenes.object <- new("Significant",
    sig.genes = sig.list
  )

  # Update the slot
  scmpObj@sig.genes <- siggenes.object

  # Return
  return(scmpObj)
}
