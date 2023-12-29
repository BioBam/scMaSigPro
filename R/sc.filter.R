#' @title Extract significant genes based on R-Square and P-values
#'
#' @description
#' `sc.filter()` creates lists of significant genes based on user-specified
#' constraints.
#'
#' @importFrom maSigPro get.siggenes
#'
#' @param scmpObj An object of class \code{\link{ScMaSigPro}}.
#' @param rsq Coefficient of determination or R-squared value threshold.
#' @param includeInflu Whether to include genes with influential observations.
#' @param p_value Overall model significance.
#' @param term_p_value Term wise significance.
#' @param vars Variables for which to extract significant genes. See details.
#' @param intercept Specify the branching path treated as reference. See details.
#' (When `vars` equals "groups").
#'
#' @details
#' `vars` Parameter can take one of the following values:
#' \itemize{
#'   \item \code{"all"}: Generates one gene list with all significant genes.
#'   \item \code{"each"}: Generates gene list for each term in the polynomial GLM.
#'   \item \code{"groups"}: Generates gene list for each branching path.
#' }
#'
#' `intercept` Parameter modulates the treatment for intercept coefficients to
#' apply for selecting significant genes when `vars` equals "groups". There are
#' three possible values:
#'  \itemize{
#'   \item \code{"none"}: No significant intercept (differences) are considered.
#'   \item \code{"dummy"}: Includes genes with significant intercept differences
#'   between branching paths.
#'   \item \code{"all"}: When both significant intercept coefficient for the
#'   reference path and significant intercept differences are considered for
#'    selecting significant genes.
#' }
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with updated `Significant`
#' slot.
#'
#' @seealso `maSigPro::get.siggenes()`
#'
#' @references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
#' maSigPro: a Method to Identify Significant Differential Expression Profiles
#' in Time-Course Microarray Experiments. Bioinformatics 22, 1096-1102}
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}, Ana Conesa and
#' Maria Jose Nueda, \email{mj.nueda@@ua.es}
#'
#' @export
sc.filter <- function(scmpObj,
                      rsq = 0.7,
                      p_value = scmpObj@Parameters@p_value,
                      vars = c("all", "each", "groups"),
                      intercept = "dummy",
                      term_p_value = 0.05,
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
  sol <- sol[sol$`p-value` <= p_value, ]
  sol <- sol[sol$`R-squared` >= rsq, ]

  # If`all`
  if (vars == "all") {
    # Return a list of gene names
    sig.list[["all"]] <- rownames(sol)
  } else if (vars == "each") {
    # Subset
    sol.sub <- sol[, !(colnames(sol) %in% c("p-value", "R-squared")), drop = FALSE]

    # Get the genes
    sig.list <- apply(sol.sub, 2, function(column) row.names(sol.sub)[which(column <= term_p_value)])
    names(sig.list) <- stringr::str_remove(names(sig.list), "p.valor_")
  } else if (vars == "groups") {
    # Subset
    sol.sub <- sol[, !(colnames(sol) %in% c("p-value", "R-squared")), drop = FALSE]

    # Get group_vector, from T fit
    group_vector <- scmpObj@Estimate@path

    # Based on the dummy, none and all
    if (intercept == "all") {
      select_cols <- seq(from = 1, to = ncol(sol.sub))
    } else if (intercept == "dummy") {
      select_cols <- seq(from = 2, to = ncol(sol.sub))
    } else if (intercept == "none") {
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
      sig.list.tmp <- row.names(sol.sub.term)[apply(sol.sub.term, 1, function(x) any(x <= term_p_value))]

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
  #   intercept = intercept,
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
    genes = sig.list
  )

  # Update the slot
  scmpObj@Significant <- siggenes.object

  # Return
  return(scmpObj)
}
