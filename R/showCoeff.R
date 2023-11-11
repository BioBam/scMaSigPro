#' Show or Return the Solution
#'
#' This function is used to view or return the solution of the provided scMaSigPro object.
#'
#' @param scmpObj an object of class 'scMaSigProClass'. This object should contain the computed solution.
#' @param view logical, whether to view the solution. If TRUE (default), the solution is displayed.
#' @param includeInflu description
#' @param return logical, whether to return the solution. If FALSE (default), the solution is not returned.
#'
#' @return The computed solution as a data.frame if return is set to TRUE.
#' If return is FALSE, the function does not return anything.
#'
#' @export
showCoeff <- function(scmpObj, view = FALSE, return = TRUE, includeInflu = FALSE) {
  # Check Object Validity
  assert_that(is(scmpObj, "scMaSigProClass"),
    msg = "Please provide object of class 'scMaSigPro'"
  )

  # Check if the sol exist
  assert_that(!all(dim(scmpObj@scTFit@sol) == c(0, 0)),
    msg = "Coeff is not computed yet"
  )

  # Extract
  coefficients <- scmpObj@scTFit@coefficients %>% as.data.frame()

  if (!includeInflu) {
    influ.gene <- colnames(showInflu(scmpObj, return = TRUE, view = FALSE))
    coefficients <- coefficients[!(rownames(coefficients) %in% influ.gene), ]
  }

  # If viewing is requested
  if (view) {
    View(coefficients)
  }

  # If requested
  if (return) {
    return(coefficients)
  }
}
