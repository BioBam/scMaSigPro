#' Show Estimated Coefficients Function
#'
#' This function extracts the estimated coefficients (beta estimates) from the
#' scMaSigPro object and returns them as a transposed matrix. The function
#' provides the option to view the matrix using the "View" function.
#'
#' @param scMaSigPro.obj An object of class "scMaSigPro" representing the output
#'                       of the `scMaSigPro` function.
#'
#' @param view Logical value indicating whether to view the matrix of beta
#'             estimates in the RStudio Viewer pane. Default is TRUE.
#'
#' @return A transposed matrix containing the estimated coefficients (beta
#'         estimates) from the scMaSigPro object.
#'
#' @examples
#' # Assuming scMaSigPro.obj is an object returned from the scMaSigPro function
#' # View the estimated coefficients
#' showEstimates(scMaSigPro.obj)
#'
#' # Get the estimated coefficients without viewing
#' estimates <- showEstimates(scMaSigPro.obj, view = FALSE)
#'
#' @importFrom stats View
#'
#' @keywords internal
#' @export
showEstimates <- function(scMaSigPro.obj, view = TRUE) {
  # Extract Coefficients
  beta_estimates <- sapply(scMaSigPro.obj@tFit@model.attributes, simplify = TRUE, function(x) {
    return(x[["beta_estimates"]])
  })

  # Transpose
  beta_estimates <- t(beta_estimates)

  # If viewing is requested
  if (view) {
    View(beta_estimates)
  }

  # Return
  return(beta_estimates)
}
