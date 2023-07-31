#' Show T-scores Matrix
#'
#' This function extracts and displays the T-scores matrix from the output of the ScMaSigPro
#' analysis (\code{scMaSigPro_t_fit} object).
#'
#' @param scMaSigPro.obj A \code{scMaSigPro_t_fit} object containing the results of the ScMaSigPro analysis.
#' @param view Logical value indicating whether to view the T-scores matrix using a data viewer (e.g., RStudio Viewer).
#'
#' @return The T-scores matrix extracted from the ScMaSigPro analysis results.
#'
#' @importFrom stats sapply t
#'
#' @keywords internal
#' @export
showTscores <- function(scMaSigPro.obj, view = TRUE) {
  # Extract Coefficients
  t_value <- sapply(scMaSigPro.obj@tFit@model.attributes, simplify = TRUE, function(x) {
    return(x[["t_value"]])
  })

  # Transpose
  t_value <- t(t_value)

  # If viewing is requested
  if (view) {
    View(t_value)
  }

  # Return
  return(t_value)
}
