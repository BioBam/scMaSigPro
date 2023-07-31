#' Show Solution Matrix
#'
#' This function extracts and displays the solution matrix from the output of the ScMaSigPro
#' analysis (\code{scMaSigPro_t_fit} object).
#'
#' @param scMaSigPro.obj A \code{scMaSigPro_t_fit} object containing the results of the ScMaSigPro analysis.
#' @param view Logical value indicating whether to view the solution matrix using a data viewer (e.g., RStudio Viewer).
#'
#' @return The solution matrix extracted from the ScMaSigPro analysis results.
#'
#' @importFrom stats sapply t
#'
#' @keywords internal
#' @export
showSol <- function(scMaSigPro.obj, view = TRUE) {
  # Extract Coefficients
  sol <- sapply(scMaSigPro.obj@tFit@model.attributes, simplify = TRUE, function(x) {
    return(x[["sol"]])
  })

  # Transpose
  sol <- t(sol)

  # If viewing is requested
  if (view) {
    View(sol)
  }

  # Return
  return(sol)
}
