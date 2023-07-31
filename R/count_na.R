#' Count NA Function
#'
#' Description:
#' This function calculates the number of non-missing (non-NA) elements in a vector.
#'
#' Usage:
#' count.na(x)
#'
#' Arguments:
#'   \describe{
#'     \item{x}{A vector for which you want to count non-missing (non-NA) elements.}
#'   }
#'
#' Details:
#' The function calculates the number of non-missing elements in the input vector
#' by subtracting the count of missing (NA) elements from the total length of the vector.
#'
#' Value:
#' The function returns a numeric value representing the number of non-missing elements
#' in the input vector.
#'
#' @examples
#' # Sample vector
#' x <- c(1, 2, NA, 4, NA, 6, 7, NA)
#'
#' # Count non-missing elements in the vector
#' result <- count.na(x)
#' print(result)
#'
#' @keywords internal
count.na <- function(x) {
  return(length(x) - length(x[is.na(x)]))
}
