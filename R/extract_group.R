#' Extract Group Components Function
#'
#' This function takes a coefficient vector and groups as input and extracts the
#' group components A and B based on the provided grouping information. It
#' assigns appropriate names to the components and combines them into a single
#' vector.
#'
#' @param coeff.vec A numeric vector representing the coefficient estimates.
#' @param groups A character vector of length 2, containing the names of the
#'               two groups.
#'
#' @return A named numeric vector containing the extracted group components A
#'         and B, with appropriate names.
#'
#' @examples
#' # Assuming coeff.vec is a numeric vector of coefficient estimates, and
#' # groups is a character vector of length 2.
#' # Extract the group components
#' components <- extract.group.components(coeff.vec, groups = c("GroupA", "GroupB"))
#' print(components)
#'
#' @keywords internal
#' @export
extract.group.components <- function(coeff.vec, groups) {
  # Convert NA to zero
  coeff.vec[is.na(coeff.vec)] <- 0

  # Extract A
  A <- coeff.vec[c(TRUE, FALSE)]

  # Extract B
  B <- A + coeff.vec[c(FALSE, TRUE)]

  # Names
  names(A) <- paste(groups[1], "beta", 0:(length(A) - 1), sep = "_")
  names(B) <- paste(groups[2], "beta", 0:(length(B) - 1), sep = "_")

  # Cbind
  group.coeff <- c(A, B)

  # Return
  return(group.coeff)
}
