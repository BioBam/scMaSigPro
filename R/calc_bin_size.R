#' @title Calculate Bin Size Function
#'
#' @description
#' This function calculates the size of a bin based on the number of elements in
#' the "cluster.members" column of the input data frame.
#'
#' @param x A data frame containing the "cluster.members" column.
#'
#' @return A numeric value representing the size of the bin (number of elements
#' in the "cluster.members" column).
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom stringr str_split
#'
#' @keywords internal
# Define a function 'calc_bin_size' which takes a data frame 'x' as input
calc_bin_size <- function(x) {
  # Use the 'str_split' function from the 'stringr' package to split the 'cluster.members' column
  # of the input data frame 'x' by the '|' character.
  # This returns a list where each element is a vector of the split strings.
  # 'c()' is used to concatenate these vectors into a single vector.
  # Finally, 'length' is used to get the length of this vector (i.e., the number of split strings),
  # which is stored in the 'size' variable.
  size <- length(c(str_split(x[["cluster.members"]], "\\|"))[[1]])

  # Convert the 'size' variable to a numeric value and return it as the result of the function
  return(as.numeric(size))
}
