#' Calculate Bin Size Function
#'
#' This function calculates the size of a bin based on the number of elements in
#' the "cluster.members" column of the input data frame.
#'
#' @param x A data frame containing the "cluster.members" column.
#'
#' @return A numeric value representing the size of the bin (number of elements
#' in the "cluster.members" column).
#'
#' @examples
#' # Sample data frame
#' data <- data.frame(
#'   cluster.members = c("A|B|C", "D|E", "F|G|H|I")
#' )
#'
#' # Calculate the bin size using the calc_bin_size function
#' result <- calc_bin_size(data)
#' print(result)
#'
#' @importFrom stringr str_split
#'
#' @keywords internal
calc_bin_size <- function(x) {
  size <- length(
    c(str_split(x[["cluster.members"]], "\\|"))[[1]]
  )
  return(as.numeric(size))
}
