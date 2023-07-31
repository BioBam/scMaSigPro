#' Average rows by match and index
#'
#' \code{average.rows} matches rownames of a matrix to a \code{match} vector and performs
#' averaging of the rows by the index provided by an \code{index} vector.
#'
#' @title Average rows by match and index
#'
#' @description
#' \code{average.rows} matches rownames of a matrix to a \code{match} vector and performs
#' averaging of the rows by the index provided by an \code{index} vector.
#'
#' @usage
#' average.rows(x, index, match, r = 0.7)
#'
#' @arguments
#' \item{x}{A matrix.}
#' \item{index}{Index vector indicating how rows must be averaged.}
#' \item{match}{Match vector for indexing rows.}
#' \item{r}{Minimal correlation value between rows to compute average.}
#'
#' @details
#' Rows will be averaged only if the Pearson correlation coefficient between all rows of each given
#' index is greater than 'r'. If not, that group of rows is discarded in the result matrix.
#'
#' @value
#' A matrix of averaged rows.
#'
#' @author Ana Conesa and Maria Jose Nueda
#' @email mj.nueda@ua.es
#'
#' @examples
#' # Create data matrix for row averaging
#' x <- matrix(rnorm(30), nrow = 6, ncol = 5)
#' rownames(x) <- paste("ID", c(1, 2, 11, 12, 19, 20), sep = "")
#' i <- paste("g", rep(c(1:10), each = 2), sep = "") # index vector
#' m <- paste("ID", c(1:20), sep = "") # match vector
#' average.rows(x, i, m, r = 0)
#'
#' @keyword arith
#'
#' @export

# average.rows: Compute row-wise averages based on specified conditions.
#
# This function takes a matrix 'x' and performs row-wise averaging based on the 'index' and 'match' vectors.
# Rows in 'x' are grouped together according to their matching row names specified in the 'index' vector,
# and the row-wise averages are computed for each group.
#
# Arguments:
#   x: A numeric matrix for which the row-wise averages are calculated.
#   index: A vector used to group rows of 'x' for averaging based on matching row names.
#   match: A character or factor vector representing the row names to match in 'x' with 'index'.
#   r: A numeric value (default is 0.7) representing the correlation threshold used to check for
#      correlation among rows in 'x'. Rows with a correlation coefficient less than 'r' will not be
#      included in the average calculation and will result in NA values in the output.
#
# Returns:
#   A data frame containing the row-wise averages of 'x' after applying the specified conditions.
#   Rows with correlation less than 'r' are omitted.
#
average.rows <- function(x, index, match, r = 0.7) {
  # Extract matching row names from 'x' based on 'index' and 'match' vectors.
  name <- index[match(rownames(x), as.character(match))]

  # Get unique group names.
  uninames <- unique(name)

  # Initialize a matrix to store the row-wise averages.
  prom <- matrix(0, nrow = length(uninames), ncol = ncol(x))
  colnames(prom) <- colnames(x)
  rownames(prom) <- uninames

  # Loop through each group and calculate the row-wise averages.
  for (i in 1:length(uninames)) {
    # Extract the rows belonging to the current group.
    c <- x[name == uninames[i], ]

    # Check if the correlation between rows is less than 'r'.
    # If so, set the entire row to NA in the 'prom' matrix.
    if (is.element(TRUE, cor(t(c), use = "pairwise.complete.obs") < r)) {
      prom[i, ] <- c(rep(NA, ncol(c)))
      rownames(prom)[i] <- NA
    } else {
      # If the correlation is greater than or equal to 'r', calculate the row-wise average.
      prom[i, ] <- apply(as.matrix(c), 2, mean, na.rm = TRUE)
    }
  }

  # Convert the resulting matrix to a data frame.
  prom <- as.data.frame(prom)

  # Remove rows with NA row names (rows that were discarded due to low correlation).
  prom <- prom[!is.na(rownames(prom)), ]

  # Return the data frame containing the row-wise averages.
  prom
}
