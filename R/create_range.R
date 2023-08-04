#' @title Create Range Function
#'
#' @description
#' This function converts a factor column "bin" into a character vector, extracts
#' numeric range values from the character vector, and combines them with
#' additional columns "bin_size" and "binned_time" to return a numeric vector
#' representing the range information.
#'
#' @param x A data frame that should contain the following columns:
#'   \describe{
#'     \item{bin}{A factor column representing the bin intervals in the format "[x, y]".}
#'     \item{bin_size}{A numeric column representing the bin size.}
#'     \item{binned_time}{A numeric column representing the binned time.}
#'   }
#'
#' @return A numeric vector containing four elements:
#'   \describe{
#'     \item{lower_bound}{The lower bound of the bin interval.}
#'     \item{upper_bound}{The upper bound of the bin interval.}
#'     \item{bin_size}{The bin size.}
#'     \item{binned_time}{The binned time.}
#'   }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#'
#' @importFrom stringr str_remove_all
#'
#' @keywords internal
create_range <- function(x) {
  # Convert the factor column "bin" to character
  y <- as.character(x[["bin"]])

  # Remove square and round brackets from the character string
  y <- y %>% str_remove_all(pattern = "\\[|\\]|\\(|\\)")

  # Split the character string by comma and extract the first element (lower bound of the range)
  y1 <- as.numeric(sapply(strsplit(y, ","), "[", 1))

  # Split the character string by comma and extract the second element (upper bound of the range)
  y2 <- as.numeric(sapply(strsplit(y, ","), "[", 2))

  # Combine the lower bound, upper bound, bin size, and binned time into a numeric vector
  rangeVec <- c(y1, y2, x[["bin_size"]], x[["binned_time"]])

  # Return the numeric vector
  return(as.numeric(rangeVec))
}
