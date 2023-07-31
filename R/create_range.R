#' Create Range Function
#'
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
#' @examples
#' # Sample data frame
#' data <- data.frame(
#'   bin = factor(c("[1, 5]", "[6, 10]", "[11, 15]")),
#'   bin_size = c(5, 5, 5),
#'   binned_time = c(10, 15, 20)
#' )
#' 
#' # Calculate the range using the create_range function
#' result <- create_range(data)
#' print(result)
#'
#' @keywords internal
#' @export
#' 
create_range <- function(x) {
    # Convert the factor into Character
    y <- as.character(x[["bin"]])
    y <- y %>% str_remove_all(pattern = "\\[|\\]|\\(|\\)")
    y1 <- as.numeric(sapply(strsplit(y, ","), "[", 1))
    y2 <- as.numeric(sapply(strsplit(y, ","), "[", 2))
    rangeVec <- c(y1, y2, x[["bin_size"]], x[["binned_time"]])
    return(as.numeric(rangeVec))
}
