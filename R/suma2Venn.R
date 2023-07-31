#' Creates a Venn Diagram from a matrix of characters
#'
#' \code{suma2Venn} transforms a matrix or a data frame with characters into a list to draw and display a Venn diagram with up to 7 sets.
#'
#' @title Creates a Venn Diagram from a matrix of characters
#'
#' @description
#' \code{suma2Venn} transforms a matrix or a data frame with characters into a list to draw and display a Venn diagram with up to 7 sets.
#'
#' @usage
#' suma2Venn(x, size = 30, cexil = 0.9, cexsn = 1, zcolor = heat.colors(ncol(x)), ...)
#'
#' @arguments
#' \item{x}{matrix or data frame of character values}
#' \item{size}{Plot size, in centimeters}
#' \item{cexil}{Character expansion for the intersection labels}
#' \item{cexsn}{Character expansion for the set names}
#' \item{zcolor}{A vector of colors for the custom zones}
#' \item{\dots}{Additional plotting arguments for the venn function}
#'
#' @details
#' \code{suma2Venn} creates a list with the columns of a matrix or a data frame of characters which can be taken by the \code{\link[venn:venn]{venn}} function to generate a Venn Diagram.
#'
#' @value
#' \code{suma2Venn} returns a Venn Plot such as that created by the \code{\link[venn:venn]{venn}} function.
#'
#' @author Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}
#'
#' @seealso \code{\link[venn:venn]{venn}}
#'
#' @examples
#' A <- c("a", "b", "c", "d", "e", NA, NA)
#' B <- c("a", "b", "f", NA, NA, NA, NA)
#' C <- c("a", "b", "e", "f", "h", "i", "j", "k")
#' x <- cbind(A, B, C)
#' suma2Venn(x)
#'
#' @keyword aplot

# Define the function "suma2Venn" with arguments x, size, cexil, cexsn, zcolor, and ...
"suma2Venn" <- function(x, size = 30, cexil = 0.9, cexsn = 1, zcolor = heat.colors(ncol(x)), ...) {
  # Get the number of columns in the input matrix or data frame
  G <- ncol(x)
  # Create an empty list "L" to store the columns of the matrix or data frame
  L <- vector("list", G)
  # Set the names of the list "L" to be the column names of the matrix or data frame
  names(L) <- colnames(x)
  # Loop through each column of the matrix or data frame
  for (i in 1:G)
  {
    # Get the values of the current column and convert them to character
    y <- as.character(x[, i])
    # Remove any empty character values from the column
    y <- y[y != " "]
    # Store the processed column values in the list "L"
    L[[i]] <- y
  }

  # Call the "venn" function to draw the Venn diagram using the list "L" and other provided arguments
  venn(L, size = size, cexil = cexil, cexsn = cexsn, zcolor = zcolor, ...)
}
