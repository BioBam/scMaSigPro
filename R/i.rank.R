#' Ranks a vector to index
#'
#' \code{i.rank} ranks the values in a vector to successive values. Ties are given the same value.
#'
#' @title Ranks a vector to index
#'
#' @description
#' \code{i.rank} ranks the values in a vector to successive values. Ties are given the same value.
#'
#' @usage
#' i.rank(x)
#'
#' @arguments
#' \item{x}{vector} The input vector that will be ranked.
#'
#' @value
#' Vector of ranked values.
#'
#' @author Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}
#'
#' @seealso \code{\link{rank}},\code{\link{order}} Functions related to ranking and ordering vectors.
#'
#' @examples
#' i.rank(c(1, 1, 1, 3, 3, 5, 7, 7, 7))
#'
#' @keyword arith

i.rank <- function(x) {
  xx <- x
  for (i in 1:length(xx)) {
    xx[i] <- c(1:length(unique(x)))[x[i] == unique(x)]
  }
  xx
}
