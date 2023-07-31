\name{suma2Venn}
\alias{suma2Venn}

\title{ Creates a Venn Diagram from a matrix of characters}
\description{
\code{suma2Venn} transforms a matrix or a data frame with characters into a list to draw and display a Venn diagram with up to 7 sets
}

\usage{
suma2Venn(x, size = 30, cexil = 0.9, cexsn = 1, zcolor = heat.colors(ncol(x)), ...)
}

\arguments{
  \item{x}{ matrix or data frame of character values}
  \item{size}{ Plot size, in centimeters}
  \item{cexil}{ Character expansion for the intersection labels}
  \item{cexsn}{ Character expansion for the set names}
  \item{zcolor}{A vector of colors for the custom zones}
  \item{\dots}{ Additional plotting arguments for the venn function }
}
\details{
  \code{suma2Venn} creates a list with the columns of a matrix or a data frame of characters which can be taken by the
\code{\link[venn:venn]{venn}} to generate a Venn Diagram
}

\value{
\code{suma2Venn} returns a Venn Plot such as that created by the
\code{\link[venn:venn]{venn}} function

}

\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}

\seealso{ \code{\link[venn:venn]{venn}}}
\examples{
A <- c("a","b","c", "d", "e", NA, NA)
B <- c("a","b","f", NA, NA, NA, NA)
C <- c("a","b","e","f", "h", "i", "j", "k")
x <- cbind(A, B, C)
suma2Venn(x)

}
\keyword{ aplot }

