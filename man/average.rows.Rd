\name{average.rows}
\alias{average.rows}

\title{ Average rows by match and index}
\description{
     \code{average.rows} matches rownames of a matrix to a \code{match} vector and performs averaging of the rows
      by the index provided by an \code{index} vector.
}
\usage{
average.rows(x, index, match, r = 0.7)
}

\arguments{
  \item{x}{ a matrix}
  \item{index}{ index vector indicating how rows must be averaged }
  \item{match}{ match vector for indexing rows }
  \item{r}{minimal correlation value between rows to compute average }
}
\details{
    rows will be averaged only if the pearson correlation coefficient between all rows of each given
    index is greater than r. If not, that group of rows is discarded in the result matrix.
}
\value{
    a matrix of averaged rows
}

\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}

\examples{

## create data matrix for row averaging
x <- matrix(rnorm(30), nrow = 6,ncol = 5)
rownames(x) <- paste("ID", c(1, 2, 11, 12, 19, 20), sep = "")
i <- paste("g", rep(c(1:10), each = 2), sep = "")  # index vector
m <- paste("ID", c(1:20), sep = "")  # match vector
average.rows(x, i, m, r = 0)

}
\keyword{ arith}
