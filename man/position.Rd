\name{position}
\alias{position}

\title{ Column position of a variable in a data frame }
\description{
  Finds the column position of a character variable in the column names of a data frame.
}
\usage{
position(matrix, vari)
}

\arguments{
  \item{matrix}{ matrix or data.frame with character column names }
  \item{vari}{ character variable}
}

\value{
numerical. Column position for the given variable.

}

\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}

\examples{
x <- matrix(c(1, 1, 2, 2, 3, 3),ncol = 3,nrow = 2)
colnames(x) <- c("one", "two", "three")
position(x, "one")

}
\keyword{ arith}
