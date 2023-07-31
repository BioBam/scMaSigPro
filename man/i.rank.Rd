\name{i.rank}
\alias{i.rank}

\title{Ranks a vector to index}
\description{
  Ranks the values in a vector to sucessive values. Ties are given the same value.
}
\usage{
i.rank(x)
}

\arguments{
  \item{x}{vector}
}

\value{Vector of ranked values
}

\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}

\seealso{\code{\link{rank}},\code{\link{order}}}
\examples{
     i.rank(c(1, 1, 1, 3, 3, 5, 7, 7, 7))

}
\keyword{ arith }
