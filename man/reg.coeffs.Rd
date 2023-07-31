\name{reg.coeffs}
\alias{reg.coeffs}

\title{ Calculate true variables regression coefficients}
\description{
  \code{reg.coeffs} calculates back regression coefficients for true variables (experimental groups) from dummy variables regression coefficients.
}
\usage{
reg.coeffs(coefficients,
indepen = groups.vector[nchar(groups.vector)==min(nchar(groups.vector))][1],
groups.vector,  group)
}

\arguments{
  \item{coefficients}{ vector of regression coefficients obtained from a regression model with dummy variables }
  \item{indepen}{ idependent variable of the regression formula}
  \item{groups.vector}{ vector indicating the true variable of each variable in \code{coefficients}  }
  \item{group}{ true variable for which regression coefficients are to be computed }
}
\details{
     regression coefficients in coefficients vector should be ordered by polynomial degree in a regression formula,
     ie: intercept, $x$ term, $x^2$ term, $x^3$ term, and so on...
}
\value{

  \item{reg.coeff}{vector of calculated regression coefficients}

}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2005.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
 }
\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}

\examples{
groups.vector <-c("CT", "T1vsCT", "T2vsCT", "CT", "T1vsCT","T2vsCT", "CT", "T1vsCT", "T2vsCT")
coefficients <- c(0.1, 1.2, -0.8, 1.7, 3.3, 0.4, 0.0, 2.1, -0.9)
## calculate true regression coefficients for variable "T1"
reg.coeffs(coefficients, groups.vector = groups.vector, group = "T1")
}

\keyword{ misc }
