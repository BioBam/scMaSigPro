% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_functions.R
\name{showSol}
\alias{showSol}
\title{Show or Return the P-values after model fitting.}
\usage{
showSol(scmpObj, view = FALSE, return = TRUE, includeInflu = TRUE)
}
\arguments{
\item{scmpObj}{An object of class \code{\link{ScMaSigPro}}.}

\item{view}{Whether to view the data in the explorer. (Default: FALSE)}

\item{return}{Whether to return the data. (Default: TRUE)}

\item{includeInflu}{Whether to include genes with influential observations.}
}
\value{
The computed p-values for each term and full model as a dataframe.
}
\description{
This function is used to view or return the matrix of p-values for each term
and the full model from the provided scMaSigPro object.
}
\author{
Priyansh Srivastava \email{spriyansh29@gmail.com}
}
