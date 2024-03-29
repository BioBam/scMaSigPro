% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scmpGenerics.R
\name{eSparse}
\alias{eSparse}
\title{Get Expression Counts from Sparse Slot}
\usage{
eSparse(object, i = "counts", ...)
}
\arguments{
\item{object}{An object of class \code{\link{ScMaSigPro}}.}

\item{i}{The name of the assay to retrieve expression counts for.
Default is "counts", which will return the default expression matrix.}

\item{...}{Additional arguments (not used currently).}
}
\value{
An expression matrix for the specified assay from the Sparse slot.
}
\description{
`eSparse` is a generic function for retrieving expression counts stored in the
Sparse slot of an \code{\link{ScMaSigPro}} object.
}
\author{
Priyansh Srivastava \email{spriyansh29@gmail.com}
}
