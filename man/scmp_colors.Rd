% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{scmp_colors}
\alias{scmp_colors}
\title{Get Color Palette}
\usage{
scmp_colors(n)
}
\arguments{
\item{n}{An integer specifying the number of colors to return.}
}
\value{
A character vector of hex color codes.
}
\description{
This function returns a specified number of contrasting colors from a predefined palette.
If the requested number of colors exceeds the length of the predefined palette,
additional unique colors from the default R color set are included.
}
\keyword{internal}
