% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScMaSigProClass.R
\docType{class}
\name{ScMaSigPro}
\alias{ScMaSigPro}
\alias{ScMaSigPro-class}
\title{ScMaSigPro}
\description{
S4 class to store analysis of ScMaSigPro worflow. It stores results,
parameters and associated data. Inherits slots from \code{SingleCellExperiment}.
}
\section{Slots}{

\describe{
\item{\code{Sparse}}{Object of class \code{\link{SingleCellExperiment}}.
See \pkg{SingleCellExperiment} for more details.}

\item{\code{Dense}}{Object of class \code{\link{SingleCellExperiment}}.
See \pkg{SingleCellExperiment} for more details.}

\item{\code{Design}}{Object of Class \code{\link{MatrixDesign}}.}

\item{\code{Profile}}{Object of Class \code{\link{VariableProfiles}}.}

\item{\code{Estimate}}{Object of Class \code{\link{Estimates}}.}

\item{\code{Significant}}{Object of Class \code{\link{Significant}}.}

\item{\code{Parameters}}{Object of Class \code{\link{ParameterConfig}}.}
}}

\keyword{classes}
