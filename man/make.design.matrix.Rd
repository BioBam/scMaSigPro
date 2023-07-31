\name{make.design.matrix}
\alias{make.design.matrix}

\title{Make a design matrix for regression fit of time series gene expression experiments}
\description{
\code{make.design.matrix} creates the design matrix of dummies for fitting time series micorarray
     gene expression experiments.
}
\usage{
make.design.matrix(edesign, degree = 2, time.col = 1,
                   repl.col = 2, group.cols = c(3:ncol(edesign)))
}

\arguments{
  \item{edesign}{ matrix describing experimental design. Rows must be arrays and columns experiment descriptors}
  \item{degree}{ the degree of the regression fit polynome. \code{degree} = 1 returns linear regression, \code{degree} = 2 returns quadratic regression, etc}
  \item{time.col}{ column number in edesign containing time values. Default is first column }
  \item{repl.col}{ column number in edesign containing coding for replicate arrays. Default is second column}
  \item{group.cols}{column numbers in edesign indicating the coding for each experimental group (treatment, tissue, ...). See details }
}
\details{
     rownames of edesign object should contain the arrays naming (i.e. array1, array2, ...).
     colnames of edesign must contain the names of experiment descriptors(i.e. "Time", "Replicates", "Treatment A", "Treatment B", etc.).
     for each experimental group a different column must be present in edesign, coding with 1 and 0 whether each array belongs to that group or not.

     \code{make.design.matrix} returns a design matrix where rows represent arrays and column variables of time, dummies and their interactions for up to the \code{degree} given.
     Dummies show the relative effect of each experimental group related to the first one. Single dummies indicate the abcissa component of each group.
     $Time*dummy$ variables indicate slope changes, $Time^2*dummy$ indicates curvature changes. Higher grade values could model complex responses.
     In case experimental groups share a initial state (i.e. common time 0), no single dummies are modeled.
}
\value{

  \item{dis }{design matrix of dummies for fitting time series}
  \item{groups.vector }{vector coding the experimental group to which each variable belongs to}
  \item{edesign }{\code{edesign} value passed as argument}

}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
Bioinformatics 22, 1096-1102}
\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}


\examples{
data(edesign.abiotic, edesignCT)
make.design.matrix(edesign.abiotic)  # quadratic model
make.design.matrix(edesignCT, degree = 3)  # cubic model with common starting time point
}

\keyword{ design }
\keyword{ regression }
