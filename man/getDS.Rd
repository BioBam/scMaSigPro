\name{getDS}
\alias{getDS}

\title{Extract lists of significant isoforms from Differentially Spliced Genes (DSG)}
\description{
  \code{getDS} creates lists of significant isoforms from Differentially Spliced Genes (DSG)
  }
\usage{
  getDS(Model, vars="all", rsq=0.4)
}

\arguments{
  \item{Model}{a \code{IsoModel} object }
  \item{vars}{argument of the \code{\link{get.siggenes}} function applied to isoforms}
  \item{rsq}{cut-off level at the R-squared value for the stepwise regression fit. Only isoforms with R-squared more than rsq are selected}
  }

\details{
There are 3 possible values for the vars argument: "all", "each" and "groups". See \code{\link{get.siggenes}}.

}
\value{
In the console a summary of the selection is printed.
  \item{Model}{a \code{IsoModel} object to be used in the following steps}
  \item{get2}{a \code{get.siggenes} object to be used in the following steps}
  \item{DSG}{Names of the selected genes: Differentially Spliced Genes}
  \item{DET}{Names of the selected Isoforms: Differentally Expressed Transcripts}
  \item{List0}{a list with the names of Differentially Spliced Genes without Isoforms with R-squared higher than rsq}
  \item{NumIso.by.gene}{Number of selected Isoforms for each Differentially Spliced Gene}
}
\references{
  Nueda, M.J., Martorell, J., Marti, C., Tarazona, S., Conesa, A. 2018. Identification and visualization of differential isoform expression in RNA-seq time series. Bioinformatics. 34, 3, 524-526.

  Nueda, M.J., Tarazona, S., Conesa, A. 2014.
  Next maSigPro: updating maSigPro bioconductor package for RNA-seq time series.
  Bioinformatics, 30, 2598-602.

  Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
  maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
  Bioinformatics 22, 1096-1102.

}
\author{Maria Jose Nueda, \email{mj.nueda@ua.es}}

\seealso{ \code{\link{get.siggenes}}, \code{\link{IsoModel}} }
\examples{
data(ISOdata)
data(ISOdesign)
mdis <- make.design.matrix(ISOdesign)
MyIso <- IsoModel(data=ISOdata[,-1], gen=ISOdata[,1], design=mdis, counts=TRUE)

Myget <- getDS(MyIso)
Myget$DSG
Myget$DET

see <- seeDS(Myget, cluster.all=FALSE, k=6)
table <- tableDS(see)
table$IsoTable
}
