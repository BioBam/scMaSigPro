\name{PodiumChange}
\alias{PodiumChange}

\title{Detection of Genes with switchs of their major isoforms}
\description{
     This function provides lists of genes that have different Major isoforms (most expressed) when different intervals of the experimental conditions are considered.

     The subrange of the experimental conditions can be chosen as a specific point, all the points of a specific experimental group or at any point.
}
\usage{
PodiumChange(get, only.sig.iso=FALSE,  comparison=c("any",
"groups","specific"), group.name="Ctr", time.points=0)
}

\arguments{
  \item{get}{a \code{getDS} object a cluster of flat Isoform}
  \item{only.sig.iso}{TRUE when changes are looked for only through statistically significant isoforms.}
  \item{comparison}{Type of search to do: any, groups or specific (see details).}
  \item{group.name}{required when \code{comparison} is \code{"specific"}.}
  \item{time.points}{required when \code{comparison} is \code{"specific"}.}
 }

\details{
    There are 3 possible values for the comparison argument:

    \code{"any"}: Detects genes with Major Isoform changes in at least one experimental condition.

    \code{"groups"}: Detects genes with different Major Isoform for different experimental groups.

    \code{"specific"}: Detects genes with Major Isoform changes in a specific time interval, especified in       \code{time.points} argument and a specific experimental group, especified in \code{group.name} argument.
}
\value{
  \item{L}{Names of the genes with PodiumChange Isoforms }
  \item{data.L}{Data values of all the isoforms belonging to the genes in L}
  \item{gen.L}{gen vector with the name of the gene of each isoform}
  \item{edesign}{matrix describing experimental design needed to visualize PodiumChange selection with \code{IsoPlot} function. It is the input of \code{make.design.matrix}.}
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

\seealso{ \code{\link{see.genes}}, \code{\link{IsoModel}} }


\examples{
data(ISOdata)
data(ISOdesign)
mdis <- make.design.matrix(ISOdesign)
MyIso <- IsoModel(data=ISOdata[,-1], gen=ISOdata[,1], design=mdis,
counts=TRUE)
Myget <- getDS(MyIso)

PC <- PodiumChange(Myget, only.sig.iso=TRUE,
comparison="specific", group.name="Group2", time.points=c(18,24))
PC$L

}
