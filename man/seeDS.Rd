\name{seeDS}
\alias{seeDS}

\title{Wrapper function for visualization of significant isoforms from Differentially Spliced Genes}
\description{
  \code{seeDS}
     This function provides visualisation tools for Significant Isoforms in a time course
     experiment. The function calls the \code{see.genes} function for selected Isoforms. This cluster will be the reference in \code{tableDS} function to identify the trends that follows the isoforms of a specific gene.}
\usage{
 seeDS(get, rsq=0.4, cluster.all=TRUE, plot.mDSG=FALSE, k=6,
 cluster.method="hclust", k.mclust=FALSE, ...)
}

\arguments{
  \item{get}{a \code{getDS} object a cluster of flat Isoform}
  \item{rsq}{Required when cluster.all=TRUE. It is the cut-off level at the R-squared value for detecting significant isoforms of all the genome.}
  \item{cluster.all}{TRUE to make the cluster with significant isoforms of all the genome.
  FALSE to make the cluster with significant isoforms of Differentially Spliced Genes.}
\item{plot.mDSG}{TRUE to make a cluster with the Isoforms belonging to monoIsoform genes}
  \item{k}{number of clusters for data partioning }
 \item{cluster.method}{clustering method for data partioning. Currently \code{"hclust"}, \code{"kmeans"} and  \code{"Mclust"} are supported}
   \item{k.mclust}{TRUE for computing the optimal number of clusters with Mclust algorithm}
  \item{\dots}{other graphical function argument }
  }

\details{The cluster reference can be made with significant isoforms of all the genome or with the isoforms belonging to the Differentially Spliced Genes.

     Alternatively a cluster of monoIsoforms can be asked.

     Next a partioning of the data is generated using a clustering method.

     The results of the clustering are visualized in two plots as in \code{see.genes}.
}
\value{
	Experiment wide Isoform profiles and by group profiles plots are generated for each data cluster in the graphical device.
  \item{Model}{a \code{IsoModel} object to be used in the following steps}
  \item{get}{a \code{get.siggenes} object to be used in the following steps}
  \item{NumIso.by.gene}{Number of selected Isoforms for each Differentially Spliced Gene}
  \item{cut}{vector indicating gene partioning into clusters}
   \item{names.genes}{vector with the name of the gene each selected isoform belongs to}
}
\references{
  Nueda, M.J., Martorell, J., Marti, C., Tarazona, S., Conesa, A. 2018. Identification and visualization of differential isoform expression in RNA-seq time series. Bioinformatics. 34, 3, 524-526.

  Nueda, M.J., Tarazona, S., Conesa, A. 2014.
  Next maSigPro: updating maSigPro bioconductor package for RNA-seq time series. Bioinformatics, 30, 2598-602.

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
MyIso <- IsoModel(data=ISOdata[,-1], gen=ISOdata[,1], design=mdis, counts=TRUE)
Myget <- getDS(MyIso)
see <- seeDS(Myget, cluster.all=FALSE, k=6)

table <- tableDS(see)
table$IsoTable
}
