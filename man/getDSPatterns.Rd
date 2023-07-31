\name{getDSPatterns}
\alias{getDSPatterns}

\title{Lists of genes with Isoforms in different clusters}
\description{
  \code{getDSPatterns} is a function that makes a list with the names of genes identified with \code{tableDS} function.
  }
\usage{
getDSPatterns(tableDS, Cluster.Major, Cluster.minor)
}

\arguments{
  \item{tableDS}{a \code{tableDS} object }
  \item{Cluster.Major}{Number of the cluster where the major isoform belongs to}
  \item{Cluster.minor}{Number(s) of the cluster(s) where the minor isoform(s) belongs to (see details)}
  }

\details{
When minor isoforms belong to different clusters, \code{tableDS} codifies them using "&". For instance: clusters 1 and 2,
will be represented as "1&2". In such cases quotation marks must be used (see examples).
When minor isoforms are only in one cluster there is no need to use quotation marks.
}

\value{
 A vector with the names of the genes.
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

\seealso{ \code{\link{tableDS}}, \code{\link{IsoModel}} }

\examples{
data(ISOdata)
data(ISOdesign)
mdis <- make.design.matrix(ISOdesign)
MyIso <- IsoModel(data=ISOdata[,-1], gen=ISOdata[,1], design=mdis, counts=TRUE)
Myget <- getDS(MyIso)
see <- seeDS(Myget, cluster.all=FALSE, k=6)
table <- tableDS(see)
table$IsoTable

getDSPatterns(table, 1, 4)
getDSPatterns(table, "1", "4") #will give the same result.

getDSPatterns(table, 1, "1&5")
}
