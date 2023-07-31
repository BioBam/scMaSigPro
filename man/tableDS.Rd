\name{tableDS}
\alias{tableDS}

\title{Identification of Mayor and minor Isoforms in the clusters}
\description{
  \code{tableDS} identifies for each Differentialy Spliced Gene (DSG) the clusters where their isoforms belong to, labelling gene transcripts as mayor (or most expressed) and minor.
  }
\usage{
 tableDS(seeDS)
}

\arguments{
  \item{seeDS}{a \code{seeDS} object }
    }

\details{
This table includes DSG with 2 or more Isoforms. Mono isoform genes are useful to determine the trends of the cluster. However, as they have only one Isoform, there is not the possibility of comparing minor and major DETs.

}

\value{
  \item{IsoTable}{A classification table that indicates the distribution of isoforms across diferent clusters}
  \item{IsoClusters}{A data.frame with genes in rows and two columns: first indicates the number of cluster of the major isoform and second the number(s) of cluster(s) of the minor isoforms.}
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

\seealso{ \code{\link{seeDS}}, \code{\link{IsoModel}} }
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
