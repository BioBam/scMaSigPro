\name{ISOdata}
\alias{ISOdata}
\docType{data}
\title{ RNA-Seq dataset example for isoforms}
\description{

 \code{ISOdata} contains an example of RNA-Seq data at Isoform level.

}
\usage{data(ISOdata)}
\format{
  A data frame with 2782 rows and 37 columns with RNA-Seq data.
}
\details{
  Rows correspond to 2782 isoforms belonging to 1000 gene.
  
  First column is the name of the gene each isoform belongs to.
  
  Remaining columns are the RNA-Seq data samples asociated to 3 replicates of 12 experimental conditions.
}


\examples{
data(ISOdata)
data(ISOdesign)

mdis <- make.design.matrix(ISOdesign)
MyIso <- IsoModel(data=ISOdata[,-1], gen=ISOdata[,1], design=mdis, counts=TRUE)
}
\keyword{datasets}
