\name{ISOdesign}
\alias{ISOdesign}
\docType{data}
\title{Experimental design for ISOdata dataset example}
\description{
  \code{ISOdesign} is the experimental design to apply ISOmaSigPro to ISOdata dataset example.
 
}

\usage{data(ISOdesign)}
\format{
  A matrix with 36 rows and 4 colums
 
rownames(ISOdesign)
 "Gr1_0h_1"  "Gr1_0h_2"  "Gr1_0h_3"  "Gr1_2h_1"  "Gr1_2h_2"  "Gr1_2h_3"  "Gr1_6h_1"  "Gr1_6h_2"  "Gr1_6h_3" 
 "Gr1_12h_1" "Gr1_12h_2" "Gr1_12h_3" "Gr1_18h_1" "Gr1_18h_2" "Gr1_18h_3" "Gr1_24h_1" "Gr1_24h_2" "Gr1_24h_3"     "Gr2_0h_1"  "Gr2_0h_2"  "Gr2_0h_3"  "Gr2_2h_1"  "Gr2_2h_2"  "Gr2_2h_3"  "Gr2_6h_1"  "Gr2_6h_2"  "Gr2_6h_3" 
 "Gr2_12h_1" "Gr2_12h_2" "Gr2_12h_3" "Gr2_18h_1" "Gr2_18h_2" "Gr2_18h_3" "Gr2_24h_1" "Gr2_24h_2" "Gr2_24h_3"
 
 colnames(ISOdesign)
 "time"      "replicate" "Group1"    "Group2"   
 
}

\details{
 Samples are given in rows and experiment descriptors are given in columns. Row
 names contain sample names.
 

 \code{"time"} indicates the values that variable Time takes in each
 experimental condition. There are 6 time points.

 
 \code{"replicate"} is an index indicating the same experimental condition.

  \code{ "Group1"} and \code{"Group2"} columns indicate assigment to experimental groups,
  coding with 1 and 0 whether each sample belongs to that group or not.
 

}

\examples{
data(ISOdata)
data(ISOdesign)

mdis <- make.design.matrix(ISOdesign)
MyIso <- IsoModel(data=ISOdata[,-1], gen=ISOdata[,1], design=mdis, counts=TRUE)

}
\keyword{datasets}
