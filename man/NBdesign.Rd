\name{NBdesign}
\alias{NBdesign}
\docType{data}
\title{Experimental design for RNA-Seq example}
\description{
  \code{NBdesign} contains a subset of a bigger normalized negative binomial simulated dataset. 
}

\usage{data(NBdesign)}
\format{
  A matrix with 36 rows and 4 colums
 
      rows [1:36]  "G1.T1.1" "G1.T1.2" "G1.T1.3" "G1.T2.1" ...
 
      columns [1:6] [1] "Time"       "Replicates" "Group.1"    "Group.2"  
}

\details{
 Samples are given in rows and experiment descriptors are given in columns. Row
 names contain sample names. 

 \code{"Time"} indicates the values that variable Time takes in each
 experimental condition. There are 6 time points.
 

 \code{"Replicates"} is an index indicating the same experimental condition.

  \code{ "Group.1"} and \code{"Group.2"} columns indicate assigment to experimental groups,
  coding with 1 and 0 whether each sample belongs to that group or not. 

}

\examples{
data(NBdesign)

}
\keyword{datasets}
