\name{edesignCT}
\alias{edesignCT}
\docType{data}
\title{Experimental design with a shared time}
\description{
  \code{edesignCT} contains the experimental set up of a time course
  microarray experiment where there is a common starting point for the
  different experimental groups. 
}

\usage{data(edesignCT)}
\format{
  A matrix with 32 rows and 7 colums
 
      rows [1:32] "Array1" "Array2" "Array3" "Array4" ...
 
      columns [1:7] "Time" "Replicates" "Control" "Tissue1" "Tissue2"
      "Tissue3" "Tissue4" 
}

\details{
 Arrays are given in rows and experiment descriptors are given in columns. Row
 names contain array names. 

 \code{"Time"} indicates the values that variable Time takes in each
 hybridization. There are 4 time points, which allows an up to 3 degree
 regression polynome. 

 \code{"Replicates"} is an index indicating replicate hyridizations, i.e. hybridizations are numbered, giving replicates the same number.

  \code{ "Control"}, \code{"Tissue1"}, \code{"Tissue2"}, \code{"Tissue3"} and
  \code{"Tissue4"}  columns indicate array assigment to experimental groups,
  coding with 1 and 0 whether each array belongs to that group or not. 

}

\examples{
data(edesignCT)

}
\keyword{datasets}
