\name{edesign.abiotic}
\alias{edesign.abiotic}
\docType{data}
\title{Experimental design potato abiotic stress}
\description{
 \code{edesign.abiotic} contains experimental set up of a time course microarray experiment where potato plants were submitted to 3 different abiotic stresses.
}
\usage{data(edesign.abiotic)}
\format{
   A matrix with 36 rows and 6 columns
\describe{
      rows [1:36] "Control 3h 1" "Control 3h 2" "Control 3h 3" "Control
      9h 1" ...

      columns [1:6] "Time" "Replicates" "Control" "Cold" "Heat" "Salt" 
  }
}
\details{
 Arrays are given in rows and experiment descriptors are given in columns. Row
 names contain array names.

    \code{"Time"} indicates the values that variable Time takes in each
    hybridization. 

    \code{"Replicates"} is an index indicating replicate hyridizations,
    i.e. hybridizations are numbered, giving replicates the same number. 

    \code{"Control"}, \code{"Cold"}, \code{"Heat"} and \code{"Salt"} columns
    indicate array assigment to experimental groups, coding with 1 and 0
    whether each array belongs to that group or not. 
}

\references{
Rensink WA, Iobst S, Hart A, Stegalkina S, Liu J, Buell CR. Gene expression 
profiling of potato responses to cold, heat, and salt stress. Funct Integr 
Genomics. 2005 Apr 22.
}
\examples{
data(edesignCR)
}
\keyword{datasets}
