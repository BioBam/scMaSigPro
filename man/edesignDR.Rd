\name{edesignDR}
\alias{edesignDR}
\docType{data}
\title{Experimental design with different replicates}
\description{
\code{edesignDR} contains experimental set up of a replicated time course microarray experiment where rats were submitted to 3 different dosis of a toxic compound. A control and an placebo treatments are also present in the experiment.
}
\usage{data(edesignDR)}
\format{
  A matrix with 54 rows and 7 columns

      rows [1:54] "Array1" "Array2" "Array3" "Array4" ...

      columns [1:7] "Time" "Replicates" "Control" "Placebo" "Low" "Medium" "High" 

}

\details{
 Arrays are given in rows and experiment descriptors are given in columns. Row
 names contain array names. 

\code{"Time"} indicates the values that variable Time takes in each hybridization.

\code{"Replicates"} is an index indicating replicate hyridizations, i.e. hybridizations are numbered, giving replicates the same number.

\code{"Control"}, \code{"Placebo"}, \code{"Low"}, \code{"Medium"} and \code{"High"} columns indicate array assigment to experimental groups, coding with 1 and 0 whether each array belongs to that group or not.
}

\references{
Heijne, W.H.M.; Stierum, R.; Slijper, M.; van Bladeren P.J. and van Ommen B.(2003).
Toxicogenomics of bromobenzene hepatotoxicity: a combined transcriptomics and proteomics approach.
Biochemical Pharmacology 65 857-875.
}
\examples{
data(edesignDR)

}
\keyword{datasets}
