\name{PlotProfiles}
\alias{PlotProfiles}

\title{
Function for visualization of gene expression profiles
}
\description{
  \code{PlotProfiles} displays the expression profiles of a group of genes.
}
\usage{
PlotProfiles(data, cond, cex.axis = 0.5, ylim = NULL, repvect,
main = NULL, sub = NULL, color.mode = "rainbow", item = NULL)
}

\arguments{
  \item{data}{ a matrix containing the gene expression data}
  \item{cond}{ vector for x axis labeling, typically array names }
  \item{cex.axis}{ graphical parameter maginfication to be used for x axis in plotting functions }
  \item{ylim}{ index vector indicating experimental replicates }
  \item{repvect}{ index vector indicating experimental replicates}
  \item{main}{ plot main title}
  \item{sub}{ plot subtitle  }
  \item{color.mode}{ color scale for plotting profiles. Can be either \code{"rainblow"} or \code{"gray"}}
  \item{item}{Name of the analysed items to show}
}
\details{
The \code{repvect} argument is used to indicate with vertical lines groups of replicated arrays.
}
\value{Plot of experiment-wide gene expression profiles.
}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2005.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
 }
\author{ Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}

\seealso{ \code{\link{PlotGroups}} }
\examples{
#### GENERATE TIME COURSE DATA
## generate n random gene expression profiles of a data set with
## one control plus 3 treatments, 3 time points and r replicates per time point.

tc.GENE <- function(n, r,
             var11 = 0.01, var12 = 0.01,var13 = 0.01,
             var21 = 0.01, var22 = 0.01, var23 =0.01,
             var31 = 0.01, var32 = 0.01, var33 = 0.01,
             var41 = 0.01, var42 = 0.01, var43 = 0.01,
             a1 = 0, a2 = 0, a3 = 0, a4 = 0,
             b1 = 0, b2 = 0, b3 = 0, b4 = 0,
             c1 = 0, c2 = 0, c3 = 0, c4 = 0)
{

  tc.dat <- NULL
  for (i in 1:n) {
    Ctl <- c(rnorm(r, a1, var11), rnorm(r, b1, var12), rnorm(r, c1, var13))  # Ctl group
    Tr1 <- c(rnorm(r, a2, var21), rnorm(r, b2, var22), rnorm(r, c2, var23))  # Tr1 group
    Tr2 <- c(rnorm(r, a3, var31), rnorm(r, b3, var32), rnorm(r, c3, var33))  # Tr2 group
    Tr3 <- c(rnorm(r, a4, var41), rnorm(r, b4, var42), rnorm(r, c4, var43))  # Tr3 group
    gene <- c(Ctl, Tr1, Tr2, Tr3)
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat
}

## create 10 genes with profile differences between Ctl, Tr2, and Tr3 groups
tc.DATA <- tc.GENE(n = 10,r = 3, b3 = 0.8, c3 = -1, a4 = -0.1, b4 = -0.8, c4 = -1.2)
rownames(tc.DATA) <- paste("gene", c(1:10), sep = "")
colnames(tc.DATA) <- paste("Array", c(1:36), sep = "")

PlotProfiles (tc.DATA, cond = colnames(tc.DATA), main = "Time Course",
              repvect = rep(c(1:12), each = 3))

}
\keyword{ aplot }
