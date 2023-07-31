\name{PlotGroups}
\alias{PlotGroups}

\title{ Function for plotting gene expression profile at different experimental groups}
\description{
  This function displays the gene expression profile for each experimental group in a time series
     gene expression experiment.
}
\usage{
PlotGroups(data, edesign = NULL, time = edesign[, 1],
groups = edesign[,c(3:ncol(edesign))], repvect = edesign[, 2],
show.lines = TRUE, show.fit = FALSE, dis = NULL,
step.method = "backward", min.obs = 2, alfa = 0.05,
nvar.correction = FALSE, summary.mode = "median",
groups.vector = NULL, main = NULL, sub = NULL, xlab = "Time",
ylab = "Expression value", item = NULL, ylim = NULL, pch = 21,
col = NULL, legend = TRUE, cex.legend = 1,lty.legend = NULL,... )
}

\arguments{
  \item{data}{ vector or matrix containing the gene expression data }
  \item{edesign}{ matrix describing experimental design. Rows must be arrays and columns experiment descriptors}
  \item{time}{ vector indicating time assigment for each array }
  \item{groups}{ matrix indicating experimental group to which each array is assigned }
  \item{repvect}{ index vector indicating experimental replicates}
  \item{show.lines}{ logical indicating whether a line must be drawn joining plotted data points for reach group}
  \item{show.fit}{ logical indicating whether regression fit curves must be plotted}
  \item{dis}{ regression design matrix  }
  \item{step.method}{stepwise regression method to fit models for cluster mean profiles. It can be either \code{"backward"}, \code{"forward"}, \code{"two.ways.backward"} or \code{"two.ways.forward"} }
  \item{min.obs}{minimal number of observations for a gene to be included in the analysis}
  \item{alfa}{ significance level used for variable selection in the stepwise regression }
  \item{nvar.correction}{argument for correcting stepwise regression significance level. See \code{\link{T.fit}} }
  \item{summary.mode}{ the method to condensate expression information when more than one gene is present in the data. Possible values are \code{"representative"} and \code{"median"} }
\item{groups.vector}{ vector indicating experimental group to which each variable belongs }
\item{main}{plot main title }
\item{sub}{ plot subtitle }
\item{xlab}{ label for the x axis }
\item{ylab}{ label for the y axis }
\item{item}{ name of the analysed items to show}
\item{ylim}{ range of the y axis }
\item{pch}{ integer specifying type of points to plot}
\item{col}{ a vector specifying colours to plot. If missing first naturals will be used}
\item{legend}{ logical indicating whether legend must be added when plotting profiles}
\item{cex.legend}{Expansion factor for legend}
\item{lty.legend}{To add a coloured line in the legend}
\item{\dots}{other graphical function argument }
}
\details{
     To compute experimental groups either a edesign object must be provided,
     or separate values must be given for the \code{time}, \code{repvect} and
     \code{groups} arguments.


     When data is a matrix, the average expression value is displayed.

     When there are array replicates in the data (as indicated by
     \code{repvect}), values are averaged by \code{repvect}.

     PlotGroups plots one single expression profile for each experimental
     group even if there are more that one genes in the data set. The way
     data is condensated for this is given by \code{summary.mode}. When this
     argument takes the value \code{"representative"}, the gene with the
     lowest distance to all genes in the cluster will be plotted. When the
     argument is  \code{"median"}, then median expression value is computed.


     When \code{show.fit} is \code{TRUE} the stepwise regression fit for the
     data will be computed and the regression curves will be displayed.

     If data is a matrix of genes and \code{summary.mode} is \code{"median"},
     the regression fit will be computed for the median expression value.

}
\value{
  Plot of gene expression profiles by-group.
}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2005.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
 }
\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}

\seealso{\code{\link{PlotProfiles}}}

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

#### CREATE EXPERIMENTAL DESIGN
Time <- rep(c(rep(c(1:3), each = 3)), 4)
Replicates <- rep(c(1:12), each = 3)
Ctl <- c(rep(1, 9), rep(0, 27))
Tr1 <- c(rep(0, 9), rep(1, 9), rep(0, 18))
Tr2 <- c(rep(0, 18), rep(1, 9), rep(0, 9))
Tr3 <- c(rep(0, 27), rep(1, 9))

PlotGroups (tc.DATA, time = Time, repvect = Replicates, groups = cbind(Ctl, Tr1, Tr2, Tr3))

}
\keyword{ aplot}
