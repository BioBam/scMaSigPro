\name{IsoPlot}
\alias{IsoPlot}

\title{Plotting the isoform profiles of a specific gene by groups}
\description{
This function makes a plot with the isoforms of a specific gene splitting the different experimental groups.
}
\usage{
IsoPlot(get, name, only.sig.iso=FALSE, ylim=NULL, xlab = "Time",
ylab = "Expression value", points=TRUE, cex.main=3,cex.legend=1.5)
}

\arguments{
  \item{get}{a \code{getDS} object a cluster of flat Isoform}
  \item{name}{Name of the specific gen to show in the plot}
  \item{only.sig.iso}{TRUE when the plot is made only with statistically significant isoforms.}
  \item{ylim}{Range of the y axis of the desired plot. If it is NULL it will be computed automatically. }
  \item{xlab}{label for the x axis}
  \item{ylab}{label for the y axis}
  \item{points}{ TRUE to plot points and lines. FALSE to plot only lines. }
  \item{cex.main}{  graphical parameter magnification to be used for main}
  \item{cex.legend}{ graphical parameter magnification to be used for legend }
 }

\details{
 The plot can be made with all the available isoforms or only with the statistilly significant ones.
}
\value{
    Plot of isoform profiles of a specific gene by groups.
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

\seealso{ \code{\link{getDS}}, \code{\link{IsoModel}} }


\examples{
data(ISOdata)
data(ISOdesign)
mdis <- make.design.matrix(ISOdesign)
MyIso <- IsoModel(data=ISOdata[,-1], gen=ISOdata[,1], design=mdis, counts=TRUE)
Myget <- getDS(MyIso)

IsoPlot(Myget,"Gene1005",only.sig.iso=FALSE,cex.main=2,cex.legend=1)

}
