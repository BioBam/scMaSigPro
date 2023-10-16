\name{see.genes}
\alias{see.genes}

\title{Wrapper function for visualization of gene expression values of time course experiments}

\description{
     This function provides visualisation tools for gene expression values in a time course 
     experiment. The function first calls the heatmap function for a general overview of 
     experiment results. Next a partioning of the data is generated using a clustering method.
     The results of the clustering are visualized both as gene expression profiles extended along all
     arrays in the experiment, as provided by the plot.profiles function, and as summary expression 
     profiles for  comparison among experimental groups.
}
\usage{
see.genes(data, edesign = data$edesign, time.col = 1, repl.col = 2, 
    group.cols = c(3:ncol(edesign)), names.groups = colnames(edesign)[3:ncol(edesign)], 
    cluster.data = 1, groups.vector = data$groups.vector, k = 9, k.mclust=FALSE,   
    cluster.method = "hclust", distance = "cor", agglo.method = "ward.D",
    show.lines = TRUE, show.fit = FALSE, dis = NULL, step.method = "backward", 
    min.obs = 3, alfa = 0.05, nvar.correction = FALSE, iter.max = 500, 
    summary.mode = "median", color.mode = "rainbow", ylim = NULL, item = "genes", 
    legend = TRUE, cex.legend = 1, lty.legend = NULL,...) 
}

\arguments{
  \item{data}{ either matrix or a list containing the gene expression data, typically a \code{\link{get.siggenes}} object}
  \item{edesign}{ matrix of experimental design }
  \item{time.col}{  column in edesign containing time values. Default is first column}
  \item{repl.col}{ column in edesign containing coding for replicates arrays. Default is second column}
  \item{group.cols}{ columns indicating the coding for each group (treatment, tissue,...) in the experiment (see details) }
  \item{names.groups}{names for experimental groups}
  \item{cluster.data}{type of data used by the cluster algorithm (see details) }
  \item{groups.vector}{vector indicating the experimental group to which each variable belongs }
  \item{k}{number of clusters for data partioning }
  \item{k.mclust}{TRUE for computing the optimal number of clusters with Mclust algorithm}
  \item{cluster.method}{clustering method for data partioning. Currently \code{"hclust"}, \code{"kmeans"} and  \code{"Mclust"} are supported}
  \item{distance}{distance measurement function when \code{cluster.method} is \code{\link{hclust}}}
  \item{agglo.method}{aggregation method used when \code{cluster.method} is \code{\link{hclust}} }
  \item{show.lines}{logical indicating whether a line must be drawn joining plotted data points for reach group}
  \item{show.fit}{logical indicating whether regression fit curves must be plotted }
  \item{dis}{regression design matrix }
  \item{step.method}{stepwise regression method to fit models for cluster mean profiles. Can be either \code{"backward"}, \code{"forward"}, \code{"two.ways.backward"} or \code{"two.ways.forward"} }
  \item{min.obs}{minimal number of observations for a gene to be included in the analysis}
  \item{alfa}{significance level used for variable selection in the stepwise regression }
  \item{nvar.correction}{argument for correcting \code{T.fit}significance level. See \code{\link{T.fit}} }
  \item{iter.max}{maximum number of iterations when \code{cluster.method} is \code{\link{kmeans}} }
  \item{summary.mode}{ the method \code{\link{PlotGroups}} takes to condensate expression information when more than one gene is present in the data. Possible values are \code{"representative"} and \code{"median"} }  
  \item{color.mode}{ color scale for plotting profiles. Can be either \code{"rainblow"} or \code{"gray"}}
 \item{ylim}{range of the y axis to be used by \code{\link{PlotProfiles}} and \code{\link{PlotGroups}} }    
 \item{item}{Name of the analysed items to show}
 \item{legend}{logical indicating whether legend must be added when plotting profiles  }
 \item{cex.legend}{Expansion factor for legend}
 \item{lty.legend}{To add a coloured line in the legend}
 \item{\dots}{other graphical function argument }
}
\details{
    Data can be provided either as a single data matrix of expression values, or a \code{\link{get.siggenes}} object. In the later case
    the other argument of the fuction can be taken directly from \code{data}.
     
    Data clustering can be done on the basis of either the original expression values, the regression coefficients, 
    or the t.scores. In case \code{data} is a \code{get.siggenes} object, this is given by providing the element names of the list
    \code{c("sig.profiles","coefficients","t.score")} of their list position (1,2 or 3).
}
\value{
	Experiment wide gene profiles and by group profiles plots are generated for each data cluster in the graphical device.
  \item{cut}{vector indicating gene partioning into clusters}
  \item{c.algo.used}{clustering algorith used for data partioning}
  \item{groups}{groups matrix used for plotting functions}
}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments. 
Bioinformatics 22, 1096-1102
}
\author{ Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}
\seealso{ \code{\link{PlotProfiles}}, \code{\link{PlotGroups}}}
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

## Create 270 flat profiles
flat <- tc.GENE(n = 270, r = 3)
## Create 10 genes with profile differences between Ctl and Tr1 groups
twodiff <- tc.GENE (n = 10, r = 3, b2 = 0.5, c2 = 1.3)
## Create 10 genes with profile differences between Ctl, Tr2, and Tr3 groups
threediff <- tc.GENE(n = 10, r = 3, b3 = 0.8, c3 = -1, a4 = -0.1, b4 = -0.8, c4 = -1.2)
## Create 10 genes with profile differences between Ctl and Tr2 and different variance
vardiff <- tc.GENE(n = 10, r = 3, a3 = 0.7, b3 = 1, c3 = 1.2, var32 = 0.03, var33 = 0.03)
## Create dataset
tc.DATA <- rbind(flat, twodiff, threediff, vardiff)
rownames(tc.DATA) <- paste("feature", c(1:300), sep = "")
colnames(tc.DATA) <- paste("Array", c(1:36), sep = "")
tc.DATA [sample(c(1:(300*36)), 300)] <- NA  # introduce missing values

#### CREATE EXPERIMENTAL DESIGN
Time <- rep(c(rep(c(1:3), each = 3)), 4)
Replicates <- rep(c(1:12), each = 3)
Control <- c(rep(1, 9), rep(0, 27))
Treat1 <- c(rep(0, 9), rep(1, 9), rep(0, 18))
Treat2 <- c(rep(0, 18), rep(1, 9), rep(0,9))
Treat3 <- c(rep(0, 27), rep(1, 9))
edesign <- cbind(Time, Replicates, Control, Treat1, Treat2, Treat3)
rownames(edesign) <- paste("Array", c(1:36), sep = "")

see.genes(tc.DATA, edesign = edesign, k = 4)

# This will show the regression fit curve
dise <- make.design.matrix(edesign)
see.genes(tc.DATA, edesign = edesign, k = 4, show.fit = TRUE, 
          dis = dise$dis, groups.vector = dise$groups.vector, distance = "euclidean")

}
\keyword{ manip }
\keyword{ aplot }
