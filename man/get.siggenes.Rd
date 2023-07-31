\name{get.siggenes}
\alias{get.siggenes}

\title{ Extract significant genes for sets of variables in  time series gene expression experiments }
\description{
  This function creates lists of significant genes for a set of variables whose significance value has been computed with the \code{T.fit} function.
}
\usage{
get.siggenes(tstep, rsq = 0.7, add.IDs = FALSE, IDs = NULL, matchID.col = 1,
             only.names = FALSE, vars = c("all", "each", "groups"),
	     significant.intercept = "dummy",

             groups.vector = NULL, trat.repl.spots = "none",
             index = IDs[, (matchID.col + 1)], match = IDs[, matchID.col],
	     r = 0.7)

}

\arguments{
  \item{tstep}{ a \code{T.fit} object }
  \item{rsq}{ cut-off level at the R-squared value for the stepwise regression fit. Only genes with R-squared more than rsq are selected}
  \item{add.IDs}{ logical indicating whether to include additional gene id's in the result}
  \item{IDs}{ matrix contaning additional gene id information (required when \code{add.IDs} is \code{TRUE}) }
  \item{matchID.col}{ number of matching column in matrix IDs for adding genes ids }
  \item{only.names}{ logical. If \code{TRUE}, expression values are ommited in the results }
  \item{vars}{ variables for which to extract significant genes (see details) }
  \item{significant.intercept}{ experimental groups for which significant intercept coefficients are considered (see details) }
  \item{groups.vector}{ required when \code{vars} is \code{"groups"}.}
  \item{trat.repl.spots}{ treatment given to replicate spots. Possible values are \code{"none"} and \code{"average"}}
  \item{index}{ argument of the \code{\link{average.rows}} function to use when \code{ trat.repl.spots} is \code{"average"} }
  \item{match}{ argument of the \code{\link{average.rows}} function to use when \code{trat.repl.spots} is \code{"average"} }
  \item{r}{ minimun pearson correlation coefficient for replicated spots profiles to be averaged }
}
\details{
    There are 3 possible values for the vars argument:

    \code{"all"}: generates one single matrix or gene list with all
 significant genes.


   \code{"each"}: generates as many significant genes extractions as
 variables in the general regression model. Each extraction contains the
 significant genes for that variable.


   \code{"groups"}: generates a significant genes extraction for each
 experimental group.


  The difference between \code{"each"} and \code{"groups"} is that in the
     first case the variables of the same group (e.g.  \code{"TreatmentA"} and
     \code{"time*TreatmentA" }) will be extracted separately and in the second
     case jointly.


    When \code{add.IDs} is \code{TRUE}, a matrix of gene ids must be provided
    as argument of IDs, the \code{matchID.col}
    column of which having same levels as in the row names of
    \code{sig.profiles}. The option \code{only.names} is \code{TRUE} will
    generate a vector of significant genes or a matrix when \code{add.IDs} is
    set also to \code{TRUE}.

    When \code{trat.repl.spots} is \code{"average"}, \code{match} and \code{index} vectors are required for the \code{\link{average.rows}} function.
    In gene expression data context, the \code{index} vector would contain geneIDs and indicate which spots
    are replicates. The \code{match} vector is used to match these genesIDs to rows in the significant genes
    matrix, and must have the same levels as the row names of \code{sig.profiles}.

    The argument \code{significant.intercept} modulates the treatment for intercept coefficients to apply for selecting significant genes
    when \code{vars} equals \code{"groups"}. There are three possible values: \code{"none"}, no significant intercept (differences) are
    considered for significant gene selection, \code{"dummy"}, includes genes with significant intercept differences between control and experimental
    groups, and \code{"all"} when both significant intercept coefficient for the control group and significant intercept
    differences are considered for selecting significant genes.

    \code{add.IDs} = TRUE and \code{trat.repl.spots} = \code{"average"} are not compatible argumet values.
    \code{add.IDs} = TRUE and \code{only.names} = \code{TRUE} are  compatible argumet values.
}

\value{
  \item{summary}{a vector or matrix listing significant genes for the variables given by the function parameters}
  \item{sig.genes}{a list with detailed information on the significant genes
    found for the variables given by the function parameters. Each element of
    the list is also a list containing:
    \describe{
      \item{\code{sig.profiles}:}{ expression values of significant genes}
      \item{\code{coefficients}:}{ regression coefficients of the adjusted
      models}
      \item{\code{groups.coeffs}:}{ regression coefficients of the impiclit
         models of each experimental group}
      \item{\code{sig.pvalues}:}{ p-values of the regression coefficients for
         significant genes}
      \item{\code{g}:}{ number of genes}
      \item{\code{...}:}{ arguments passed by previous functions}
    }
  }
}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2006.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
Bioinformatics 22, 1096-1102
}
\author{Ana Conesa and Maria Jose Nueda, \email{mj.nueda@ua.es}}
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

tc.p <- p.vector(tc.DATA, design = make.design.matrix(edesign), Q = 0.01)
tc.tstep <- T.fit(data = tc.p , alfa = 0.05)

## This will obtain sigificant genes per experimental group
## which have a regression model Rsquared > 0.9
tc.sigs <- get.siggenes (tc.tstep, rsq = 0.9, vars = "groups")

## This will obtain all sigificant genes regardless the Rsquared value.
## Replicated genes are averaged.
IDs <- rbind(paste("feature", c(1:300), sep = ""),
       rep(paste("gene", c(1:150), sep = ""), each = 2))
tc.sigs.ALL <- get.siggenes (tc.tstep, rsq = 0, vars = "all", IDs = IDs)
tc.sigs.groups <- get.siggenes (tc.tstep, rsq = 0, vars = "groups", significant.intercept="dummy")

}
\keyword{ manip }
