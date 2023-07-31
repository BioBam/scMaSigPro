\name{stepfor}
\alias{stepfor}

\title{ Fitting a linear model by forward-stepwise regression}
\description{
  \code{stepfor} fits a linear regression model applying forward-stepwise strategy.
}
\usage{
stepfor(y = y, d = d, alfa = 0.05, family = gaussian(), epsilon=0.00001  )
}

\arguments{
  \item{y}{ dependent variable }
  \item{d}{ data frame containing by columns the set of variables that could be in the selected model }
  \item{alfa}{ significance level to decide if a variable stays or not in the model}
  \item{family}{the distribution function to be used in the glm model}
  \item{epsilon}{argument to pass to \code{glm.control}, convergence tolerance in the iterative process to estimate de glm model}
}

\details{
The strategy begins analysing all the possible models with only one of the variables included in \code{d}. 
The most statistically significant variable (with the lowest p-value) is included in the model and 
then it is considered to introduce in the model another variable analysing all the possible models 
with two variables (the selected variable in the previous step plus a new variable). Again the most 
statistically significant variable (with lowest p-value) is included in the model. The process is 
repeated till there are no more statistically significant variables to include.       
}
\value{
     \code{stepfor} returns an object of the class \code{\link{lm}}, where the model uses 
     \code{y} as dependent variable and all the selected variables from \code{d} as independent variables. 
             
     The function \code{\link{summary}} are used to obtain a summary and analysis of variance table of the results. 
     The generic accessor functions \code{\link{coefficients}}, \code{\link{effects}},
     \code{\link{fitted.values}} and \code{\link{residuals}} extract various useful features of the value returned by \code{\link{lm}}.
}
\references{Conesa, A., Nueda M.J., Alberto Ferrer, A., Talon, T. 2005.
maSigPro: a Method to Identify Significant Differential Expression Profiles in Time-Course Microarray Experiments.
 }
\author{ Ana Conesa, aconesa@cipf.es; Maria Jose Nueda, mj.nueda@ua.es}
\seealso{ \code{\link{lm}}, \code{\link{step}}, \code{\link{stepback}}, \code{\link{two.ways.stepback}}, \code{\link{two.ways.stepfor}}}
\examples{

## create design matrix
Time <- rep(c(rep(c(1:3), each = 3)), 4)
Replicates <- rep(c(1:12), each = 3)
Control <- c(rep(1, 9), rep(0, 27))
Treat1 <- c(rep(0, 9), rep(1, 9), rep(0, 18))
Treat2 <- c(rep(0, 18), rep(1, 9), rep(0,9))
Treat3 <- c(rep(0, 27), rep(1, 9))
edesign <- cbind(Time, Replicates, Control, Treat1, Treat2, Treat3)
rownames(edesign) <- paste("Array", c(1:36), sep = "")
dise <- make.design.matrix(edesign)
dis <- as.data.frame(dise$dis)


## expression vector
y <- c(0.082, 0.021, 0.010, 0.113, 0.013, 0.077, 0.068, 0.042, -0.056, -0.232, -0.014, -0.040, 
-0.055, 0.150, -0.027, 0.064, -0.108, -0.220, 0.275, -0.130, 0.130, 1.018, 1.005, 0.931,
 -1.009, -1.101, -1.014, -0.045, -0.110, -0.128, -0.643, -0.785, -1.077, -1.187, -1.249, -1.463)

s.fit <- stepfor(y = y, d = dis)
summary(s.fit)

}
\keyword{ regression}